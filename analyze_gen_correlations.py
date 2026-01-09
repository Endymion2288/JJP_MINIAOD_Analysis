#!/usr/bin/env python3
"""
Analyze MiniAOD prunedGenParticles to study correlations between
Jpsi1, Jpsi2, and phi particles.

Selection modes:
- SPS:   Jpsi1, Jpsi2, Phi come from the same parton-level interaction
         Phi comes from gluon decay in this interaction
- DPS_1: Jpsi1 and Phi come from the same parton-level interaction (Phi from gluon)
         Jpsi2 comes from a different parton-level interaction
- DPS_2: Jpsi1 and Jpsi2 come from the same parton-level interaction
         Phi comes from a different parton-level interaction (from gluon decay)
- TPS:   Jpsi1, Jpsi2, Phi come from different parton-level interactions
         Phi comes from gluon decay

Output: ROOT file with histograms and 2D plots
"""

import ROOT
import glob
import os
import sys
from array import array
import multiprocessing as mp
from functools import partial
import tempfile
import argparse

# PDG IDs
PDG_JPSI = 443
PDG_PHI = 333
PDG_GLUON = 21

def delta_phi(phi1, phi2):
    """Calculate delta phi wrapped to [-pi, pi]"""
    dphi = phi1 - phi2
    while dphi > ROOT.TMath.Pi():
        dphi -= 2.0 * ROOT.TMath.Pi()
    while dphi < -ROOT.TMath.Pi():
        dphi += 2.0 * ROOT.TMath.Pi()
    return dphi

def rapidity(particle):
    """Calculate rapidity y = 0.5 * ln((E + pz) / (E - pz))"""
    E = particle.energy()
    pz = particle.pz()
    if E - pz <= 0 or E + pz <= 0:
        # Fallback to pseudorapidity if rapidity is undefined
        return particle.eta()
    return 0.5 * ROOT.TMath.Log((E + pz) / (E - pz))

def get_mother_chain(particle, cache=None):
    """Get the chain of mothers for a particle; cache to avoid re-traversal."""
    if cache is not None:
        key = id(particle)
        if key in cache:
            return cache[key]

    mothers = []
    p = particle
    while p.numberOfMothers() > 0:
        mother = p.mother(0)
        mothers.append(mother)
        p = mother

    if cache is not None:
        cache[id(particle)] = mothers
    return mothers


def get_ancestor_set(particle, chain_cache, ancestor_cache):
    """Get the set of ancestors for a particle as tuples for comparison."""
    key = id(particle)
    if key in ancestor_cache:
        return ancestor_cache[key]
    moms = get_mother_chain(particle, chain_cache)
    aset = set((m.pdgId(), round(m.pt(), 2), round(m.eta(), 3), round(m.phi(), 3)) for m in moms)
    ancestor_cache[key] = aset
    return aset


def find_common_ancestor(particle1, particle2, chain_cache, ancestor_cache):
    """Find if two particles share a common ancestor using cached chains."""
    a1 = get_ancestor_set(particle1, chain_cache, ancestor_cache)
    a2 = get_ancestor_set(particle2, chain_cache, ancestor_cache)
    common = a1.intersection(a2)
    return len(common) > 0, common


def particles_share_ancestor(p1, p2, p3, chain_cache, ancestor_cache):
    """Check if three particles all share a common ancestor."""
    a1 = get_ancestor_set(p1, chain_cache, ancestor_cache)
    a2 = get_ancestor_set(p2, chain_cache, ancestor_cache)
    a3 = get_ancestor_set(p3, chain_cache, ancestor_cache)
    common = a1.intersection(a2).intersection(a3)
    return len(common) > 0, common


def phi_from_gluon(phi_particle, chain_cache):
    """Check if phi comes from gluon decay and return the gluon (cached)."""
    mothers = get_mother_chain(phi_particle, chain_cache)
    for m in mothers:
        if abs(m.pdgId()) == PDG_GLUON:
            return True, m
    return False, None


def phi_shares_gluon_with(phi_particle, other_particle, chain_cache, ancestor_cache):
    """
    Check if phi comes from a gluon that is in the same parton chain as other_particle.
    Returns True if the gluon ancestor of phi is also an ancestor of other_particle.
    """
    is_from_gluon, gluon = phi_from_gluon(phi_particle, chain_cache)
    if not is_from_gluon:
        return False
    
    # Check if the gluon (or its ancestors) are also ancestors of other_particle
    other_ancestors = get_ancestor_set(other_particle, chain_cache, ancestor_cache)
    
    # Check if gluon is in the ancestor chain
    gluon_tuple = (gluon.pdgId(), round(gluon.pt(), 2), round(gluon.eta(), 3), round(gluon.phi(), 3))
    if gluon_tuple in other_ancestors:
        return True
    
    # Also check if they share common ancestors through the gluon
    gluon_ancestors = get_ancestor_set(gluon, chain_cache, ancestor_cache)
    common = gluon_ancestors.intersection(other_ancestors)
    return len(common) > 0


def create_histograms(prefix=""):
    """Create all histograms with optional prefix for naming."""
    histograms = {}
    
    # 1D histograms for delta y (rapidity difference)
    histograms['h_dy_jpsi1_jpsi2'] = ROOT.TH1F(f"{prefix}h_dy_jpsi1_jpsi2", 
        "Delta y (J/#psi_{1} - J/#psi_{2});|#Delta y|;Events", 50, 0, 5)
    histograms['h_dy_jpsi1_phi'] = ROOT.TH1F(f"{prefix}h_dy_jpsi1_phi", 
        "Delta y (J/#psi_{1} - #phi);|#Delta y|;Events", 50, 0, 5)
    histograms['h_dy_jpsi2_phi'] = ROOT.TH1F(f"{prefix}h_dy_jpsi2_phi", 
        "Delta y (J/#psi_{2} - #phi);|#Delta y|;Events", 50, 0, 5)
    
    # 1D histograms for delta phi
    histograms['h_dphi_jpsi1_jpsi2'] = ROOT.TH1F(f"{prefix}h_dphi_jpsi1_jpsi2", 
        "Delta #phi (J/#psi_{1} - J/#psi_{2});|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    histograms['h_dphi_jpsi1_phi'] = ROOT.TH1F(f"{prefix}h_dphi_jpsi1_phi", 
        "Delta #phi (J/#psi_{1} - #phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    histograms['h_dphi_jpsi2_phi'] = ROOT.TH1F(f"{prefix}h_dphi_jpsi2_phi", 
        "Delta #phi (J/#psi_{2} - #phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    
    # 2D histograms: delta y vs delta phi
    histograms['h2_dy_dphi_jpsi1_jpsi2'] = ROOT.TH2F(f"{prefix}h2_dy_dphi_jpsi1_jpsi2", 
        "J/#psi_{1} - J/#psi_{2}: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
        50, 0, 5, 50, 0, ROOT.TMath.Pi())
    histograms['h2_dy_dphi_jpsi1_phi'] = ROOT.TH2F(f"{prefix}h2_dy_dphi_jpsi1_phi", 
        "J/#psi_{1} - #phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
        50, 0, 5, 50, 0, ROOT.TMath.Pi())
    histograms['h2_dy_dphi_jpsi2_phi'] = ROOT.TH2F(f"{prefix}h2_dy_dphi_jpsi2_phi", 
        "J/#psi_{2} - #phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
        50, 0, 5, 50, 0, ROOT.TMath.Pi())
    
    # Kinematic distributions
    histograms['h_jpsi1_pt'] = ROOT.TH1F(f"{prefix}h_jpsi1_pt", 
        "J/#psi_{1} p_{T};p_{T} [GeV];Events", 100, 0, 50)
    histograms['h_jpsi2_pt'] = ROOT.TH1F(f"{prefix}h_jpsi2_pt", 
        "J/#psi_{2} p_{T};p_{T} [GeV];Events", 100, 0, 50)
    histograms['h_phi_pt'] = ROOT.TH1F(f"{prefix}h_phi_pt", 
        "#phi p_{T};p_{T} [GeV];Events", 100, 0, 50)
    
    histograms['h_jpsi1_eta'] = ROOT.TH1F(f"{prefix}h_jpsi1_eta", 
        "J/#psi_{1} #eta;#eta;Events", 60, -3, 3)
    histograms['h_jpsi2_eta'] = ROOT.TH1F(f"{prefix}h_jpsi2_eta", 
        "J/#psi_{2} #eta;#eta;Events", 60, -3, 3)
    histograms['h_phi_eta'] = ROOT.TH1F(f"{prefix}h_phi_eta", 
        "#phi #eta;#eta;Events", 60, -3, 3)
    
    histograms['h_jpsi1_y'] = ROOT.TH1F(f"{prefix}h_jpsi1_y", 
        "J/#psi_{1} y;y;Events", 60, -3, 3)
    histograms['h_jpsi2_y'] = ROOT.TH1F(f"{prefix}h_jpsi2_y", 
        "J/#psi_{2} y;y;Events", 60, -3, 3)
    histograms['h_phi_y'] = ROOT.TH1F(f"{prefix}h_phi_y", 
        "#phi y;y;Events", 60, -3, 3)
    
    # Invariant mass distributions
    histograms['h_mass_jpsi1_jpsi2'] = ROOT.TH1F(f"{prefix}h_mass_jpsi1_jpsi2",
        "M(J/#psi_{1} + J/#psi_{2});M [GeV];Events", 100, 6, 30)
    histograms['h_mass_jpsi1_phi'] = ROOT.TH1F(f"{prefix}h_mass_jpsi1_phi",
        "M(J/#psi_{1} + #phi);M [GeV];Events", 100, 3, 20)
    histograms['h_mass_jpsi2_phi'] = ROOT.TH1F(f"{prefix}h_mass_jpsi2_phi",
        "M(J/#psi_{2} + #phi);M [GeV];Events", 100, 3, 20)
    histograms['h_mass_all'] = ROOT.TH1F(f"{prefix}h_mass_all",
        "M(J/#psi_{1} + J/#psi_{2} + #phi);M [GeV];Events", 100, 7, 40)
    
    return histograms


def fill_histograms(histograms, jpsi1, jpsi2, phi):
    """Fill all histograms with the selected particles."""
    # Calculate rapidities
    y_jpsi1 = rapidity(jpsi1)
    y_jpsi2 = rapidity(jpsi2)
    y_phi = rapidity(phi)
    
    # Delta y
    dy_jpsi1_jpsi2 = abs(y_jpsi1 - y_jpsi2)
    dy_jpsi1_phi = abs(y_jpsi1 - y_phi)
    dy_jpsi2_phi = abs(y_jpsi2 - y_phi)
    
    # Delta phi
    dphi_jpsi1_jpsi2 = abs(delta_phi(jpsi1.phi(), jpsi2.phi()))
    dphi_jpsi1_phi = abs(delta_phi(jpsi1.phi(), phi.phi()))
    dphi_jpsi2_phi = abs(delta_phi(jpsi2.phi(), phi.phi()))
    
    # Fill 1D histograms
    histograms['h_dy_jpsi1_jpsi2'].Fill(dy_jpsi1_jpsi2)
    histograms['h_dy_jpsi1_phi'].Fill(dy_jpsi1_phi)
    histograms['h_dy_jpsi2_phi'].Fill(dy_jpsi2_phi)
    
    histograms['h_dphi_jpsi1_jpsi2'].Fill(dphi_jpsi1_jpsi2)
    histograms['h_dphi_jpsi1_phi'].Fill(dphi_jpsi1_phi)
    histograms['h_dphi_jpsi2_phi'].Fill(dphi_jpsi2_phi)
    
    # Fill 2D histograms
    histograms['h2_dy_dphi_jpsi1_jpsi2'].Fill(dy_jpsi1_jpsi2, dphi_jpsi1_jpsi2)
    histograms['h2_dy_dphi_jpsi1_phi'].Fill(dy_jpsi1_phi, dphi_jpsi1_phi)
    histograms['h2_dy_dphi_jpsi2_phi'].Fill(dy_jpsi2_phi, dphi_jpsi2_phi)
    
    # Fill kinematic histograms
    histograms['h_jpsi1_pt'].Fill(jpsi1.pt())
    histograms['h_jpsi2_pt'].Fill(jpsi2.pt())
    histograms['h_phi_pt'].Fill(phi.pt())
    
    histograms['h_jpsi1_eta'].Fill(jpsi1.eta())
    histograms['h_jpsi2_eta'].Fill(jpsi2.eta())
    histograms['h_phi_eta'].Fill(phi.eta())
    
    histograms['h_jpsi1_y'].Fill(y_jpsi1)
    histograms['h_jpsi2_y'].Fill(y_jpsi2)
    histograms['h_phi_y'].Fill(y_phi)
    
    # Calculate invariant masses using TLorentzVector
    v_jpsi1 = ROOT.TLorentzVector()
    v_jpsi2 = ROOT.TLorentzVector()
    v_phi = ROOT.TLorentzVector()
    
    v_jpsi1.SetPtEtaPhiM(jpsi1.pt(), jpsi1.eta(), jpsi1.phi(), jpsi1.mass())
    v_jpsi2.SetPtEtaPhiM(jpsi2.pt(), jpsi2.eta(), jpsi2.phi(), jpsi2.mass())
    v_phi.SetPtEtaPhiM(phi.pt(), phi.eta(), phi.phi(), phi.mass())
    
    histograms['h_mass_jpsi1_jpsi2'].Fill((v_jpsi1 + v_jpsi2).M())
    histograms['h_mass_jpsi1_phi'].Fill((v_jpsi1 + v_phi).M())
    histograms['h_mass_jpsi2_phi'].Fill((v_jpsi2 + v_phi).M())
    histograms['h_mass_all'].Fill((v_jpsi1 + v_jpsi2 + v_phi).M())


def select_sps(jpsis, phis, chain_cache, ancestor_cache):
    """
    SPS selection: Jpsi1, Jpsi2, Phi come from the same parton-level interaction
    Phi comes from gluon decay in this interaction.
    Returns: (jpsi1, jpsi2, phi) or (None, None, None)
    """
    # Need at least 2 J/psi
    if len(jpsis) < 2:
        return None, None, None
    
    # Try all combinations of two J/psi
    for i, jpsi1 in enumerate(jpsis):
        for jpsi2 in jpsis[i+1:]:
            # Check if Jpsi1 and Jpsi2 share common ancestor
            has_common_jj, _ = find_common_ancestor(jpsi1, jpsi2, chain_cache, ancestor_cache)
            if not has_common_jj:
                continue
            
            for phi_cand in phis:
                # Check if phi comes from gluon
                is_from_gluon, _ = phi_from_gluon(phi_cand, chain_cache)
                if not is_from_gluon:
                    continue
                
                # Check if all three share common ancestor
                all_common, _ = particles_share_ancestor(jpsi1, jpsi2, phi_cand, chain_cache, ancestor_cache)
                if all_common:
                    return jpsi1, jpsi2, phi_cand
    
    return None, None, None


def select_dps1(jpsis, phis, chain_cache, ancestor_cache):
    """
    DPS_1 selection: Jpsi1 and Phi come from the same parton-level interaction (Phi from gluon)
                     Jpsi2 comes from a different parton-level interaction
    Returns: (jpsi1, jpsi2, phi) or (None, None, None)
    """
    # Need at least 2 J/psi
    if len(jpsis) < 2:
        return None, None, None
    
    for jpsi1 in jpsis:
        for phi_cand in phis:
            # Check if phi comes from gluon
            is_from_gluon, _ = phi_from_gluon(phi_cand, chain_cache)
            if not is_from_gluon:
                continue
            
            # Check if Jpsi1 and Phi share common ancestor (same parton interaction)
            has_common_j1p, _ = find_common_ancestor(jpsi1, phi_cand, chain_cache, ancestor_cache)
            if not has_common_j1p:
                continue
            
            # Find Jpsi2 from a different parton interaction
            for jpsi2 in jpsis:
                if jpsi2 is jpsi1:
                    continue
                    
                # Jpsi2 should NOT share common ancestor with Jpsi1
                has_common_jj, _ = find_common_ancestor(jpsi1, jpsi2, chain_cache, ancestor_cache)
                # Also check Jpsi2 doesn't share ancestor with Phi
                has_common_j2p, _ = find_common_ancestor(jpsi2, phi_cand, chain_cache, ancestor_cache)
                
                if not has_common_jj and not has_common_j2p:
                    return jpsi1, jpsi2, phi_cand
    
    return None, None, None


def select_dps2(jpsis, phis, chain_cache, ancestor_cache):
    """
    DPS_2 selection: Jpsi1 and Jpsi2 come from the same parton-level interaction
                     Phi comes from a different parton-level interaction (from gluon decay)
    Returns: (jpsi1, jpsi2, phi) or (None, None, None)
    """
    # Need at least 2 J/psi
    if len(jpsis) < 2:
        return None, None, None
    
    for i, jpsi1 in enumerate(jpsis):
        for jpsi2 in jpsis[i+1:]:
            # Check if Jpsi1 and Jpsi2 share common ancestor (same parton interaction)
            has_common_jj, _ = find_common_ancestor(jpsi1, jpsi2, chain_cache, ancestor_cache)
            if not has_common_jj:
                continue
            
            # Find Phi from a different parton interaction (and from gluon)
            for phi_cand in phis:
                # Check if phi comes from gluon
                is_from_gluon, _ = phi_from_gluon(phi_cand, chain_cache)
                if not is_from_gluon:
                    continue
                
                # Phi should NOT share common ancestor with Jpsi1 or Jpsi2
                has_common_j1p, _ = find_common_ancestor(jpsi1, phi_cand, chain_cache, ancestor_cache)
                has_common_j2p, _ = find_common_ancestor(jpsi2, phi_cand, chain_cache, ancestor_cache)
                
                if not has_common_j1p and not has_common_j2p:
                    return jpsi1, jpsi2, phi_cand
    
    return None, None, None


def select_tps(jpsis, phis, chain_cache, ancestor_cache):
    """
    TPS selection: only require Phi comes from gluon decay.
    Choose highest-pT Jpsi1, Jpsi2, and highest-pT Phi from gluon.
    Returns: (jpsi1, jpsi2, phi) or (None, None, None)
    """
    # Need at least 2 J/psi
    if len(jpsis) < 2:
        return None, None, None

    # jpsis and phis are already sorted by pT (descending) upstream
    jpsi1 = jpsis[0]
    jpsi2 = jpsis[1]

    for phi_cand in phis:
        is_from_gluon, _ = phi_from_gluon(phi_cand, chain_cache)
        if is_from_gluon:
            return jpsi1, jpsi2, phi_cand

    return None, None, None


def get_selection_function(mode):
    """Return the appropriate selection function based on mode."""
    selection_functions = {
        'SPS': select_sps,
        'DPS_1': select_dps1,
        'DPS_2': select_dps2,
        'TPS': select_tps,
    }
    return selection_functions.get(mode.upper())


def process_file_batch(file_batch, batch_id, output_dir, mode, max_events_per_batch=-1):
    """
    Process a batch of files in a worker process.
    Returns the path to temporary output ROOT file.
    """
    # Each worker must initialize ROOT/FWLite independently
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    
    from DataFormats.FWLite import Events, Handle
    
    # Create temporary output file for this batch
    temp_output = os.path.join(output_dir, f"temp_batch_{batch_id}.root")
    fout = ROOT.TFile(temp_output, "RECREATE")
    
    # Create histograms for this batch
    histograms = create_histograms()
    
    # Get selection function
    select_func = get_selection_function(mode)
    
    handle_genParticles = Handle("std::vector<reco::GenParticle>")
    label_genParticles = ("prunedGenParticles", "", "PAT")
    
    n_total = 0
    n_selected = 0
    
    for fname in file_batch:
        try:
            events = Events(fname)
        except Exception as e:
            print(f"  [Batch {batch_id}] Error opening {fname}: {e}")
            continue
        
        for i, event in enumerate(events):
            if max_events_per_batch > 0 and n_total >= max_events_per_batch:
                break
                
            n_total += 1
            
            event.getByLabel(label_genParticles, handle_genParticles)
            if not handle_genParticles.isValid():
                continue
            
            genParticles = handle_genParticles.product()
            
            jpsis = []
            phis = []
            
            for gp in genParticles:
                pdg = abs(gp.pdgId())
                
                if pdg == PDG_JPSI:
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        jpsis.append(gp)
                    elif gp.status() == 2:
                        jpsis.append(gp)
                
                if pdg == PDG_PHI:
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        phis.append(gp)
                    elif gp.status() == 2:
                        phis.append(gp)
            
            # Need at least 2 J/psi and 1 phi
            if len(jpsis) < 2 or len(phis) < 1:
                continue

            # Sort by pT (highest first)
            jpsis.sort(key=lambda x: x.pt(), reverse=True)
            phis.sort(key=lambda x: x.pt(), reverse=True)

            # Per-event caches
            mother_chain_cache = {}
            ancestor_cache = {}

            # Apply selection
            jpsi1, jpsi2, phi = select_func(jpsis, phis, mother_chain_cache, ancestor_cache)
            
            if jpsi1 is None or jpsi2 is None or phi is None:
                continue
            
            n_selected += 1
            
            # Fill histograms
            fill_histograms(histograms, jpsi1, jpsi2, phi)
        
        if max_events_per_batch > 0 and n_total >= max_events_per_batch:
            break
    
    fout.Write()
    fout.Close()
    
    print(f"  [Batch {batch_id}] Processed {len(file_batch)} files, {n_total} events, {n_selected} selected")
    
    return temp_output, n_total, n_selected


def merge_histograms(temp_files, output_file):
    """Merge histograms from temporary files into final output."""
    print(f"\nMerging {len(temp_files)} temporary files...")
    
    # Use ROOT's hadd-like functionality
    merger = ROOT.TFileMerger(False)
    merger.SetFastMethod(True)
    merger.OutputFile(output_file)
    
    for tf in temp_files:
        if os.path.exists(tf):
            merger.AddFile(tf)
    
    success = merger.Merge()
    
    # Clean up temporary files
    for tf in temp_files:
        if os.path.exists(tf):
            os.remove(tf)
    
    return success


def analyze_miniaod_files(input_files, output_file, mode, max_events=-1, n_workers=1):
    """
    Main analysis function to process MiniAOD files.
    Supports parallel processing with n_workers > 1.
    """
    print(f"\n=== Running {mode} Selection ===")
    
    if n_workers > 1:
        # Parallel processing mode
        print(f"Running in parallel mode with {n_workers} workers")
        
        # Create temporary directory for intermediate files
        temp_dir = tempfile.mkdtemp(prefix="jjp_gen_corr_")
        print(f"Temporary directory: {temp_dir}")
        
        # Split files into batches
        n_files = len(input_files)
        files_per_worker = (n_files + n_workers - 1) // n_workers
        file_batches = []
        for i in range(0, n_files, files_per_worker):
            file_batches.append(input_files[i:i+files_per_worker])
        
        print(f"Split {n_files} files into {len(file_batches)} batches")
        
        # Calculate max events per batch if total limit is set
        max_events_per_batch = -1
        if max_events > 0:
            max_events_per_batch = (max_events + len(file_batches) - 1) // len(file_batches)
        
        # Process batches in parallel using multiprocessing Pool
        ctx = mp.get_context('spawn')
        
        results = []
        with ctx.Pool(processes=n_workers) as pool:
            tasks = []
            for batch_id, batch in enumerate(file_batches):
                task = pool.apply_async(
                    process_file_batch,
                    args=(batch, batch_id, temp_dir, mode, max_events_per_batch)
                )
                tasks.append(task)
            
            # Collect results
            for task in tasks:
                try:
                    result = task.get(timeout=7200)  # 2 hour timeout per batch
                    results.append(result)
                except Exception as e:
                    print(f"Error in worker: {e}")
        
        # Gather statistics
        temp_files = [r[0] for r in results]
        total_events = sum(r[1] for r in results)
        selected_events = sum(r[2] for r in results)
        
        # Merge all temporary files
        merge_histograms(temp_files, output_file)
        
        # Clean up temp directory
        try:
            os.rmdir(temp_dir)
        except:
            pass
        
        print(f"\n=== Analysis Summary ({mode}) ===")
        print(f"Total events processed: {total_events}")
        print(f"Events with valid J/psi1, J/psi2, Phi selection: {selected_events}")
        print(f"Selection efficiency: {100.0*selected_events/total_events:.2f}%" if total_events > 0 else "N/A")
        print(f"\nOutput saved to: {output_file}")
        
    else:
        # Sequential processing mode
        analyze_miniaod_files_sequential(input_files, output_file, mode, max_events)


def analyze_miniaod_files_sequential(input_files, output_file, mode, max_events=-1):
    """
    Main analysis function to process MiniAOD files (sequential mode)
    """
    # Load FWLite
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    
    from DataFormats.FWLite import Events, Handle
    
    # Create output file and histograms
    fout = ROOT.TFile(output_file, "RECREATE")
    histograms = create_histograms()
    
    # Get selection function
    select_func = get_selection_function(mode)
    if select_func is None:
        print(f"Error: Unknown selection mode '{mode}'")
        sys.exit(1)
    
    # Handle for prunedGenParticles
    handle_genParticles = Handle("std::vector<reco::GenParticle>")
    label_genParticles = ("prunedGenParticles", "", "PAT")
    
    n_total = 0
    n_selected = 0
    
    for fname in input_files:
        print(f"Processing: {fname}")
        
        try:
            events = Events(fname)
        except Exception as e:
            print(f"  Error opening file: {e}")
            continue
        
        for i, event in enumerate(events):
            if max_events > 0 and n_total >= max_events:
                break
                
            n_total += 1
            if n_total % 1000 == 0:
                print(f"  Processed {n_total} events, selected {n_selected}...")
            
            # Get gen particles
            event.getByLabel(label_genParticles, handle_genParticles)
            if not handle_genParticles.isValid():
                continue
            
            genParticles = handle_genParticles.product()
            
            # Find J/psi and phi particles
            jpsis = []
            phis = []
            
            for gp in genParticles:
                pdg = abs(gp.pdgId())
                
                if pdg == PDG_JPSI:
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        jpsis.append(gp)
                    elif gp.status() == 2:
                        jpsis.append(gp)
                
                if pdg == PDG_PHI:
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        phis.append(gp)
                    elif gp.status() == 2:
                        phis.append(gp)
            
            # Need at least 2 J/psi and 1 phi
            if len(jpsis) < 2 or len(phis) < 1:
                continue

            # Sort by pT (highest first)
            jpsis.sort(key=lambda x: x.pt(), reverse=True)
            phis.sort(key=lambda x: x.pt(), reverse=True)

            # Per-event caches to speed up mother/ancestor traversals
            mother_chain_cache = {}
            ancestor_cache = {}

            # Apply selection based on mode
            jpsi1, jpsi2, phi = select_func(jpsis, phis, mother_chain_cache, ancestor_cache)
            
            if jpsi1 is None or jpsi2 is None or phi is None:
                continue
            
            n_selected += 1
            
            # Fill histograms
            fill_histograms(histograms, jpsi1, jpsi2, phi)
        
        if max_events > 0 and n_total >= max_events:
            break
    
    print(f"\n=== Analysis Summary ({mode}) ===")
    print(f"Total events processed: {n_total}")
    print(f"Events with valid J/psi1, J/psi2, Phi selection: {n_selected}")
    print(f"Selection efficiency: {100.0*n_selected/n_total:.2f}%" if n_total > 0 else "N/A")
    
    # Write and close
    fout.Write()
    fout.Close()
    print(f"\nOutput saved to: {output_file}")


def get_input_files_from_xrootd(base_dir, xrootd_server="root://cceos.ihep.ac.cn"):
    """
    Get list of MINIAOD files from xrootd directory structure.
    Expected structure: base_dir/*/output_MINIAOD.root
    """
    import subprocess
    
    files = []
    
    # First, list subdirectories
    cmd = f"xrdfs {xrootd_server} ls {base_dir}"
    try:
        result = subprocess.run(cmd.split(), capture_output=True, text=True, timeout=60)
        subdirs = result.stdout.strip().split('\n')
        
        for subdir in subdirs:
            if not subdir:
                continue
            # Check for MINIAOD file
            miniaod_path = f"{subdir}/output_MINIAOD.root"
            full_path = f"{xrootd_server}/{miniaod_path}"
            files.append(full_path)
        
    except Exception as e:
        print(f"Error listing xrootd directory: {e}")
    
    return files


def main():
    parser = argparse.ArgumentParser(
        description='Analyze MiniAOD gen particles for J/psi + J/psi + Phi correlations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Selection modes:
  SPS   - All three particles from same parton interaction, Phi from gluon
  DPS_1 - J/psi1 + Phi from same interaction (Phi from gluon), J/psi2 separate
  DPS_2 - J/psi1 + J/psi2 from same interaction, Phi separate (from gluon)
  TPS   - All three from different interactions, Phi from gluon

Examples:
  # Process local files
  python3 analyze_gen_correlations.py -i /path/to/miniaod/ -o output.root -m DPS_1
  
  # Process from IHEP xrootd
  python3 analyze_gen_correlations.py --xrootd-base /eos/ihep/cms/store/user/xcheng/MC_Production/output/JJP_DPS1 -o dps1_output.root -m DPS_1
""")
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input-dir', 
                        help='Input directory containing MiniAOD files (local or EOS)')
    input_group.add_argument('--xrootd-base',
                        help='Base directory on xrootd (e.g., /eos/ihep/cms/store/user/...)')
    input_group.add_argument('--file-list',
                        help='Text file containing list of input files (one per line)')
    
    parser.add_argument('-o', '--output', 
                        default='gen_correlation_histograms.root',
                        help='Output ROOT file (default: gen_correlation_histograms.root)')
    parser.add_argument('-m', '--mode', 
                        choices=['SPS', 'DPS_1', 'DPS_2', 'TPS'],
                        default='DPS_1',
                        help='Selection mode (default: DPS_1)')
    parser.add_argument('-n', '--max-events', type=int, default=-1,
                        help='Maximum number of events to process (-1 for all)')
    parser.add_argument('--max-files', type=int, default=-1,
                        help='Maximum number of files to process (-1 for all)')
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help='Number of parallel workers (default: 1, sequential)')
    parser.add_argument('--xrootd-server', 
                        default='root://cceos.ihep.ac.cn',
                        help='xrootd server address (default: root://cceos.ihep.ac.cn)')
    
    args = parser.parse_args()
    
    # Get input files based on input method
    input_files = []
    
    if args.input_dir:
        # Local or EOS directory with *.root pattern
        input_pattern = os.path.join(args.input_dir, '*.root')
        input_files = sorted(glob.glob(input_pattern))
        
        # Also try subdirectory pattern
        if not input_files:
            input_pattern = os.path.join(args.input_dir, '*', 'output_MINIAOD.root')
            input_files = sorted(glob.glob(input_pattern))
    
    elif args.xrootd_base:
        # Get files from xrootd
        input_files = get_input_files_from_xrootd(args.xrootd_base, args.xrootd_server)
    
    elif args.file_list:
        # Read from file list
        with open(args.file_list, 'r') as f:
            input_files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    if not input_files:
        print(f"Error: No ROOT files found")
        sys.exit(1)
    
    print(f"Found {len(input_files)} input files")
    
    if args.max_files > 0:
        input_files = input_files[:args.max_files]
        print(f"Processing first {len(input_files)} files")
    
    # Set X509 proxy for xrootd access if needed
    x509_proxy = os.environ.get('X509_USER_PROXY')
    if not x509_proxy and args.xrootd_base:
        # Try default location
        default_proxy = f"/afs/cern.ch/user/x/xcheng/x509up_u180107"
        if os.path.exists(default_proxy):
            os.environ['X509_USER_PROXY'] = default_proxy
            print(f"Set X509_USER_PROXY to {default_proxy}")
    
    analyze_miniaod_files(input_files, args.output, args.mode, args.max_events, args.jobs)


if __name__ == '__main__':
    main()
