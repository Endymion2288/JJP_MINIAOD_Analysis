#!/usr/bin/env python3
"""
Analyze MiniAOD prunedGenParticles to study correlations between
Jpsi1, Jpsi2, and phi particles.

Selection criteria:
- Jpsi2 and phi must come from the same parton-level interaction
- phi comes from a gluon decay
- The gluon and Jpsi2 are produced in the same parton-parton interaction

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

def get_hard_process_ancestors(particle):
    """Get ancestors from hard process (status 21-29 or initial status)"""
    ancestors = []
    mothers = get_mother_chain(particle)
    for m in mothers:
        # Hard process particles typically have status 21-29 in Pythia8
        # or status 3 in older generators
        if m.status() == 21 or m.status() == 22 or m.status() == 23 or m.status() == 3:
            ancestors.append(m)
    return ancestors

def find_common_ancestor(particle1, particle2, chain_cache, ancestor_cache):
    """Find if two particles share a common ancestor using cached chains."""

    def ancestor_set(p):
        key = id(p)
        if key in ancestor_cache:
            return ancestor_cache[key]
        moms = get_mother_chain(p, chain_cache)
        aset = set((m.pdgId(), round(m.pt(), 2), round(m.eta(), 3), round(m.phi(), 3)) for m in moms)
        ancestor_cache[key] = aset
        return aset

    a1 = ancestor_set(particle1)
    a2 = ancestor_set(particle2)
    common = a1.intersection(a2)
    return len(common) > 0, common

def phi_from_gluon(phi_particle, chain_cache):
    """Check if phi comes from gluon decay and return the gluon (cached)."""
    mothers = get_mother_chain(phi_particle, chain_cache)
    for m in mothers:
        if abs(m.pdgId()) == PDG_GLUON:
            return True, m
    return False, None


def process_file_batch(file_batch, batch_id, output_dir, max_events_per_batch=-1):
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
    h_dy_12 = ROOT.TH1F("h_dy_12", "Delta y (Jpsi1-Jpsi2);|#Delta y|;Events", 50, 0, 5)
    h_dy_1phi = ROOT.TH1F("h_dy_1phi", "Delta y (Jpsi1-Phi);|#Delta y|;Events", 50, 0, 5)
    h_dy_2phi = ROOT.TH1F("h_dy_2phi", "Delta y (Jpsi2-Phi);|#Delta y|;Events", 50, 0, 5)
    
    h_dphi_12 = ROOT.TH1F("h_dphi_12", "Delta #phi (Jpsi1-Jpsi2);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    h_dphi_1phi = ROOT.TH1F("h_dphi_1phi", "Delta #phi (Jpsi1-Phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    h_dphi_2phi = ROOT.TH1F("h_dphi_2phi", "Delta #phi (Jpsi2-Phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    
    h2_12 = ROOT.TH2F("h2_dy_dphi_12", "Jpsi1-Jpsi2: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                       50, 0, 5, 50, 0, ROOT.TMath.Pi())
    h2_1phi = ROOT.TH2F("h2_dy_dphi_1phi", "Jpsi1-Phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                         50, 0, 5, 50, 0, ROOT.TMath.Pi())
    h2_2phi = ROOT.TH2F("h2_dy_dphi_2phi", "Jpsi2-Phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                         50, 0, 5, 50, 0, ROOT.TMath.Pi())
    
    h_jpsi1_pt = ROOT.TH1F("h_jpsi1_pt", "Jpsi1 p_{T};p_{T} [GeV];Events", 100, 0, 50)
    h_jpsi2_pt = ROOT.TH1F("h_jpsi2_pt", "Jpsi2 p_{T};p_{T} [GeV];Events", 100, 0, 50)
    h_phi_pt = ROOT.TH1F("h_phi_pt", "Phi p_{T};p_{T} [GeV];Events", 100, 0, 50)
    
    h_jpsi1_eta = ROOT.TH1F("h_jpsi1_eta", "Jpsi1 #eta;#eta;Events", 60, -3, 3)
    h_jpsi2_eta = ROOT.TH1F("h_jpsi2_eta", "Jpsi2 #eta;#eta;Events", 60, -3, 3)
    h_phi_eta = ROOT.TH1F("h_phi_eta", "Phi #eta;#eta;Events", 60, -3, 3)
    
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
            
            if len(jpsis) < 2 or len(phis) < 1:
                continue

            jpsis.sort(key=lambda x: x.pt(), reverse=True)
            phis.sort(key=lambda x: x.pt(), reverse=True)

            selected_jpsi1 = None
            selected_jpsi2 = None
            selected_phi = None

            mother_chain_cache = {}
            ancestor_cache = {}

            for phi_cand in phis:
                if phi_cand.pt() <= 4.0:
                    continue
                is_from_gluon, gluon = phi_from_gluon(phi_cand, mother_chain_cache)
                if not is_from_gluon:
                    continue

                for jpsi_cand in jpsis:
                    has_common, _ = find_common_ancestor(jpsi_cand, gluon, mother_chain_cache, ancestor_cache)
                    if has_common:
                        selected_jpsi1 = jpsi_cand
                        selected_phi = phi_cand
                        break

                if selected_jpsi1 is not None:
                    break

            if selected_jpsi1 is not None:
                for jpsi_cand in jpsis:
                    if jpsi_cand == selected_jpsi1:
                        continue
                    selected_jpsi2 = jpsi_cand
                    break

            if selected_jpsi1 is None or selected_jpsi2 is None or selected_phi is None:
                continue
            
            n_selected += 1
            
            jpsi1 = selected_jpsi1
            jpsi2 = selected_jpsi2
            phi = selected_phi
            
            y_jpsi1 = rapidity(jpsi1)
            y_jpsi2 = rapidity(jpsi2)
            y_phi = rapidity(phi)
            
            dy_12 = abs(y_jpsi1 - y_jpsi2)
            dy_1phi = abs(y_jpsi1 - y_phi)
            dy_2phi = abs(y_jpsi2 - y_phi)
            
            dphi_12 = abs(delta_phi(jpsi1.phi(), jpsi2.phi()))
            dphi_1phi = abs(delta_phi(jpsi1.phi(), phi.phi()))
            dphi_2phi = abs(delta_phi(jpsi2.phi(), phi.phi()))
            
            h_dy_12.Fill(dy_12)
            h_dy_1phi.Fill(dy_1phi)
            h_dy_2phi.Fill(dy_2phi)
            
            h_dphi_12.Fill(dphi_12)
            h_dphi_1phi.Fill(dphi_1phi)
            h_dphi_2phi.Fill(dphi_2phi)
            
            h2_12.Fill(dy_12, dphi_12)
            h2_1phi.Fill(dy_1phi, dphi_1phi)
            h2_2phi.Fill(dy_2phi, dphi_2phi)
            
            h_jpsi1_pt.Fill(jpsi1.pt())
            h_jpsi2_pt.Fill(jpsi2.pt())
            h_phi_pt.Fill(phi.pt())
            
            h_jpsi1_eta.Fill(jpsi1.eta())
            h_jpsi2_eta.Fill(jpsi2.eta())
            h_phi_eta.Fill(phi.eta())
        
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

def analyze_miniaod_files(input_files, output_file, max_events=-1, n_workers=1):
    """
    Main analysis function to process MiniAOD files.
    Supports parallel processing with n_workers > 1.
    """
    if n_workers > 1:
        # Parallel processing mode
        print(f"Running in parallel mode with {n_workers} workers")
        
        # Create temporary directory for intermediate files
        temp_dir = tempfile.mkdtemp(prefix="gen_corr_")
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
        # Note: we use spawn to avoid issues with ROOT in forked processes
        ctx = mp.get_context('spawn')
        
        results = []
        with ctx.Pool(processes=n_workers) as pool:
            tasks = []
            for batch_id, batch in enumerate(file_batches):
                task = pool.apply_async(
                    process_file_batch,
                    args=(batch, batch_id, temp_dir, max_events_per_batch)
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
        
        print(f"\n=== Analysis Summary ===")
        print(f"Total events processed: {total_events}")
        print(f"Events with valid Jpsi1, Jpsi2, Phi selection: {selected_events}")
        print(f"Selection efficiency: {100.0*selected_events/total_events:.2f}%" if total_events > 0 else "N/A")
        print(f"\nOutput saved to: {output_file}")
        
    else:
        # Sequential processing mode (original behavior)
        analyze_miniaod_files_sequential(input_files, output_file, max_events)


def analyze_miniaod_files_sequential(input_files, output_file, max_events=-1):
    """
    Main analysis function to process MiniAOD files
    """
    # Load FWLite
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    
    from DataFormats.FWLite import Events, Handle
    
    # Create output file and histograms
    fout = ROOT.TFile(output_file, "RECREATE")
    
    # 1D histograms for delta y (rapidity difference)
    h_dy_12 = ROOT.TH1F("h_dy_12", "Delta y (Jpsi1-Jpsi2);|#Delta y|;Events", 50, 0, 5)
    h_dy_1phi = ROOT.TH1F("h_dy_1phi", "Delta y (Jpsi1-Phi);|#Delta y|;Events", 50, 0, 5)
    h_dy_2phi = ROOT.TH1F("h_dy_2phi", "Delta y (Jpsi2-Phi);|#Delta y|;Events", 50, 0, 5)
    
    # 1D histograms for delta phi
    h_dphi_12 = ROOT.TH1F("h_dphi_12", "Delta #phi (Jpsi1-Jpsi2);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    h_dphi_1phi = ROOT.TH1F("h_dphi_1phi", "Delta #phi (Jpsi1-Phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    h_dphi_2phi = ROOT.TH1F("h_dphi_2phi", "Delta #phi (Jpsi2-Phi);|#Delta#phi|;Events", 50, 0, ROOT.TMath.Pi())
    
    # 2D histograms: delta y vs delta phi
    h2_12 = ROOT.TH2F("h2_dy_dphi_12", "Jpsi1-Jpsi2: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                       50, 0, 5, 50, 0, ROOT.TMath.Pi())
    h2_1phi = ROOT.TH2F("h2_dy_dphi_1phi", "Jpsi1-Phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                         50, 0, 5, 50, 0, ROOT.TMath.Pi())
    h2_2phi = ROOT.TH2F("h2_dy_dphi_2phi", "Jpsi2-Phi: #Delta y vs #Delta#phi;|#Delta y|;|#Delta#phi|", 
                         50, 0, 5, 50, 0, ROOT.TMath.Pi())
    
    # Kinematic distributions
    h_jpsi1_pt = ROOT.TH1F("h_jpsi1_pt", "Jpsi1 p_{T};p_{T} [GeV];Events", 100, 0, 50)
    h_jpsi2_pt = ROOT.TH1F("h_jpsi2_pt", "Jpsi2 p_{T};p_{T} [GeV];Events", 100, 0, 50)
    h_phi_pt = ROOT.TH1F("h_phi_pt", "Phi p_{T};p_{T} [GeV];Events", 100, 0, 50)
    
    h_jpsi1_eta = ROOT.TH1F("h_jpsi1_eta", "Jpsi1 #eta;#eta;Events", 60, -3, 3)
    h_jpsi2_eta = ROOT.TH1F("h_jpsi2_eta", "Jpsi2 #eta;#eta;Events", 60, -3, 3)
    h_phi_eta = ROOT.TH1F("h_phi_eta", "Phi #eta;#eta;Events", 60, -3, 3)
    
    # Handle for prunedGenParticles
    handle_genParticles = Handle("std::vector<reco::GenParticle>")
    label_genParticles = ("prunedGenParticles", "", "PAT")
    
    n_total = 0
    n_selected = 0
    
    for fname in input_files:
        print(f"Processing: {fname}")
        
        events = Events(fname)
        
        for i, event in enumerate(events):
            if max_events > 0 and n_total >= max_events:
                break
                
            n_total += 1
            if n_total % 1000 == 0:
                print(f"  Processed {n_total} events...")
            
            # Get gen particles
            event.getByLabel(label_genParticles, handle_genParticles)
            if not handle_genParticles.isValid():
                continue
            
            genParticles = handle_genParticles.product()
            
            # Find J/psi and phi particles (status 2 = final state before decay)
            jpsis = []
            phis = []
            
            for gp in genParticles:
                pdg = abs(gp.pdgId())
                
                # Look for J/psi: status can be 2 (decayed) or we need to check isLastCopy()
                if pdg == PDG_JPSI:
                    # In MiniAOD, use isLastCopy to get the final J/psi
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        jpsis.append(gp)
                    elif gp.status() == 2:
                        jpsis.append(gp)
                
                # Look for phi meson
                if pdg == PDG_PHI:
                    if hasattr(gp, 'isLastCopy') and gp.isLastCopy():
                        phis.append(gp)
                    elif gp.status() == 2:
                        phis.append(gp)
            
            # Need at least 2 J/psi and 1 phi
            if len(jpsis) < 2 or len(phis) < 1:
                continue

            # Sort by pT (highest first) so we pick the hardest objects when multiple match
            jpsis.sort(key=lambda x: x.pt(), reverse=True)
            phis.sort(key=lambda x: x.pt(), reverse=True)

            # New definition:
            #   Jpsi1 = a J/psi that shares a hard-process ancestor with a phi from gluon decay (phi pT > 4 GeV)
            #   Jpsi2 = any other J/psi (different from Jpsi1), prefer highest pT remaining
            selected_jpsi1 = None
            selected_jpsi2 = None
            selected_phi = None

            # Per-event caches to speed up mother/ancestor traversals
            mother_chain_cache = {}
            ancestor_cache = {}

            # First, find a phi from gluon with pT > 4 and a J/psi sharing the gluon ancestry
            for phi_cand in phis:
                if phi_cand.pt() <= 4.0:
                    continue
                is_from_gluon, gluon = phi_from_gluon(phi_cand, mother_chain_cache)
                if not is_from_gluon:
                    continue

                for jpsi_cand in jpsis:
                    has_common, _ = find_common_ancestor(jpsi_cand, gluon, mother_chain_cache, ancestor_cache)
                    if has_common:
                        selected_jpsi1 = jpsi_cand
                        selected_phi = phi_cand
                        break

                if selected_jpsi1 is not None:
                    break

            # Pick Jpsi2 as a different J/psi (highest pT among remaining)
            if selected_jpsi1 is not None:
                for jpsi_cand in jpsis:
                    if jpsi_cand == selected_jpsi1:
                        continue
                    selected_jpsi2 = jpsi_cand
                    break

            if selected_jpsi1 is None or selected_jpsi2 is None or selected_phi is None:
                continue
            
            n_selected += 1
            
            # Calculate kinematic variables
            jpsi1 = selected_jpsi1
            jpsi2 = selected_jpsi2
            phi = selected_phi
            
            # Delta y (rapidity difference, absolute value)
            y_jpsi1 = rapidity(jpsi1)
            y_jpsi2 = rapidity(jpsi2)
            y_phi = rapidity(phi)
            
            dy_12 = abs(y_jpsi1 - y_jpsi2)
            dy_1phi = abs(y_jpsi1 - y_phi)
            dy_2phi = abs(y_jpsi2 - y_phi)
            
            # Delta phi (absolute value, wrapped)
            dphi_12 = abs(delta_phi(jpsi1.phi(), jpsi2.phi()))
            dphi_1phi = abs(delta_phi(jpsi1.phi(), phi.phi()))
            dphi_2phi = abs(delta_phi(jpsi2.phi(), phi.phi()))
            
            # Fill 1D histograms
            h_dy_12.Fill(dy_12)
            h_dy_1phi.Fill(dy_1phi)
            h_dy_2phi.Fill(dy_2phi)
            
            h_dphi_12.Fill(dphi_12)
            h_dphi_1phi.Fill(dphi_1phi)
            h_dphi_2phi.Fill(dphi_2phi)
            
            # Fill 2D histograms
            h2_12.Fill(dy_12, dphi_12)
            h2_1phi.Fill(dy_1phi, dphi_1phi)
            h2_2phi.Fill(dy_2phi, dphi_2phi)
            
            # Fill kinematic distributions
            h_jpsi1_pt.Fill(jpsi1.pt())
            h_jpsi2_pt.Fill(jpsi2.pt())
            h_phi_pt.Fill(phi.pt())
            
            h_jpsi1_eta.Fill(jpsi1.eta())
            h_jpsi2_eta.Fill(jpsi2.eta())
            h_phi_eta.Fill(phi.eta())
        
        if max_events > 0 and n_total >= max_events:
            break
    
    print(f"\n=== Analysis Summary ===")
    print(f"Total events processed: {n_total}")
    print(f"Events with valid Jpsi1, Jpsi2, Phi selection: {n_selected}")
    print(f"Selection efficiency: {100.0*n_selected/n_total:.2f}%" if n_total > 0 else "N/A")
    
    # Write and close
    fout.Write()
    fout.Close()
    print(f"\nOutput saved to: {output_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze MiniAOD gen particles')
    parser.add_argument('-i', '--input-dir', 
                        default='/eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/MINIAOD/',
                        help='Input directory containing MiniAOD files')
    parser.add_argument('-o', '--output', 
                        default='gen_correlation_histograms.root',
                        help='Output ROOT file')
    parser.add_argument('-n', '--max-events', type=int, default=-1,
                        help='Maximum number of events to process (-1 for all)')
    parser.add_argument('--max-files', type=int, default=-1,
                        help='Maximum number of files to process (-1 for all)')
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help='Number of parallel workers (default: 1, sequential)')
    
    args = parser.parse_args()
    
    # Get input files
    input_pattern = os.path.join(args.input_dir, '*.root')
    input_files = sorted(glob.glob(input_pattern))
    
    if not input_files:
        print(f"Error: No ROOT files found in {args.input_dir}")
        sys.exit(1)
    
    print(f"Found {len(input_files)} input files")
    
    if args.max_files > 0:
        input_files = input_files[:args.max_files]
        print(f"Processing first {len(input_files)} files")
    
    analyze_miniaod_files(input_files, args.output, args.max_events, args.jobs)


if __name__ == '__main__':
    main()
