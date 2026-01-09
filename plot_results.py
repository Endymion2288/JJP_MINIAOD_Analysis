#!/usr/bin/env python3
"""
Plot the results from JJP (J/psi + J/psi + Phi) gen correlation analysis.
Creates publication-quality plots for delta eta and delta phi correlations.
"""

import ROOT
import os
import sys
import argparse

# Set ROOT style
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetPadLeftMargin(0.12)
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.12)
ROOT.gStyle.SetPadTopMargin(0.08)
ROOT.gStyle.SetTitleOffset(1.2, "X")
ROOT.gStyle.SetTitleOffset(1.4, "Y")
ROOT.gStyle.SetPalette(ROOT.kViridis)


def plot_1d_comparison(fin, hist_names, labels, colors, output_name, title, xlabel, ylabel="Events", logy=False):
    """Plot multiple 1D histograms on the same canvas"""
    c = ROOT.TCanvas("c", "", 800, 600)
    c.SetLogy(logy)
    
    legend = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    
    max_val = 0
    hists = []
    
    for i, (hname, label, color) in enumerate(zip(hist_names, labels, colors)):
        h = fin.Get(hname)
        if not h:
            print(f"Warning: histogram {hname} not found")
            continue
        
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20 + i)
        
        max_val = max(max_val, h.GetMaximum())
        hists.append((h, label))
    
    if not hists:
        return
    
    # Draw
    first = True
    for h, label in hists:
        if logy:
            h.SetMaximum(max_val * 5)
            h.SetMinimum(0.5)
        else:
            h.SetMaximum(max_val * 1.3)
            h.SetMinimum(0)
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle(ylabel)
        h.SetTitle(title)
        
        if first:
            h.Draw("HIST E")
            first = False
        else:
            h.Draw("HIST E SAME")
        
        legend.AddEntry(h, label, "l")
    
    legend.Draw()
    
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def plot_2d(fin, hist_name, output_name, title):
    """Plot a 2D histogram"""
    c = ROOT.TCanvas("c2d", "", 800, 700)
    c.SetLogz(0)
    
    h = fin.Get(hist_name)
    if not h:
        print(f"Warning: histogram {hist_name} not found")
        return
    
    h.SetTitle(title)
    h.GetXaxis().SetTitle("|#Delta y|")
    h.GetYaxis().SetTitle("|#Delta#phi|")
    h.GetZaxis().SetTitle("Events")
    
    h.Draw("COLZ")
    
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def plot_2d_all(fin, output_dir, mode):
    """Plot all three 2D correlations side by side"""
    c = ROOT.TCanvas("c2d_all", "", 1800, 500)
    c.Divide(3, 1)
    
    hist_names = ["h2_dy_dphi_jpsi1_jpsi2", "h2_dy_dphi_jpsi1_phi", "h2_dy_dphi_jpsi2_phi"]
    titles = ["J/#psi_{1} - J/#psi_{2}", "J/#psi_{1} - #phi", "J/#psi_{2} - #phi"]
    
    for i, (hname, title) in enumerate(zip(hist_names, titles)):
        c.cd(i + 1)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetLeftMargin(0.12)
        
        h = fin.Get(hname)
        if not h:
            print(f"Warning: histogram {hname} not found")
            continue
        
        h.SetTitle(f"{title} ({mode})")
        h.GetXaxis().SetTitle("|#Delta y|")
        h.GetYaxis().SetTitle("|#Delta#phi|")
        h.Draw("COLZ")
    
    output_name = os.path.join(output_dir, f"correlation_2d_all_{mode}")
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def plot_kinematics(fin, output_dir, mode):
    """Plot kinematic distributions"""
    # pT distributions
    c = ROOT.TCanvas("c_pt", "", 800, 600)
    
    legend = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    
    h_jpsi1_pt = fin.Get("h_jpsi1_pt")
    h_jpsi2_pt = fin.Get("h_jpsi2_pt")
    h_phi_pt = fin.Get("h_phi_pt")
    
    if h_jpsi1_pt and h_jpsi2_pt and h_phi_pt:
        h_jpsi1_pt.SetLineColor(ROOT.kRed)
        h_jpsi2_pt.SetLineColor(ROOT.kBlue)
        h_phi_pt.SetLineColor(ROOT.kGreen + 2)
        
        for h in [h_jpsi1_pt, h_jpsi2_pt, h_phi_pt]:
            h.SetLineWidth(2)
        
        max_val = max(h_jpsi1_pt.GetMaximum(), h_jpsi2_pt.GetMaximum(), h_phi_pt.GetMaximum())
        h_jpsi1_pt.SetMaximum(max_val * 1.3)
        h_jpsi1_pt.SetMinimum(0)
        h_jpsi1_pt.SetTitle(f"Transverse Momentum Distributions ({mode})")
        h_jpsi1_pt.GetXaxis().SetTitle("p_{T} [GeV]")
        h_jpsi1_pt.GetYaxis().SetTitle("Events")
        
        h_jpsi1_pt.Draw("HIST")
        h_jpsi2_pt.Draw("HIST SAME")
        h_phi_pt.Draw("HIST SAME")
        
        legend.AddEntry(h_jpsi1_pt, "J/#psi_{1}", "l")
        legend.AddEntry(h_jpsi2_pt, "J/#psi_{2}", "l")
        legend.AddEntry(h_phi_pt, "#phi", "l")
        legend.Draw()
        
        output_name = os.path.join(output_dir, f"pt_distributions_{mode}")
        c.SaveAs(output_name + ".pdf")
        c.SaveAs(output_name + ".png")
        print(f"Saved: {output_name}.pdf/.png")
    
    # Eta distributions
    c2 = ROOT.TCanvas("c_eta", "", 800, 600)
    legend2 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    legend2.SetTextSize(0.035)
    
    h_jpsi1_eta = fin.Get("h_jpsi1_eta")
    h_jpsi2_eta = fin.Get("h_jpsi2_eta")
    h_phi_eta = fin.Get("h_phi_eta")
    
    if h_jpsi1_eta and h_jpsi2_eta and h_phi_eta:
        h_jpsi1_eta.SetLineColor(ROOT.kRed)
        h_jpsi2_eta.SetLineColor(ROOT.kBlue)
        h_phi_eta.SetLineColor(ROOT.kGreen + 2)
        
        for h in [h_jpsi1_eta, h_jpsi2_eta, h_phi_eta]:
            h.SetLineWidth(2)
        
        max_val = max(h_jpsi1_eta.GetMaximum(), h_jpsi2_eta.GetMaximum(), h_phi_eta.GetMaximum())
        h_jpsi1_eta.SetMaximum(max_val * 1.3)
        h_jpsi1_eta.SetMinimum(0)
        h_jpsi1_eta.SetTitle(f"Pseudorapidity Distributions ({mode})")
        h_jpsi1_eta.GetXaxis().SetTitle("#eta")
        h_jpsi1_eta.GetYaxis().SetTitle("Events")
        
        h_jpsi1_eta.Draw("HIST")
        h_jpsi2_eta.Draw("HIST SAME")
        h_phi_eta.Draw("HIST SAME")
        
        legend2.AddEntry(h_jpsi1_eta, "J/#psi_{1}", "l")
        legend2.AddEntry(h_jpsi2_eta, "J/#psi_{2}", "l")
        legend2.AddEntry(h_phi_eta, "#phi", "l")
        legend2.Draw()
        
        output_name = os.path.join(output_dir, f"eta_distributions_{mode}")
        c2.SaveAs(output_name + ".pdf")
        c2.SaveAs(output_name + ".png")
        print(f"Saved: {output_name}.pdf/.png")
    
    # Rapidity distributions
    c3 = ROOT.TCanvas("c_y", "", 800, 600)
    legend3 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0)
    legend3.SetTextSize(0.035)
    
    h_jpsi1_y = fin.Get("h_jpsi1_y")
    h_jpsi2_y = fin.Get("h_jpsi2_y")
    h_phi_y = fin.Get("h_phi_y")
    
    if h_jpsi1_y and h_jpsi2_y and h_phi_y:
        h_jpsi1_y.SetLineColor(ROOT.kRed)
        h_jpsi2_y.SetLineColor(ROOT.kBlue)
        h_phi_y.SetLineColor(ROOT.kGreen + 2)
        
        for h in [h_jpsi1_y, h_jpsi2_y, h_phi_y]:
            h.SetLineWidth(2)
        
        max_val = max(h_jpsi1_y.GetMaximum(), h_jpsi2_y.GetMaximum(), h_phi_y.GetMaximum())
        h_jpsi1_y.SetMaximum(max_val * 1.3)
        h_jpsi1_y.SetMinimum(0)
        h_jpsi1_y.SetTitle(f"Rapidity Distributions ({mode})")
        h_jpsi1_y.GetXaxis().SetTitle("y")
        h_jpsi1_y.GetYaxis().SetTitle("Events")
        
        h_jpsi1_y.Draw("HIST")
        h_jpsi2_y.Draw("HIST SAME")
        h_phi_y.Draw("HIST SAME")
        
        legend3.AddEntry(h_jpsi1_y, "J/#psi_{1}", "l")
        legend3.AddEntry(h_jpsi2_y, "J/#psi_{2}", "l")
        legend3.AddEntry(h_phi_y, "#phi", "l")
        legend3.Draw()
        
        output_name = os.path.join(output_dir, f"rapidity_distributions_{mode}")
        c3.SaveAs(output_name + ".pdf")
        c3.SaveAs(output_name + ".png")
        print(f"Saved: {output_name}.pdf/.png")


def plot_invariant_mass(fin, output_dir, mode):
    """Plot invariant mass distributions"""
    c = ROOT.TCanvas("c_mass", "", 1600, 500)
    c.Divide(4, 1)
    
    mass_hists = [
        ("h_mass_jpsi1_jpsi2", "M(J/#psi_{1} + J/#psi_{2})"),
        ("h_mass_jpsi1_phi", "M(J/#psi_{1} + #phi)"),
        ("h_mass_jpsi2_phi", "M(J/#psi_{2} + #phi)"),
        ("h_mass_all", "M(J/#psi_{1} + J/#psi_{2} + #phi)")
    ]
    
    for i, (hname, title) in enumerate(mass_hists):
        c.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.05)
        
        h = fin.Get(hname)
        if not h:
            print(f"Warning: histogram {hname} not found")
            continue
        
        h.SetLineColor(ROOT.kBlue)
        h.SetLineWidth(2)
        h.SetTitle(f"{title} ({mode})")
        h.GetXaxis().SetTitle("M [GeV]")
        h.GetYaxis().SetTitle("Events")
        h.Draw("HIST")
    
    output_name = os.path.join(output_dir, f"invariant_mass_{mode}")
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def main():
    parser = argparse.ArgumentParser(description='Plot JJP gen correlation results')
    parser.add_argument('-i', '--input', 
                        default='gen_correlation_histograms.root',
                        help='Input ROOT file with histograms')
    parser.add_argument('-o', '--output-dir', 
                        default='plots',
                        help='Output directory for plots')
    parser.add_argument('-m', '--mode',
                        default='DPS_1',
                        help='Selection mode label for plot titles')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Open input file
    fin = ROOT.TFile.Open(args.input, "READ")
    if not fin or fin.IsZombie():
        print(f"Error: Cannot open {args.input}")
        sys.exit(1)
    
    print(f"Processing: {args.input}")
    print(f"Output directory: {args.output_dir}")
    print(f"Mode: {args.mode}")
    
    mode = args.mode
    
    # Plot delta y (rapidity) comparisons
    plot_1d_comparison(
        fin,
        ["h_dy_jpsi1_jpsi2", "h_dy_jpsi1_phi", "h_dy_jpsi2_phi"],
        ["J/#psi_{1} - J/#psi_{2}", "J/#psi_{1} - #phi", "J/#psi_{2} - #phi"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        os.path.join(args.output_dir, f"delta_y_comparison_{mode}"),
        f"#Delta y Distributions ({mode})",
        "|#Delta y|"
    )
    
    # Plot delta phi comparisons
    plot_1d_comparison(
        fin,
        ["h_dphi_jpsi1_jpsi2", "h_dphi_jpsi1_phi", "h_dphi_jpsi2_phi"],
        ["J/#psi_{1} - J/#psi_{2}", "J/#psi_{1} - #phi", "J/#psi_{2} - #phi"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        os.path.join(args.output_dir, f"delta_phi_comparison_{mode}"),
        f"#Delta#phi Distributions ({mode})",
        "|#Delta#phi|"
    )
    
    # Plot individual 2D correlations
    plot_2d(fin, "h2_dy_dphi_jpsi1_jpsi2", 
            os.path.join(args.output_dir, f"correlation_2d_jpsi1_jpsi2_{mode}"),
            f"J/#psi_{{1}} - J/#psi_{{2}}: #Delta y vs #Delta#phi ({mode})")
    
    plot_2d(fin, "h2_dy_dphi_jpsi1_phi",
            os.path.join(args.output_dir, f"correlation_2d_jpsi1_phi_{mode}"),
            f"J/#psi_{{1}} - #phi: #Delta y vs #Delta#phi ({mode})")
    
    plot_2d(fin, "h2_dy_dphi_jpsi2_phi",
            os.path.join(args.output_dir, f"correlation_2d_jpsi2_phi_{mode}"),
            f"J/#psi_{{2}} - #phi: #Delta y vs #Delta#phi ({mode})")
    
    # Plot all 2D correlations together
    plot_2d_all(fin, args.output_dir, mode)
    
    # Plot kinematic distributions
    plot_kinematics(fin, args.output_dir, mode)
    
    # Plot invariant mass distributions
    plot_invariant_mass(fin, args.output_dir, mode)
    
    fin.Close()
    print("\nDone!")


if __name__ == '__main__':
    main()
