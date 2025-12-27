#!/usr/bin/env python3
"""
Plot the results from gen_correlation_histograms.root
Creates publication-quality plots for delta eta and delta phi correlations.
"""

import ROOT
import os
import sys

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


def plot_1d_comparison(fin, hist_names, labels, colors, output_name, title, xlabel, ylabel="Events"):
    """Plot multiple 1D histograms on the same canvas"""
    c = ROOT.TCanvas("c", "", 800, 600)
    c.SetLogy(0)
    
    legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
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
        h.SetMaximum(max_val * 1.3)
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
    h.GetXaxis().SetTitle("|#Delta#eta|")
    h.GetYaxis().SetTitle("|#Delta#phi|")
    h.GetZaxis().SetTitle("Events")
    
    h.Draw("COLZ")
    
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def plot_2d_all(fin, output_dir):
    """Plot all three 2D correlations side by side"""
    c = ROOT.TCanvas("c2d_all", "", 1800, 500)
    c.Divide(3, 1)
    
    hist_names = ["h2_deta_dphi_12", "h2_deta_dphi_1phi", "h2_deta_dphi_2phi"]
    titles = ["J/#psi_{1} - J/#psi_{2}", "J/#psi_{1} - #phi", "J/#psi_{2} - #phi"]
    
    for i, (hname, title) in enumerate(zip(hist_names, titles)):
        c.cd(i + 1)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetLeftMargin(0.12)
        
        h = fin.Get(hname)
        if not h:
            print(f"Warning: histogram {hname} not found")
            continue
        
        h.SetTitle(title)
        h.GetXaxis().SetTitle("|#Delta#eta|")
        h.GetYaxis().SetTitle("|#Delta#phi|")
        h.Draw("COLZ")
    
    output_name = os.path.join(output_dir, "correlation_2d_all")
    c.SaveAs(output_name + ".pdf")
    c.SaveAs(output_name + ".png")
    print(f"Saved: {output_name}.pdf/.png")


def plot_kinematics(fin, output_dir):
    """Plot kinematic distributions"""
    # pT distributions
    c = ROOT.TCanvas("c_pt", "", 800, 600)
    
    legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
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
        h_jpsi1_pt.SetTitle("Transverse Momentum Distributions")
        h_jpsi1_pt.GetXaxis().SetTitle("p_{T} [GeV]")
        h_jpsi1_pt.GetYaxis().SetTitle("Events")
        
        h_jpsi1_pt.Draw("HIST")
        h_jpsi2_pt.Draw("HIST SAME")
        h_phi_pt.Draw("HIST SAME")
        
        legend.AddEntry(h_jpsi1_pt, "J/#psi_{1}", "l")
        legend.AddEntry(h_jpsi2_pt, "J/#psi_{2}", "l")
        legend.AddEntry(h_phi_pt, "#phi", "l")
        legend.Draw()
        
        output_name = os.path.join(output_dir, "pt_distributions")
        c.SaveAs(output_name + ".pdf")
        c.SaveAs(output_name + ".png")
        print(f"Saved: {output_name}.pdf/.png")
    
    # Eta distributions
    c2 = ROOT.TCanvas("c_eta", "", 800, 600)
    legend2 = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    
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
        h_jpsi1_eta.SetTitle("Pseudorapidity Distributions")
        h_jpsi1_eta.GetXaxis().SetTitle("#eta")
        h_jpsi1_eta.GetYaxis().SetTitle("Events")
        
        h_jpsi1_eta.Draw("HIST")
        h_jpsi2_eta.Draw("HIST SAME")
        h_phi_eta.Draw("HIST SAME")
        
        legend2.AddEntry(h_jpsi1_eta, "J/#psi_{1}", "l")
        legend2.AddEntry(h_jpsi2_eta, "J/#psi_{2}", "l")
        legend2.AddEntry(h_phi_eta, "#phi", "l")
        legend2.Draw()
        
        output_name = os.path.join(output_dir, "eta_distributions")
        c2.SaveAs(output_name + ".pdf")
        c2.SaveAs(output_name + ".png")
        print(f"Saved: {output_name}.pdf/.png")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot gen correlation results')
    parser.add_argument('-i', '--input', 
                        default='gen_correlation_histograms.root',
                        help='Input ROOT file with histograms')
    parser.add_argument('-o', '--output-dir', 
                        default='plots',
                        help='Output directory for plots')
    
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
    
    # Plot delta eta comparisons
    plot_1d_comparison(
        fin,
        ["h_deta_12", "h_deta_1phi", "h_deta_2phi"],
        ["J/#psi_{1}-J/#psi_{2}", "J/#psi_{1}-#phi", "J/#psi_{2}-#phi"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        os.path.join(args.output_dir, "delta_eta_comparison"),
        "#Delta#eta Distributions",
        "|#Delta#eta|"
    )
    
    # Plot delta phi comparisons
    plot_1d_comparison(
        fin,
        ["h_dphi_12", "h_dphi_1phi", "h_dphi_2phi"],
        ["J/#psi_{1}-J/#psi_{2}", "J/#psi_{1}-#phi", "J/#psi_{2}-#phi"],
        [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2],
        os.path.join(args.output_dir, "delta_phi_comparison"),
        "#Delta#phi Distributions",
        "|#Delta#phi|"
    )
    
    # Plot individual 2D correlations
    plot_2d(fin, "h2_deta_dphi_12", 
            os.path.join(args.output_dir, "correlation_2d_jpsi1_jpsi2"),
            "J/#psi_{1} - J/#psi_{2}: #Delta#eta vs #Delta#phi")
    
    plot_2d(fin, "h2_deta_dphi_1phi",
            os.path.join(args.output_dir, "correlation_2d_jpsi1_phi"),
            "J/#psi_{1} - #phi: #Delta#eta vs #Delta#phi")
    
    plot_2d(fin, "h2_deta_dphi_2phi",
            os.path.join(args.output_dir, "correlation_2d_jpsi2_phi"),
            "J/#psi_{2} - #phi: #Delta#eta vs #Delta#phi")
    
    # Plot all 2D correlations together
    plot_2d_all(fin, args.output_dir)
    
    # Plot kinematic distributions
    plot_kinematics(fin, args.output_dir)
    
    fin.Close()
    print("\nDone!")


if __name__ == '__main__':
    main()
