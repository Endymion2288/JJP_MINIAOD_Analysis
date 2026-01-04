#!/bin/bash
# Run the gen-level correlation analysis on MiniAOD files (SPS selection)

# Setup CMSSW environment
cd /eos/home-x/xcheng/learn_MC/loopmix_pythia/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`

# Go to analysis directory
cd JJP_SPS_MINIAOD_Analysis

# Run analysis (adjust parameters as needed)
# -n: max events (-1 for all)
# --max-files: max number of files to process (-1 for all)
# -j: number of parallel workers

echo "=========================================="
echo "Running Gen-Level Correlation Analysis (SPS)"
echo "=========================================="

python3 analyze_gen_correlations.py \
    -i /eos/user/x/xcheng/learn_MC/JJP_SPS_MC_output/MINIAOD/ \
    -o gen_correlation_histograms_sps.root \
    -n -1 \
    --max-files 10

echo ""
echo "=========================================="
echo "Creating Plots"
echo "=========================================="

python3 plot_results.py \
    -i gen_correlation_histograms_sps.root \
    -o plots

echo ""
echo "Analysis complete!"
echo "Output histograms: gen_correlation_histograms_sps.root"
echo "Output plots: plots/"
