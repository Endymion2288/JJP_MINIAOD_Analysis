#!/bin/bash
# Run all analysis modes for J/psi + J/psi + Phi

# Setup environment
cd /eos/user/x/xcheng/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/x/xcheng/x509up_u180107

cd JJPMCAnalyzer

# Create main output directory
mkdir -p output

echo "=========================================="
echo "Running ALL JJP analysis modes"
echo "=========================================="

# List of all modes
MODES="SPS DPS_1 DPS_2 TPS"

for MODE in ${MODES}; do
    echo ""
    echo "=========================================="
    echo "Processing mode: ${MODE}"
    echo "=========================================="
    
    ./run_analysis.sh -m ${MODE} --max-files 10 -n 10000 -j 4 -o output
    
    echo ""
done

echo ""
echo "=========================================="
echo "All modes complete!"
echo "=========================================="
echo "Results in: output/"
ls -la output/
