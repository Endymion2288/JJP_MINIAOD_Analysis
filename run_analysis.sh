#!/bin/bash
# Run the gen-level correlation analysis for J/psi + J/psi + Phi
# Supports SPS, DPS_1, DPS_2, TPS selection modes

# Default values
MODE="DPS_1"
MAX_FILES=-1
MAX_EVENTS=-1
N_JOBS=1
OUTPUT_DIR="output"
XROOTD_BASE=""
INPUT_DIR=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        --max-files)
            MAX_FILES="$2"
            shift 2
            ;;
        -n|--max-events)
            MAX_EVENTS="$2"
            shift 2
            ;;
        -j|--jobs)
            N_JOBS="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --xrootd-base)
            XROOTD_BASE="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  -m, --mode MODE       Selection mode: SPS, DPS_1, DPS_2, TPS (default: DPS_1)"
            echo "  --max-files N         Maximum number of files to process (-1 for all)"
            echo "  -n, --max-events N    Maximum number of events to process (-1 for all)"
            echo "  -j, --jobs N          Number of parallel workers (default: 1)"
            echo "  -o, --output-dir DIR  Output directory (default: output)"
            echo "  --xrootd-base PATH    External xrootd input path (overrides default)"
            echo "  -i, --input-dir PATH  Local/EOS input directory with MiniAOD files (overrides xrootd)"
            echo "  -h, --help            Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Setup environment
echo "=========================================="
echo "Setting up CMSSW environment..."
echo "=========================================="

cd /eos/user/x/xcheng/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`

# Set X509 proxy for xrootd access
export X509_USER_PROXY=/afs/cern.ch/user/x/xcheng/x509up_u180107

# Go to analysis directory
cd JJPMCAnalyzer

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Determine input source (local directory has highest priority)
ANALYZE_INPUT_ARGS=""
INPUT_DESC=""

if [ -n "${INPUT_DIR}" ]; then
    ANALYZE_INPUT_ARGS=( -i "${INPUT_DIR}" )
    INPUT_DESC="Local/EOS directory: ${INPUT_DIR}"
else
    # Define input base path on IHEP xrootd (if not externally provided)
    if [ -z "${XROOTD_BASE}" ]; then
        XROOTD_BASE_DEFAULT="/eos/ihep/cms/store/user/xcheng/MC_Production/output"
        
        case ${MODE} in
            SPS)
                INPUT_PATH="${XROOTD_BASE_DEFAULT}/JJP_SPS"
                ;;
            DPS_1)
                INPUT_PATH="${XROOTD_BASE_DEFAULT}/JJP_DPS1"
                ;;
            DPS_2)
                INPUT_PATH="${XROOTD_BASE_DEFAULT}/JJP_DPS2"
                ;;
            TPS)
                INPUT_PATH="${XROOTD_BASE_DEFAULT}/JJP_TPS"
                ;;
            *)
                echo "Unknown mode: ${MODE}"
                exit 1
                ;;
        esac
    else
        INPUT_PATH="${XROOTD_BASE}"
    fi
    ANALYZE_INPUT_ARGS=( --xrootd-base "${INPUT_PATH}" )
    INPUT_DESC="xrootd path: ${INPUT_PATH}"
fi

OUTPUT_ROOT="${OUTPUT_DIR}/gen_correlation_${MODE}.root"
PLOT_DIR="${OUTPUT_DIR}/plots_${MODE}"

echo ""
echo "=========================================="
echo "Running Gen-Level Correlation Analysis"
echo "=========================================="
echo "Mode: ${MODE}"
echo "Input: ${INPUT_DESC}"
echo "Output ROOT: ${OUTPUT_ROOT}"
echo "Plot directory: ${PLOT_DIR}"
echo "Max files: ${MAX_FILES}"
echo "Max events: ${MAX_EVENTS}"
echo "Parallel jobs: ${N_JOBS}"
echo "=========================================="
echo ""

# Run analysis
python3 analyze_gen_correlations.py \
    "${ANALYZE_INPUT_ARGS[@]}" \
    -o "${OUTPUT_ROOT}" \
    -m "${MODE}" \
    -n ${MAX_EVENTS} \
    --max-files ${MAX_FILES} \
    -j ${N_JOBS}

# Check if analysis succeeded
if [ $? -ne 0 ]; then
    echo "Analysis failed!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Creating Plots"
echo "=========================================="

python3 plot_results.py \
    -i "${OUTPUT_ROOT}" \
    -o "${PLOT_DIR}" \
    -m "${MODE}"

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
echo "Output histograms: ${OUTPUT_ROOT}"
echo "Output plots: ${PLOT_DIR}/"
