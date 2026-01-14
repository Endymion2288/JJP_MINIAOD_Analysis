#!/bin/bash
# ==============================================================================
# submit.sh - HTCondor job submission manager for JJPMCAnalyzer
# ==============================================================================
# This script manages HTCondor job submission for JJP Gen-Level analysis.
#
# Usage:
#   ./submit.sh --help                       # Show help
#   ./submit.sh -m DPS_1                     # Submit DPS_1 mode
#   ./submit.sh -m all                       # Submit all modes
#   ./submit.sh --status                     # Check job status
#   ./submit.sh --clean                      # Clean log files
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

msg_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
msg_ok() { echo -e "${GREEN}[OK]${NC} $1"; }
msg_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
msg_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# ==============================================================================
# Help
# ==============================================================================
print_help() {
    cat << EOF
${CYAN}JJPMCAnalyzer HTCondor Submission Manager${NC}

${YELLOW}Usage:${NC}
    $0 [options]
    $0 --status | --clean | --help

${YELLOW}Options:${NC}
    -m, --mode MODE       Selection mode: SPS, DPS_1, DPS_2, TPS
                          Use 'all' to submit all modes (default: DPS_1)
    -j, --jobs N          Number of parallel workers (default: 8)
    -n, --max-events N    Maximum events to process (-1=all)
    --max-files N         Maximum files to process (-1=all)
    -o, --output-dir DIR  Output directory (default: output)
    --dry-run             Show command without submitting
    --flavor FLAVOR       Job flavor (espresso/microcentury/longlunch/workday/tomorrow)

${YELLOW}Management Commands:${NC}
    --status              Show job status (condor_q)
    --history             Show job history
    --clean               Clean log files
    --check-proxy         Check VOMS proxy status

${YELLOW}Examples:${NC}
    $0 -m DPS_1                         # Submit DPS_1 mode
    $0 -m DPS_2 -j 16                   # Submit DPS_2 with 16 workers
    $0 -m all                           # Submit all modes (SPS, DPS_1, DPS_2, TPS)
    $0 -m SPS --max-files 100           # Submit SPS with limited files
    $0 --status                         # Check job status

EOF
}

# ==============================================================================
# Check proxy
# ==============================================================================
check_proxy() {
    msg_info "Checking VOMS proxy..."
    
    if ! command -v voms-proxy-info &>/dev/null; then
        msg_warn "voms-proxy-info not available. Skipping proxy check."
        return 0
    fi
    
    if ! voms-proxy-info --exists &>/dev/null; then
        msg_error "No valid VOMS proxy found!"
        echo ""
        echo "Please create a proxy with:"
        echo "  voms-proxy-init --voms cms --valid 168:00"
        echo ""
        echo "Then copy it for HTCondor access:"
        echo "  cp /tmp/x509up_u\$(id -u) /afs/cern.ch/user/x/xcheng/"
        return 1
    fi
    
    TIMELEFT=$(voms-proxy-info --timeleft 2>/dev/null || echo "0")
    HOURS_LEFT=$((TIMELEFT / 3600))
    
    if [ $HOURS_LEFT -lt 12 ]; then
        msg_warn "Proxy expires in ${HOURS_LEFT} hours. Consider renewing."
    else
        msg_ok "Proxy valid for ${HOURS_LEFT} hours"
    fi
    
    # Check AFS copy
    AFS_PROXY="/afs/cern.ch/user/x/xcheng/x509up_u180107"
    if [ ! -f "$AFS_PROXY" ]; then
        msg_warn "Proxy not found in AFS. Copying..."
        cp "/tmp/x509up_u$(id -u)" "$AFS_PROXY"
        msg_ok "Proxy copied to AFS"
    fi
    
    return 0
}

# ==============================================================================
# Job status
# ==============================================================================
show_status() {
    msg_info "Job status for user: $(whoami)"
    echo ""
    condor_q
}

show_history() {
    msg_info "Recent job history:"
    echo ""
    condor_history -limit 20
}

# ==============================================================================
# Clean logs
# ==============================================================================
clean_logs() {
    msg_info "Cleaning log files..."
    
    if [ -d "logs" ]; then
        rm -rf logs/*.out logs/*.err logs/*.log 2>/dev/null || true
        msg_ok "Cleaned logs directory"
    fi
    
    msg_ok "Log cleanup complete"
}

# ==============================================================================
# Submit job
# ==============================================================================
submit_job() {
    # Parse options
    local MODE="DPS_1"
    local JOBS=""
    local MAX_EVENTS=""
    local MAX_FILES=""
    local OUTPUT_DIR=""
    local DRY_RUN=false
    local FLAVOR=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -m|--mode) MODE="$2"; shift 2;;
            -j|--jobs) JOBS="$2"; shift 2;;
            -n|--max-events) MAX_EVENTS="$2"; shift 2;;
            --max-files) MAX_FILES="$2"; shift 2;;
            -o|--output-dir) OUTPUT_DIR="$2"; shift 2;;
            --dry-run) DRY_RUN=true; shift;;
            --flavor) FLAVOR="$2"; shift 2;;
            *) msg_error "Unknown option: $1"; exit 1;;
        esac
    done
    
    SUB_FILE="jjp_gen.sub"
    
    if [ ! -f "$SUB_FILE" ]; then
        msg_error "Submit file not found: $SUB_FILE"
        exit 1
    fi
    
    # Create log directory
    mkdir -p logs
    
    # Build condor_submit arguments
    local SUBMIT_ARGS=""
    
    [ -n "$JOBS" ] && SUBMIT_ARGS="$SUBMIT_ARGS JOBS=$JOBS"
    [ -n "$MAX_EVENTS" ] && SUBMIT_ARGS="$SUBMIT_ARGS MAX_EVENTS=$MAX_EVENTS"
    [ -n "$MAX_FILES" ] && SUBMIT_ARGS="$SUBMIT_ARGS MAX_FILES=$MAX_FILES"
    [ -n "$OUTPUT_DIR" ] && SUBMIT_ARGS="$SUBMIT_ARGS OUTPUT_DIR=$OUTPUT_DIR"
    [ -n "$FLAVOR" ] && SUBMIT_ARGS="$SUBMIT_ARGS '+JobFlavour=\"$FLAVOR\"'"
    
    # Handle 'all' mode (JJP has SPS, DPS_1, DPS_2, TPS - no DPS_3)
    local MODES_TO_SUBMIT=()
    if [ "${MODE,,}" = "all" ]; then
        MODES_TO_SUBMIT=(SPS DPS_1 DPS_2 TPS)
    else
        MODES_TO_SUBMIT=("$MODE")
    fi
    
    # Submit jobs
    for m in "${MODES_TO_SUBMIT[@]}"; do
        local CMD="condor_submit $SUB_FILE MODE=$m $SUBMIT_ARGS"
        
        echo ""
        msg_info "Submitting: JJP Gen-Level $m"
        echo "Command: $CMD"
        
        if [ "$DRY_RUN" = true ]; then
            msg_warn "Dry run - not submitting"
        else
            eval $CMD
            msg_ok "Job submitted"
        fi
    done
}

# ==============================================================================
# Main
# ==============================================================================
main() {
    if [ $# -eq 0 ]; then
        print_help
        exit 0
    fi
    
    case "$1" in
        -h|--help)
            print_help
            exit 0
            ;;
        --status)
            show_status
            exit 0
            ;;
        --history)
            show_history
            exit 0
            ;;
        --clean)
            clean_logs
            exit 0
            ;;
        --check-proxy)
            check_proxy
            exit $?
            ;;
        -m|--mode|-j|--jobs|-n|--max-events|--max-files|-o|--output-dir|--dry-run|--flavor)
            check_proxy || exit 1
            submit_job "$@"
            ;;
        *)
            msg_error "Unknown command: $1"
            print_help
            exit 1
            ;;
    esac
}

main "$@"
