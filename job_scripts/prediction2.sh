#!/usr/bin/env bash
# =============================================================================
# prediction2.sh
# =============================================================================
# Runs all Prediction 2 simulations for "Reptating polymers and integer
# partitions: a Bose-Einstein correspondence".
#
# Each run:
#   - pipes stdout+stderr to a uniquely named log file
#   - renames the output CSV to a uniquely named file
#   - prints progress to the terminal
#
# Usage:
#   ./prediction2.sh
#
# To run in background and log everything:
#   nohup ./prediction2.sh > master.log 2>&1 &
#
# Output files produced:
#   ./simple/           simple test run
#   ./primary/          primary scaling run
#   ./density_low/      obstacle density rho=0.02
#   ./density_high/     obstacle density rho=0.10
#   ./seed_check_1/     seed=100 reproducibility check
#   ./seed_check_2/     seed=200 reproducibility check
#   ./extended/         L=1000-2000
# =============================================================================


#SBATCH -o prediction2.%j.out
#SBATCH -e prediction2.%j.err
#SBATCH -J prediction2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=48
#SBATCH --partition=CPU

set -euo pipefail

WORKDIR="$HOME/rp_production/prediction2" 
PROGRAM="$WORKDIR/prediction2.py"
RESULTS_DIR="$WORKDIR"
LOGS_DIR="$WORKDIR"

# --- create output directories if they don't exist ---
mkdir -p "$RESULTS_DIR/prediction2_simple"
mkdir -p "$RESULTS_DIR/prediction2_primary"
mkdir -p "$RESULTS_DIR/prediction2_density_low"
mkdir -p "$RESULTS_DIR/prediction2_density_high"
mkdir -p "$RESULTS_DIR/prediction2_seed_check_1"
mkdir -p "$RESULTS_DIR/prediction2_seed_check_2"
#mkdir -p "$RESULTS_DIR/prediction2_extended"

# --- helper function ---
# run_sim LABEL OUTDIR EXTRA_ARGS...
#   LABEL    short name used in log filename and progress messages
#   OUTDIR   directory to move the CSV into after the run
#   ...      all remaining args passed to prediction2.py
run_sim() {
    local label="$1"
    local outdir="$2"
    shift 2

    local csv_out="$outdir/${label}.csv"
    local log_out="$LOGS_DIR/${label}.log"

    # Extract the --outfile value from the remaining args
    local outfile=""
    local args_copy=("$@")
    for ((i=0; i<${#args_copy[@]}; i++)); do
        if [[ "${args_copy[$i]}" == "--outfile" ]]; then
            outfile="${args_copy[$((i+1))]}"
            break
        fi
    done
    local csv_src="$(dirname "$PROGRAM")/$outfile"

    echo ""
    echo "============================================================"
    echo "Starting run: $label"
    echo "  Args: $*"
    echo "  Log:  $log_out"
    echo "  CSV:  $csv_out"
    echo "  Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"

    # Run simulation; tee mirrors output to terminal and log file
    python3 "$PROGRAM" "$@" 2>&1 | tee "$log_out"

    # Move CSV from script directory to named output location
    if [ -f "$csv_src" ]; then
        mv "$csv_src" "$csv_out"
        echo "  CSV saved to: $csv_out"
    else
        echo "  WARNING: expected CSV $csv_src not found after run."
    fi

    echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
}

# =============================================================================
# SIMPLE TEST RUN
# Core dataset: L = 20-200, rho = 0.05
# =============================================================================
run_sim "simple" "$RESULTS_DIR/prediction2_simple" \
    --Lvals   20,30,50,75,100,150,200 \
    --lattice 200 \
    --obstacle_density 0.05 \
    --n_chains 8 \
    --n_runs   4 \
    --seed     42 \
    --outfile  prediction2_results_simple.csv

# =============================================================================
# PRIMARY SCALING RUN
# Core dataset: L = 100-1000, rho = 0.05, 128 independent chains
# =============================================================================
run_sim "primary" "$RESULTS_DIR/prediction2_primary" \
    --Lvals   100,150,200,300,400,500,750,1000 \
    --lattice 2000 \
    --obstacle_density 0.05 \
    --n_chains 16 \
    --n_runs   8 \
    --seed     42 \
    --outfile  prediction2_results_primary.csv

# =============================================================================
# OBSTACLE DENSITY SENSITIVITY: rho = 0.02 (wider tube)
# Checks universality of slope with respect to tube width
# =============================================================================
run_sim "density_low" "$RESULTS_DIR/prediction2_density_low" \
    --Lvals   100,200,300,500 \
    --lattice 2000 \
    --obstacle_density 0.02 \
    --n_chains 16 \
    --n_runs   8 \
    --seed     43 \
    --outfile  prediction2_results_density_low.csv

# =============================================================================
# OBSTACLE DENSITY SENSITIVITY: rho = 0.10 (narrower tube)
# =============================================================================
run_sim "density_high" "$RESULTS_DIR/prediction2_density_high" \
    --Lvals   100,200,300,500 \
    --lattice 2000 \
    --obstacle_density 0.10 \
    --n_chains 16 \
    --n_runs   8 \
    --seed     44 \
    --outfile  prediction2_results_density_high.csv

# =============================================================================
# SEED REPRODUCIBILITY CHECK 1
# =============================================================================
run_sim "seed_check_1" "$RESULTS_DIR/prediction2_seed_check_1" \
    --Lvals   200,300,500 \
    --lattice 2000 \
    --obstacle_density 0.05 \
    --n_chains 16 \
    --n_runs   8 \
    --seed     100 \
    --outfile  prediction2_results_seed_check_1.csv

# =============================================================================
# SEED REPRODUCIBILITY CHECK 2
# =============================================================================
run_sim "seed_check_2" "$RESULTS_DIR/prediction2_seed_check_2" \
    --Lvals   200,300,500 \
    --lattice 2000 \
    --obstacle_density 0.05 \
    --n_chains 16 \
    --n_runs   8 \
    --seed     200 \
    --outfile  prediction2_results_seed_check_2.csv

# =============================================================================
# EXTENDED RUN: L = 1000-2000
# =============================================================================
#
#run_sim "extended" "$RESULTS_DIR/prediction2_extended" \
#    --Lvals   1000,1500,2000 \
#    --lattice 5000 \
#    --obstacle_density 0.05 \
#    --n_chains 8 \
#    --n_runs   4 \
#    --seed     42 \
#    --outfile  prediction2_results_extended.csv

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo "All runs complete: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
echo "CSV results:"
find "$RESULTS_DIR" -name "*.csv" | sort | while read -r f; do
    echo "  $f"
done
echo ""
echo "Log files:"
find "$LOGS_DIR" -name "*.log" | sort | while read -r f; do
    n=$(wc -l < "$f")
    echo "  $f  ($n lines)"
done
echo ""
echo "To check for errors:"
echo "  grep -l 'ERROR\|Traceback\|Error' $LOGS_DIR/*.log"
echo ""
