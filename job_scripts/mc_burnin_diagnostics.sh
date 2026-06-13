#!/bin/bash
#SBATCH --job-name=mc_burnin_diagnostics
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=mc_burnin_diagnostics_%j.out
#SBATCH --error=mc_burnin_diagnostics_%j.err
#SBATCH --partition=CPU

# ---------------------------------------------------------------------------

START_TIME=$(date +%s)

WORKDIR="$HOME/rp_production/mc_burnin_diagnostics"
PROGRAM="$WORKDIR/burnin_diagnostics.py"

echo "=================================="
echo "Host     : $(hostname)"
echo "Started  : $(date)"
echo "=================================="

python3 "$PROGRAM" \
    --L 500 \
    --burnin 100000 \
    --chains 4 \
    --outdir "$WORKDIR"

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo "================================================================"
echo "Finished   : $(date)"
echo "Exit code  : ${EXIT_CODE}"
echo "Elapsed    : ${ELAPSED_MIN}m ${ELAPSED_SEC}s (${ELAPSED}s total)"
echo "================================================================"

exit ${EXIT_CODE}
