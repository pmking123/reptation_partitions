#!/bin/bash
#SBATCH --job-name=mc
#SBATCH --array=0-5                  # one task per L in SCALING_L_VALS
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=mc_%A_%a.out
#SBATCH --error=mc_%A_%a.err
#SBATCH --partition=CPU

# ---------------------------------------------------------------------------
# SLURM job array for mc_fast.py (weighted-measure scaling study).
#
# Array index  L value
# -----------  -------
#     0          20
#     1          50
#     2         100
#     3         200
#     4         300
#     5         500
#
# Results are written to scaling_<L>.json.
# After all tasks complete, run:
#     python3 aggregate.py --mode weighted --indir "$HOME/rp_production/mc"
# ---------------------------------------------------------------------------

START_TIME=$(date +%s)

WORKDIR="$HOME/rp_production/mc"
PROGRAM="$WORKDIR/mc_fast.py"

echo "=================================="
echo "Host     : $(hostname)"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Started  : $(date)"
echo "=================================="

python3 "$PROGRAM" \
    --task   "${SLURM_ARRAY_TASK_ID}" \
    --outdir "$WORKDIR" \
    --seed   42

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
