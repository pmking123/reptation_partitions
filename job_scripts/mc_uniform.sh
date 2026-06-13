#!/bin/bash
#SBATCH --job-name=mc_uniform
#SBATCH --array=0-11                 # 9 mode tasks + 3 RLD tasks
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=mcu_%A_%a.out
#SBATCH --error=mcu_%A_%a.err
#SBATCH --partition=CPU

# ---------------------------------------------------------------------------
# SLURM job array for mc_uniform_fast.py.
#
# Tasks 0-8:  mode/mean convergence study (MODE_L_VALS)
#   Index  L
#   -----  -----
#     0      20
#     1      50
#     2     100
#     3     200
#     4     500
#     5    1000
#     6    2000
#     7    5000
#     8   10000   <-- longest task; allow 48h
#
# Tasks 9-11: run-length distribution (RLD_L_VALS)
#   Index  L
#   -----  -----
#     9     100
#    10     200
#    11     500
#
# Results are written to mode_<L>.json
# and rld_<L>.json.
#
# After all tasks complete, run:
#     python3 aggregate.py --mode uniform --indir "$HOME/rp_production/mc_uniform"
# ---------------------------------------------------------------------------

WORKDIR="$HOME/rp_production/mc_uniform"
PROGRAM="$WORKDIR/mc_uniform_fast.py"

START_TIME=$(date +%s)

echo "============================================"
echo "Host     : $(hostname)"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Started  : $(date)"
echo "============================================"

python3 "$PROGRAM" \
    --task     "${SLURM_ARRAY_TASK_ID}" \
    --outdir   "$WORKDIR" \
    --seed     42 \
    --n-steps  600000

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo "===================================================================="
echo "Finished      : $(date)"
echo "Exit code     : ${EXIT_CODE}"
echo "Elapsed       : ${ELAPSED_MIN}m ${ELAPSED_SEC}s  (${ELAPSED}s total)"
echo "===================================================================="

exit ${EXIT_CODE}
