#!/bin/bash
#SBATCH --job-name=rld
#SBATCH --output=rld_%j.out
#SBATCH --error=rld_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=CPU

# -------------------------------------------------------
# Run-length distribution analysis
# Single node, 48 cores (AMD EPYC 9454)
# Wall time: 15 min is generous; expect ~2-3 min
# -------------------------------------------------------

START_TIME=$(date +%s)

echo "======================================================"
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Node         : $(hostname)"
echo "CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "Started      : $(date)"
echo "======================================================"

WORKDIR="$HOME/rp_production/run_length_distribution"
PROGRAM="$WORKDIR/rld_fast.py"
LOG_FILE="$WORKDIR/rld_out_1.txt"

exec > "$LOG_FILE" 2>&1

echo "[$(date +%T)] Python starting"

python3 "$PROGRAM"

EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo "======================================================"
echo "Finished     : $(date)"
echo "Exit code    : ${EXIT_CODE}"
echo "Elapsed      : ${ELAPSED_MIN}m ${ELAPSED_SEC}s  (${ELAPSED}s total)"
echo "======================================================"

exit ${EXIT_CODE}
