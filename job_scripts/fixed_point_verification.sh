#!/usr/bin/env bash
#SBATCH -o fpv.%j.out
#SBATCH -e fpv.%j.err
#SBATCH -J fpv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=CPU

START_TIME=$(date +%s)

echo "===================================================================="
echo "Job ID        : ${SLURM_JOB_ID}"
echo "Node          : $(hostname)"
echo "CPUs          : ${SLURM_CPUS_PER_TASK}"
echo "Started       : $(date)"
echo "===================================================================="

echo "[$(date +%T)] Python starting"

WORKDIR="$HOME/rp_production/fixed_point_verification"
PROGRAM="$WORKDIR/fixed_point_verification.py"
OUTPUT_FILE="$WORKDIR/fpv_out_1.txt"

python3 "$PROGRAM" > "$OUTPUT_FILE"

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
