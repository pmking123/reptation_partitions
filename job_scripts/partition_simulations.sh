#!/usr/bin/env bash
#SBATCH -o partition_simulations.%j.out
#SBATCH -e partition_simulations.%j.err
#SBATCH -J partition_simulations
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --partition=CPU

set -euo pipefail

WORKDIR="$HOME/rp_production/partition_simulations"
PROGRAM="$WORKDIR/partition_simulations.py"
LOG_FILE="$WORKDIR/partition_simulations_out_1.txt"

exec > "$LOG_FILE" 2>&1

START_TIME=$(date +%s)

echo "===================================================================="
echo "Node          : $(hostname)"
echo "Started       : $(date)"
echo "===================================================================="

echo "[$(date +%T)] Python starting"

python3 "$PROGRAM" 

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
