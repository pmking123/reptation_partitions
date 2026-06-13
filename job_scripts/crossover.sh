#!/usr/bin/env bash
#SBATCH -o crossover.%j.out
#SBATCH -e crossover.%j.err
#SBATCH -J crossover
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --partition=CPU

set -euo pipefail

WORKDIR="$HOME/rp_production/crossover"
PROGRAM="$WORKDIR/crossover.py"

START_TIME=$(date +%s)

echo "===================================================================="
echo "Node          : $(hostname)"
echo "Started       : $(date)"
echo "===================================================================="

echo "[$(date +%T)] Python starting"

python3 "$PROGRAM" \
	--Ne 10,20,30,50,75,85,100,200 \
	--s 1 \
	--Lmax 20000 \
        --prefix s1 > "$WORKDIR/crossover_out_1.txt" 2>&1
python3 "$PROGRAM" \
        --Ne 10,20,30,50,75,85,100,200 \
	--s 2 \
	--Lmax 20000 \
	--prefix s2 > "$WORKDIR/crossover_out_2.txt" 2>&1
python3 "$PROGRAM" \
	--Ne 10,20,30,50,75,85,100,200 \
	--s 3 \
	--Lmax 20000 \
	--prefix s3 > "$WORKDIR/crossover_out_3.txt" 2>&1

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
