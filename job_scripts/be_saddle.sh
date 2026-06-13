#!/usr/bin/env bash

WORKDIR="$HOME/rp_production/be_saddle"
PROGRAM="$WORKDIR/be_saddle.py"
LOG_FILE="$WORKDIR/bes_out_1.txt"

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
