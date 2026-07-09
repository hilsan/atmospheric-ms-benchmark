#!/bin/bash
# check_qcxms_runs.sh
# Run getres and assess trajectory completion for each molecule directory.
# Skips molecules where plotms.res already passes quality gates.

set -euo pipefail

BASE_DIR="${1:-$(pwd)}"
SUMMARY_FILE="${2:-$BASE_DIR/analysis_decision.txt}"
MIN_SUCCESS_PERCENT=90

> "$SUMMARY_FILE"
echo "Checking MS-run directories in $BASE_DIR"
echo "Minimum success threshold: $MIN_SUCCESS_PERCENT%"
echo "Summary will be written to $SUMMARY_FILE"
echo

for i in $(seq -w 0000 0060); do

    MSRUN_DIR="$BASE_DIR/$i/GS-opt/MS-run"
    TMP_DIR="$MSRUN_DIR/TMPQCXMS"

    if [[ ! -d "$MSRUN_DIR" ]]; then
        echo "$i : MS-run directory missing -> EXCLUDE" | tee -a "$SUMMARY_FILE"
        continue
    fi

    # --- Count TMP.* dirs once, reuse below ---
    TOTAL=0
    if [[ -d "$TMP_DIR" ]]; then
        TOTAL=$(find "$TMP_DIR" -maxdepth 1 -type d -name "TMP.*" 2>/dev/null | wc -l)
    fi

    # --- Skip if plotms.res already passes quality gates ---
    PLOTMS_RES="$MSRUN_DIR/plotms.res"
    if [[ -f "$PLOTMS_RES" && "$TOTAL" -gt 0 ]]; then
        PLOTMS_RUNS=$(grep -oE '^[0-9]+\s+successfull runs' "$PLOTMS_RES" 2>/dev/null | grep -oE '^[0-9]+' || echo "")
        PLOTMS_BASE=$(grep -oE 'Theoretical counts in 100 %\s+signal:\s+[0-9]+' "$PLOTMS_RES" 2>/dev/null | grep -oE '[0-9]+$' || echo "")
        if [[ -n "$PLOTMS_RUNS" && -n "$PLOTMS_BASE" ]]; then
            THRESHOLD=$(awk "BEGIN {printf \"%d\", $TOTAL * $MIN_SUCCESS_PERCENT / 100}")
            if [[ "$PLOTMS_RUNS" -ge "$THRESHOLD" && "$PLOTMS_BASE" -ge 95 ]]; then
                echo "$i : already done ($PLOTMS_RUNS runs, base peak $PLOTMS_BASE) -> INCLUDE" | tee -a "$SUMMARY_FILE"
                continue
            fi
        fi
    fi

    if [[ "$TOTAL" -eq 0 ]]; then
        echo "$i : 0 TMP folders -> EXCLUDE" | tee -a "$SUMMARY_FILE"
        continue
    fi

    echo "Processing $i ..."

    # --- Skip getres if tmpqcxms.res exists and previous run met the success threshold ---
    TMPRES="$MSRUN_DIR/tmpqcxms.res"
    NEED_GETRES=true
    if [[ -f "$TMPRES" && -f "$MSRUN_DIR/success_stats.txt" ]]; then
        PREV_SUCCESS=$(grep -oE '^[0-9]+' "$MSRUN_DIR/success_stats.txt" 2>/dev/null || echo 0)
        THRESHOLD=$(awk "BEGIN {printf \"%d\", $TOTAL * $MIN_SUCCESS_PERCENT / 100}")
        if [[ "$PREV_SUCCESS" -ge "$THRESHOLD" ]]; then
            echo "$i : getres already good ($PREV_SUCCESS/$TOTAL), skipping"
            NEED_GETRES=false
        fi
    fi

    if [[ "$NEED_GETRES" == true ]]; then
        if [[ -d "$TMP_DIR" ]]; then
            rm -rf "$TMP_DIR"/frag*out "$TMP_DIR"/frag*err "$TMP_DIR"/MS*slurm* 2>/dev/null || true
        fi

        rm -f "$TMPRES" "$MSRUN_DIR/success_stats.txt"
        (
            cd "$MSRUN_DIR"
            if command -v getres &>/dev/null; then
                getres > success_stats.txt || echo "getres failed" > success_stats.txt
            else
                echo "getres command not found" > success_stats.txt
            fi
        )
    fi

    # --- Count successes from getres output ---
    SUCCESS=0
    if [[ -f "$MSRUN_DIR/success_stats.txt" ]]; then
        SUCCESS=$(grep -oE '^[0-9]+' "$MSRUN_DIR/success_stats.txt" 2>/dev/null || echo 0)
    fi

    PERCENT=$(awk "BEGIN {printf \"%.2f\", ($SUCCESS/$TOTAL)*100}")

    if awk "BEGIN {exit !($PERCENT >= $MIN_SUCCESS_PERCENT)}"; then
        DECISION="INCLUDE"
    else
        DECISION="EXCLUDE"
    fi

    echo "$i : $SUCCESS/$TOTAL ($PERCENT%) -> $DECISION" | tee -a "$SUMMARY_FILE"

done

echo
echo "Finished. Summary saved in $SUMMARY_FILE"
