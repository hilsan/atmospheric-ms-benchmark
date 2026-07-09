#!/bin/bash
# Quick status check for QCxMS2 DFT jobs.
# Usage: bash check_dft_jobs.sh [dataset]   (default: all datasets)

BASE=/scratch/project_2006752/hsandstr/Project/atmospheric-ms-benchmark/data/simulation_results
DATASETS="${1:-franklin franklin_tms ucb_globes_tracers}"

echo "=== Queue status ==="
pending=$(squeue -u hsandstr --format="%.10i %.15j %.8T" --noheader | grep "qcxms2_dft" | grep -c "PENDING" || true)
running=$(squeue -u hsandstr --format="%.10i %.15j %.8T" --noheader | grep "qcxms2_dft" | grep -c "RUNNING" || true)
echo "  PENDING: $pending   RUNNING: $running"

echo ""
echo "=== Recently completed/failed DFT jobs (last 12h) ==="
sacct --starttime=$(date -d '12 hours ago' '+%Y-%m-%dT%H:%M:%S') \
      --format=JobID,JobName,State,ExitCode,Elapsed,WorkDir \
      --noheader 2>/dev/null \
    | grep "qcxms2_dft" \
    | grep -v "\.batch\|\.extern" \
    | head -20 || echo "  (none)"

echo ""
echo "=== DFT completion per dataset ==="
for dataset in $DATASETS; do
    dft_dir="$BASE/$dataset/QCxMS2_dft"
    [ -d "$dft_dir" ] || continue

    total=0; done_count=0; perm_fail=0; in_progress=0; not_started=0
    for mol_dir in "$dft_dir"/*/; do
        [ -d "$mol_dir" ] || continue
        total=$((total + 1))
        wdir="$mol_dir/qcxms2_wb97x3c"
        if [ -f "$wdir/QCXMS2_PERMANENT_FAIL" ]; then
            perm_fail=$((perm_fail + 1))
        elif [ -f "$wdir/peaks.csv" ]; then
            done_count=$((done_count + 1))
        elif [ -f "$wdir/qcxms2.log" ]; then
            in_progress=$((in_progress + 1))
        else
            not_started=$((not_started + 1))
        fi
    done
    echo "  $dataset: $total total | $done_count done | $in_progress in-progress/failed | $not_started not started | $perm_fail perm-fail"
done

echo ""
echo "=== Recent SLURM output snippets (running jobs) ==="
squeue -u hsandstr --format="%i %Z" --noheader 2>/dev/null \
    | grep "QCxMS2_dft" \
    | while read jobid wdir; do
        out=$(ls "$wdir"/qcxms2_dft_${jobid}.out 2>/dev/null | head -1)
        [ -f "$out" ] && tail -3 "$out" | sed "s/^/  [$jobid] /"
    done | head -30
