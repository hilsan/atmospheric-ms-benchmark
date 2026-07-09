#!/bin/bash
#SBATCH --job-name=crest
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=crest_log.out
#SBATCH --error=crest.err

module purge
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
ulimit -s unlimited

export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so
export KMP_AFFINITY=granularity=fine,compact,1,0

echo "======================================"
echo "Running in: $(pwd)"
echo "Started: $(date)"
echo "======================================"
echo

if [ ! -f structure.sdf ]; then
    echo "structure.sdf not found! Exiting."
    exit 1
fi

############################################
# RESTART WATCHER
############################################
RESTART="$(pwd)/crest.restart"
BACKUP="$(pwd)/crest.restart.bak"
WATCH_LOG="$(pwd)/watch_restart.log"
WATCH_INTERVAL=300
MIN_SIZE=10240

watch_restart() {
    echo "Watcher started: $(date)" > "$WATCH_LOG"
    while true; do
        sleep $WATCH_INTERVAL
        if [ -f "$RESTART" ]; then
            SIZE1=$(stat -c%s "$RESTART" 2>/dev/null || echo 0)
            if (( SIZE1 < MIN_SIZE )); then
                echo "$(date): file too small (${SIZE1} bytes), skipping" >> "$WATCH_LOG"
                continue
            fi
            sleep 2
            SIZE2=$(stat -c%s "$RESTART" 2>/dev/null || echo 0)
            if [ "$SIZE1" -eq "$SIZE2" ]; then
                cp "$RESTART" "$BACKUP"
                echo "$(date): backup updated (${SIZE1} bytes)" >> "$WATCH_LOG"
            else
                echo "$(date): still being written, skipping" >> "$WATCH_LOG"
            fi
        fi
    done
}

watch_restart &
WATCHER_PID=$!
echo "Restart watcher launched (PID=$WATCHER_PID)"

############################################
# HELPER: run CREST
############################################
CREST_ARGS="-T ${SLURM_CPUS_PER_TASK} --chrg 1 --ewin 4 --rthr 0.2"

run_crest() {
    local input=$1
    local log=$2
    crest "$input" $CREST_ARGS > "$log" 2>&1
    return $?
}

crest_succeeded() {
    grep -q "CREST terminated normally." crest*.log 2>/dev/null
}

############################################
# MAIN LOGIC
############################################
if [ -f crest.restart ]; then
    SIZE=$(stat -c%s crest.restart 2>/dev/null || echo 0)
    echo "crest.restart found (${SIZE} bytes)"

    if (( SIZE >= MIN_SIZE )); then
        echo "Attempting restart from crest.restart..."
        run_crest crest.restart crest_restart.log
    else
        echo "crest.restart too small — likely corrupted."
    fi

    # Try backup if restart didn't succeed
    if ! crest_succeeded; then
        if [ -f "$BACKUP" ] && (( $(stat -c%s "$BACKUP") >= MIN_SIZE )); then
            echo "Trying backup restart file..."
            cp "$BACKUP" crest.restart
            run_crest crest.restart crest_restart.log
        fi
    fi

    # Fall back to fresh run if still not succeeded
    if ! crest_succeeded; then
        echo "Restart failed — falling back to fresh run from xtbopt.sdf."
        rm -f crest.restart crest.restart.bak

        if [ ! -f xtbopt.sdf ]; then
            echo "xtbopt.sdf not found — rerunning xTB optimization..."
            xtb structure.sdf --opt --gfn 2 --chrg 1 -P ${SLURM_CPUS_PER_TASK} > xtb_opt.out
            if [ $? -ne 0 ]; then
                kill $WATCHER_PID 2>/dev/null
                echo "xTB optimization failed. Exiting."
                exit 1
            fi
        fi

        run_crest xtbopt.sdf crest.log
    fi

else
    ############################################
    # Fresh run
    ############################################
    echo "No crest.restart found — starting fresh."
    echo "Optimizing structure with GFN2-xTB..."

    xtb structure.sdf --opt --gfn 2 --chrg 1 -P ${SLURM_CPUS_PER_TASK} > xtb_opt.out
    if [ $? -ne 0 ]; then
        kill $WATCHER_PID 2>/dev/null
        echo "xTB optimization failed. Exiting."
        exit 1
    fi

    echo "Running CREST conformer search..."
    run_crest xtbopt.sdf crest.log
fi

kill $WATCHER_PID 2>/dev/null
echo "Watcher stopped."

############################################
# CHECK NORMAL TERMINATION
############################################
if ! crest_succeeded; then
    echo "CREST did not terminate normally — keeping restart file."
    exit 1
fi

echo
echo "CREST finished successfully."

############################################
# CLEAN DIRECTORY
############################################
echo "Cleaning directory..."
find . -type f \
    ! -name "crest_best.xyz" \
    ! -name "structure.sdf" \
    ! -name "crest_conformers.xyz" \
    ! -name "xtbopt.sdf" \
    ! -name "xtb_opt.out" \
    ! -name "crest.log" \
    ! -name "crest_restart.log" \
    ! -name "crest_log.out" \
    ! -name "crest.err" \
    -delete

echo
echo "Finished: $(date)"
echo "======================================"