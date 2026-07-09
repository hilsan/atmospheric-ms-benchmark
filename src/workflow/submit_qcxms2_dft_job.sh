#!/bin/bash
#SBATCH --job-name=qcxms2_dft
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=qcxms2_dft_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --gres=nvme:400
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hildaagnesolivia@gmail.com

# QCxMS2 with recommended mixed-level settings:
#   GFN2-xTB  : geometry optimizations (-geolevel) and IP prescreening (-iplevel)
#   wB97X-3c  : barrier/energy calculations and IP refinement (-tslevel, -ip2level)

set -e
set -o pipefail

BASE_DIR=$1

if [ -z "$BASE_DIR" ]; then
    echo "Usage: sbatch submit_qcxms2_dft_job.sh <BASE_DIR>"
    exit 1
fi

WORKDIR=$(pwd)
FOLDER=$(basename "$WORKDIR")
OUT_DIR="qcxms2_wb97x3c"
NVME_WORKDIR="$LOCAL_SCRATCH/qcxms2_work"
NVME_OUT="$NVME_WORKDIR/$OUT_DIR"

echo "Processing folder: ${WORKDIR}"

if [ ! -f crest_best.xyz ]; then
    echo "crest_best.xyz not found in ${WORKDIR}"
    exit 1
fi

module purge
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
module load biopythontools/11.3.0_3.10.6

export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Trap: copy qcxms2.log and EDIST_GAUSSIAN back to Lustre on any exit so
# getieeab detection and permanent-fail detection work on the next submission.
cleanup() {
    local dest="$WORKDIR/$OUT_DIR"
    mkdir -p "$dest"
    for f in qcxms2.log EDIST_GAUSSIAN; do
        [ -f "$NVME_OUT/$f" ] && cp "$NVME_OUT/$f" "$dest/" 2>/dev/null || true
    done
}
trap cleanup EXIT

# -----------------------------
# Check for prior run state on Lustre
# -----------------------------
USE_GAUSSIAN=false

if [ -d "$WORKDIR/$OUT_DIR" ]; then
    # Already marked as permanently failed
    if [ -f "$WORKDIR/$OUT_DIR/QCXMS2_PERMANENT_FAIL" ]; then
        echo "Molecule marked as permanently failed — skipping."
        trap - EXIT
        exit 0
    fi

    if grep -q "QCxMS2 terminated normally" "$WORKDIR/$OUT_DIR/qcxms2.log" 2>/dev/null; then
        echo "Completed run found in $OUT_DIR — skipping."
        trap - EXIT
        exit 0
    elif grep -q "fragment generation failed" "$WORKDIR/$OUT_DIR/qcxms2.log" 2>/dev/null; then
        echo "Fragment generation failed (msreact error) — marking as permanently failed."
        echo "fragment generation failed detected in qcxms2.log on $(date)" \
            > "$WORKDIR/$OUT_DIR/QCXMS2_PERMANENT_FAIL"
        trap - EXIT
        exit 0
    elif grep -q "internal error inside getieeab" "$WORKDIR/$OUT_DIR/qcxms2.log" 2>/dev/null; then
        if [ -f "$WORKDIR/$OUT_DIR/EDIST_GAUSSIAN" ]; then
            echo "getieeab error persists even with Gaussian edist — molecule unrecoverable. Skipping."
            trap - EXIT
            exit 0
        fi
        echo "IEE distribution error detected — retrying with Gaussian edist."
        rm -rf "$WORKDIR/$OUT_DIR"
        USE_GAUSSIAN=true
    else
        echo "Incomplete run found in $OUT_DIR — removing stale Lustre data and restarting on NVMe."
        rm -rf "$WORKDIR/$OUT_DIR"
    fi
fi

# -----------------------------
# Set up NVMe workspace
# -----------------------------
mkdir -p "$NVME_OUT"
cp "$WORKDIR/crest_best.xyz" "$NVME_OUT/in.xyz"
cd "$NVME_OUT"

if [ "$USE_GAUSSIAN" = true ]; then
    echo "This molecule failed with Poisson IEE distribution." > EDIST_GAUSSIAN
    echo "Rerun with -edist gaussian on $(date)"              >> EDIST_GAUSSIAN
    echo "Folder: ${WORKDIR}"                                 >> EDIST_GAUSSIAN
    qcxms2 in.xyz \
        -T ${SLURM_CPUS_PER_TASK} \
        -geolevel gfn2 \
        -iplevel gfn2 \
        -tslevel wb97x3c \
        -ip2level wb97x3c \
        -notsgeo \
        -edist gaussian > qcxms2.log 2>&1
else
    qcxms2 in.xyz \
        -T ${SLURM_CPUS_PER_TASK} \
        -geolevel gfn2 \
        -iplevel gfn2 \
        -tslevel wb97x3c \
        -ip2level wb97x3c \
        -notsgeo > qcxms2.log 2>&1
fi

# Check for msreact failure (exit code 0 but no peaks.csv)
if grep -q "fragment generation failed" qcxms2.log 2>/dev/null; then
    echo "Fragment generation failed (msreact error) — marking as permanently failed."
    mkdir -p "$WORKDIR/$OUT_DIR"
    echo "fragment generation failed detected in qcxms2.log on $(date)" \
        > "$WORKDIR/$OUT_DIR/QCXMS2_PERMANENT_FAIL"
    cp qcxms2.log "$WORKDIR/$OUT_DIR/"
    trap - EXIT
    exit 0
fi

# -----------------------------
# Archive auxiliary files on NVMe, then copy results to Lustre
# -----------------------------
KEEP_FILES=("peaks.csv" "qcxms2.log" "in.xyz" "*.in" "EDIST_GAUSSIAN")

mkdir -p keep_temp
shopt -s nullglob
for pattern in "${KEEP_FILES[@]}"; do
    for file in $pattern; do
        [[ -e "$file" ]] && mv -- "$file" keep_temp/
    done
done

find . -not -path "./keep_temp/*" \( -name "orca.out" -o -name "geo.out" -o -name "g98.out" \) -delete
tar --exclude='./keep_temp' -czf qcxms2_auxiliary.tar.gz ./*

find . -mindepth 1 -maxdepth 1 \
    ! -name 'qcxms2_auxiliary.tar.gz' \
    ! -name 'keep_temp' \
    -exec rm -rf {} +

mv keep_temp/* ./
rmdir keep_temp

# Copy final results to Lustre
trap - EXIT
mkdir -p "$WORKDIR/$OUT_DIR"
cp -r "$NVME_OUT/"* "$WORKDIR/$OUT_DIR/"

echo "QCxMS2 DFT run completed. Results copied to Lustre."
echo "Finished ${FOLDER} with wB97X-3c"
