#!/bin/bash
#SBATCH --job-name=qcxms2_gfn2
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=qcxms2_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --gres=nvme:100
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=hildaagnesolivia@gmail.com

set -e
set -o pipefail

# -----------------------------
# Arguments
# -----------------------------
BASE_DIR=$1
METHOD=$2

if [ -z "$BASE_DIR" ] || [ -z "$METHOD" ]; then
    echo "Usage: sbatch submit_qcxms2_job.sh <BASE_DIR> <METHOD>"
    exit 1
fi

WORKDIR=$(pwd)
FOLDER=$(basename "$WORKDIR")
OUT_DIR="qcxms2_${METHOD}"
NVME_WORKDIR="$LOCAL_SCRATCH/qcxms2_work"
NVME_OUT="$NVME_WORKDIR/$OUT_DIR"

echo "Processing folder: ${WORKDIR}"

if [ ! -f crest_best.xyz ]; then
    echo "crest_best.xyz not found in ${WORKDIR}"
    exit 1
fi

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
# Modules
# -----------------------------
module purge
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
module load biopythontools/11.3.0_3.10.6

export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# -----------------------------
# Define levels
# -----------------------------
if [ "$METHOD" = "gfn2" ]; then
    GEOLEVEL="gfn2"
    TSLEVEL="gfn2"
    IPLEVEL="gfn2"
    IP2LEVEL="gfn2"
elif [ "$METHOD" = "wb97x3c" ]; then
    GEOLEVEL="gfn2"
    TSLEVEL="wb97x3c"
    IPLEVEL="gfn2"
    IP2LEVEL="wb97x3c"
else
    echo "Unknown method: $METHOD"
    exit 1
fi

# -----------------------------
# Check for prior run state on Lustre
# -----------------------------
EDIST_FLAG="EDIST_GAUSSIAN"
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
        if [ -f "$WORKDIR/$OUT_DIR/$EDIST_FLAG" ]; then
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

# Leave a note if switching to Gaussian (written before run so trap captures it on failure)
if [ "$USE_GAUSSIAN" = true ]; then
    echo "This molecule failed with Poisson IEE distribution." > "$EDIST_FLAG"
    echo "Rerun with -edist gaussian on $(date)"              >> "$EDIST_FLAG"
    echo "Folder: ${WORKDIR}"                                 >> "$EDIST_FLAG"
fi

# -----------------------------
# Run QCxMS2 on NVMe
# -----------------------------
if [ "$USE_GAUSSIAN" = true ]; then
    QCXMS_CMD="qcxms2 in.xyz -T ${SLURM_CPUS_PER_TASK} \
        -geolevel $GEOLEVEL \
        -tslevel $TSLEVEL \
        -iplevel $IPLEVEL \
        -ip2level $IP2LEVEL \
        -edist gaussian"
else
    QCXMS_CMD="qcxms2 in.xyz -T ${SLURM_CPUS_PER_TASK} \
        -geolevel $GEOLEVEL \
        -tslevel $TSLEVEL \
        -iplevel $IPLEVEL \
        -ip2level $IP2LEVEL"
fi

echo "Running: $QCXMS_CMD"
$QCXMS_CMD > qcxms2.log 2>&1

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
KEEP_FILES=("peaks.csv" "qcxms2.log" "in.xyz" "*.in" "qcxms2_command.log" "$EDIST_FLAG")

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

echo "QCxMS2 run completed. Results copied to Lustre."
echo "Finished ${FOLDER} with ${METHOD}"
