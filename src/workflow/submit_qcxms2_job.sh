#!/bin/bash
#SBATCH --job-name=qcxms2_run
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=qcxms2_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=336:00:00

set -e  # Exit on error
set -o pipefail

# -----------------------------
# Modules and environment
# -----------------------------
module purge
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
module load biopythontools/11.3.0_3.10.6

export MKL_THREADING_LAYER=GNU
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# -----------------------------
# Parse flags
# -----------------------------
usage() {
    echo "Usage: $0 -i INPUT_FILE -g GEOLEVEL -t TSLEVEL -p IPLEVEL -d DIST_METHOD"
    echo "  -i INPUT_FILE   CREST best xyz input"
    echo "  -g GEOLEVEL     Geometry level (gfn2, 6, 7, ...)"
    echo "  -t TSLEVEL      Transition state level"
    echo "  -p IPLEVEL      IP2 level"
    echo "  -d DIST_METHOD  Distance/energy method (gfn2, gaussian)"
    exit 1
}

while getopts "i:g:t:p:d:" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        g) GEOLEVEL="$OPTARG" ;;
        t) TSLEVEL="$OPTARG" ;;
        p) IPLEVEL="$OPTARG" ;;
        d) DIST_METHOD="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$INPUT_FILE" || -z "$GEOLEVEL" || -z "$TSLEVEL" || -z "$IPLEVEL" || -z "$DIST_METHOD" ]] && usage

# -----------------------------
# Setup directories
# -----------------------------
OUT_DIR="qcxms2_${DIST_METHOD}"
mkdir -p "$OUT_DIR"
cp "$INPUT_FILE" "$OUT_DIR/in.xyz"

cd "$OUT_DIR"

# -----------------------------
# Construct QCxMS command
# -----------------------------
QCXMS_CMD="qcxms in.xyz -T ${SLURM_CPUS_PER_TASK} -geolevel $GEOLEVEL -tslevel $TSLEVEL -iplevel $IPLEVEL -edist $DIST_METHOD -notsgeo"
echo "Running command: $QCXMS_CMD" | tee qcxms2_command.log

# -----------------------------
# Run QCxMS
# -----------------------------
$QCXMS_CMD > qcxms2.log 2>&1

# -----------------------------
# Keep essential files
# -----------------------------
KEEP_FILES=("peaks.csv" "qcxms2.log" "in.xyz" "*.in" "qcxms2_command.log")
mkdir -p keep_temp

for f in "${KEEP_FILES[@]}"; do
    for file in $f; do
        [[ -e $file ]] && mv "$file" keep_temp/
    done
done

# Tar and clean auxiliary files
tar --exclude=keep_temp -czf qcxms2_auxiliary.tar.gz ./*
find . -mindepth 1 -maxdepth 1 ! -name 'qcxms2_auxiliary.tar.gz' ! -name 'keep_temp' -exec rm -rf {} +

mv keep_temp/* . && rmdir keep_temp

echo "QCxMS run completed. Essential files kept; auxiliary files archived."

