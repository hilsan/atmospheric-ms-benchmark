#!/bin/bash
#SBATCH --job-name=restrip_qcxms2
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --output=restrip_%A_%a.out
#SBATCH --error=restrip_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --gres=nvme:30

LIST_FILE=$1

if [ -z "$LIST_FILE" ]; then
    echo "Usage: sbatch --array=0-151 restrip_qcxms2_archives.sh <list_file>"
    exit 1
fi

ARCHIVE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST_FILE")

if [ -z "$ARCHIVE" ]; then
    echo "No archive for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

if [ ! -f "$ARCHIVE" ]; then
    echo "Archive not found: $ARCHIVE"
    exit 1
fi

echo "Processing: $ARCHIVE"
echo "Archive size: $(du -sh "$ARCHIVE" | cut -f1)"

TMPDIR="$LOCAL_SCRATCH/restrip_$$"
mkdir -p "$TMPDIR"

# Extract excluding the large ORCA output files
tar -xzf "$ARCHIVE" \
    --exclude="*/orca.out" \
    --exclude="*/geo.out" \
    --exclude="*/g98.out" \
    -C "$TMPDIR"

if [ $? -ne 0 ]; then
    echo "Extraction failed for $ARCHIVE"
    rm -rf "$TMPDIR"
    exit 1
fi

# Recompress using all available CPUs
NEW_ARCHIVE="$LOCAL_SCRATCH/qcxms2_auxiliary_new_$$.tar.gz"
tar -czf "$NEW_ARCHIVE" -C "$TMPDIR" .

if [ $? -ne 0 ]; then
    echo "Recompression failed"
    rm -rf "$TMPDIR" "$NEW_ARCHIVE"
    exit 1
fi

echo "Old size: $(du -sh "$ARCHIVE" | cut -f1)  New size: $(du -sh "$NEW_ARCHIVE" | cut -f1)"

# Replace original
mv "$NEW_ARCHIVE" "$ARCHIVE"

rm -rf "$TMPDIR"

echo "Done: $ARCHIVE"
