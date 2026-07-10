#!/bin/bash
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=compress_tmp
#SBATCH --output=compress_tmp_%A_%a.out
#SBATCH --error=compress_tmp_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00

LIST=/scratch/project_2006752/hsandstr/Project/atmospheric-ms-benchmark/src/workflow/compress_tmpqcxms_list.txt

TARGET_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [ -z "$TARGET_DIR" ]; then
    echo "No entry for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

MSRUN=$(dirname "$TARGET_DIR")

if [ ! -d "$TARGET_DIR" ]; then
    echo "Already compressed or missing: $TARGET_DIR"
    exit 0
fi

if [ -f "$MSRUN/TMPQCXMS.tar.gz" ]; then
    echo "Archive exists — testing integrity: $MSRUN/TMPQCXMS.tar.gz"
    if tar -tzf "$MSRUN/TMPQCXMS.tar.gz" > /dev/null 2>&1; then
        echo "Archive OK — removing raw dir: $TARGET_DIR"
        rm -rf "$TARGET_DIR"
        exit 0
    else
        echo "Archive corrupt — removing and re-compressing"
        rm -f "$MSRUN/TMPQCXMS.tar.gz"
    fi
fi

echo "Compressing: $TARGET_DIR"
cd "$MSRUN"
tar -czf TMPQCXMS.tar.gz TMPQCXMS/
if [ $? -eq 0 ]; then
    rm -rf TMPQCXMS/
    echo "Done: $TARGET_DIR"
else
    echo "tar failed for $TARGET_DIR — leaving uncompressed"
    rm -f TMPQCXMS.tar.gz
    exit 1
fi
