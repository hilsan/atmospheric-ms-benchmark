#!/bin/bash

# Base directory
BASE_DIR="/scratch/project_2006752/hsandstr/Project/EIMS/Dataset_analysis/data/NEIMS/Franklin_TMS-v2"

# Loop over directories 0002 to 0068 (zero-padded to 4 digits)
for dir in $(seq -w 0016 0068); do
    RUN_DIR="$BASE_DIR/$dir/QCxMS/MS-run"

    echo "Processing directory $RUN_DIR"

    if [ ! -d "$RUN_DIR" ]; then
        echo "❌ Failed to enter directory $RUN_DIR"
        continue
    fi

    cd "$RUN_DIR" || continue

    # Run the submission script and capture job ID
    output=$(sh $SCRIPTS_NEIMS/submit_frag_parallell.sh)
    echo "$output"

    job_id=$(echo "$output" | grep -oP 'Submitted batch job \K[0-9]+')

    if [ -z "$job_id" ]; then
        echo "❌ Failed to get SLURM job ID for $dir, skipping..."
        continue
    fi

    echo "🕒 Waiting for SLURM job $job_id to finish..."

    # Poll SLURM queue every 10 minutes
    while squeue -j "$job_id" | grep -q "$job_id"; do
        echo "Job $job_id still running... checking again in 10 minutes."
        sleep 1200
    done

    echo "✅ Job $job_id finished for directory $dir."
done

echo "🎉 All submissions completed."

