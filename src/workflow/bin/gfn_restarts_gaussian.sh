#!/bin/bash

# Save the current directory
initial_dir=$(pwd)

# Loop over directories 0000 to 0068
for dir in {0000..0068}; do
    echo "Checking directory $dir"

    # Check if the QCxMS2 directory exists
    if [ -d "$dir/QCxMS2" ]; then
        echo "$dir/QCxMS2 exists."

        # Check if CREST has finished successfully in $dir/QCxMS2
        if ! grep -q "CREST terminated normally" "$dir/QCxMS2"/*.out; then
            echo "CREST did not terminate normally in $dir/QCxMS2. Skipping directory."
            continue
        else
            echo "CREST terminated normally in $dir/QCxMS2."
        fi
    else
        echo "$dir/QCxMS2 does not exist or cannot be accessed. Skipping."
        continue
    fi

    # Check log files in qcxms2_gfn2 directory
    log_files=("$dir/QCxMS2/qcxms2_gfn2"/*.log)

    # If no log files found, skip
    if [ ! -e "${log_files[0]}" ]; then
        echo "No log files found in $dir/QCxMS2/qcxms2_gfn2. Skipping."
        continue
    fi

    # Flag to determine if job should be submitted
    submit_job=false

    for log_file in "${log_files[@]}"; do
        if ! grep -q "QCxMS2 terminated normally" "$log_file"; then
            if grep -q "internal error inside getieeab" "$log_file"; then
                submit_job=true
                break
            fi
        fi
    done

    if [ "$submit_job" = true ]; then
        echo "QCxMS2 failed with internal error in $dir. Cleaning and submitting the job..."
        rm -rf "$dir/QCxMS2/qcxms2_gfn2"/*

        cd "$dir/QCxMS2" || exit 1
        sbatch "$SCRIPTS_NEIMS/qcxms2_gfn2_gaussian.sh"
        cd "$initial_dir"

        echo "Job submitted for $dir after cleaning. Moving to the next directory."
    else
        echo "No matching error found or job terminated normally in $dir. Skipping submission."
    fi
done

