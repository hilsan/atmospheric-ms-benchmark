#!/bin/bash

initial_dir=$(pwd)

for dir in {0000..0068}; do
    echo "Checking directory $dir"

    if [ -d "$dir/QCxMS2" ]; then
        echo "$dir/QCxMS2 exists."

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

    shopt -s nullglob
    log_files=("$dir/QCxMS2/qcxms2_wb97x3c/"*.log)
    shopt -u nullglob

    if [ ${#log_files[@]} -eq 0 ]; then
        echo "No log files found in $dir/QCxMS2/qcxms2_wb97x3c. Skipping."
        continue
    fi

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
        rm -rf "$dir/QCxMS2/qcxms2_wb97x3c/"*

        cd "$dir/QCxMS2" || exit 1
        sbatch "$SCRIPTS_NEIMS/qcxms2_wB97x3c_gaussian.sh"
        cd "$initial_dir"

        echo "Job submitted for $dir after cleaning. Moving to the next directory."
    else
        echo "No matching error found or job terminated normally in $dir. Skipping submission."
    fi
done
