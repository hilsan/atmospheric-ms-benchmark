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

    # Check if QCxMS2 has terminated normally in the qcxms2_gfn2 log file
    if grep -q "QCxMS2 terminated normally" "$dir/QCxMS2/qcxms2_gfn2"/*.log; then
        echo "QCxMS2 terminated normally in $dir. Skipping job submission."
        continue
    fi

    # If QCxMS2 didn't terminate normally, clean the qcxms2_gfn2 directory
    echo "QCxMS2 did not terminate normally in $dir. Cleaning directory and submitting the job..."
    rm -rf "$dir/QCxMS2/qcxms2_gfn2"/*

    # Submit the job using sbatch
    cd "$dir/QCxMS2" || exit 1
    sbatch $SCRIPTS_NEIMS/qcxms2_gfn2_4.sh

    # Return to the original directory
    cd "$initial_dir"

    echo "Job submitted for $dir after cleaning. Moving to the next directory."
done

# If all directories processed successfully
echo "All directories processed. No further issues found."

