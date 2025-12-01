#!/bin/bash

# Loop over directories 0000 to 0068
for dir in {0000..0068}; do
    # Check if QCxMS2 directory exists
    if [[ -d "$dir/QCxMS2" ]]; then
        # Look for "CREST terminated normally" in any G*out file
        if ! grep -q "CREST terminated normally" "$dir/QCxMS2/G"*out 2>/dev/null; then
            echo "CREST did not terminate normally in $dir. Submitting batch_crest.sh..."
            
            # Enter the directory
            cd "$dir/QCxMS2" || { echo "Failed to enter $dir/QCxMS2"; continue; }
	    cp ../smiles.sdf .            
            # Submit the job
            sbatch "$SCRIPTS_NEIMS/batch_crest.sh"
            
            # Exit back to the original directory
            cd - > /dev/null
        else
            echo "CREST completed normally in $dir. Skipping."
        fi
    else
        echo "Directory $dir/QCxMS2 does not exist. Skipping."
    fi
done

