#!/bin/bash

# Check argument
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <submit_script_path>"
    exit 1
fi

SUBMIT_SCRIPT="$1"

# Verify the submit script exists
if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    echo "Error: submit script '$SUBMIT_SCRIPT' not found."
    exit 1
fi

dirs=()

for dir in 0*; do
    [[ -d "$dir" ]] || continue

    msrun="$dir/GS-opt/MS-run"
    outfiles=("$msrun"/*.out)

    failed=false

    # No MS-run directory or no .out files → failed
    if [[ ! -d "$msrun" ]] || [[ ! -e "${outfiles[0]}" ]]; then
        failed=true
    else
        # No successful termination in any .out file → failed
        if ! grep -q "normal termination of QCxMS" "${outfiles[@]}"; then
            failed=true
        fi
    fi

    if [[ "$failed" == true ]]; then
        echo "FAILED: $dir — cleaning GS-opt and MS-run"

        # Safety: clear files
        [[ -d "$dir/GS-opt" ]] && rm -rf "$dir/GS-opt"/*
        [[ -d "$msrun" ]] && rm -rf "$msrun"/*

        dirs+=("$dir")
    fi
done

# Submit array job
if (( ${#dirs[@]} > 0 )); then
    IFS=,
    sbatch --array="${dirs[*]}" "$SUBMIT_SCRIPT"
    unset IFS
else
    echo "No failed QCxMS runs found."
fi
