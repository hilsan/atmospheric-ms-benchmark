#!/bin/bash

# -------------------------------------------------------------------
#   Initial setup
# -------------------------------------------------------------------

d=$(pwd)
bname=$(basename "$d")

DIR="$d/TMPQCXMS"
if [ ! -d "$DIR" ]; then 
    echo "Run mpspred first!"
    exit 1
fi

echo "Starting parallel qcxms run on $DIR"

cd "$DIR" || exit 1

# -------------------------------------------------------------------
#   Collect TMPQCXMS subdirectories
# -------------------------------------------------------------------

dir_list=()

for vz in "$DIR"/*/; do
    [ -d "$vz" ] || continue

    vz_name=$(basename "$vz")

    # Check if already finished
    if [ -f "$vz/qcxms.out" ] && grep -q "normal termination of QCxMS" "$vz/qcxms.out"; then
        echo "Skipping $vz_name (already finished successfully)"
        continue
    fi

    # Clean previous outputs without touching SLURM files
    rm -f "$vz"/*.err "$vz"/*.out "$vz"/*.log

    dir_list+=("$vz_name")
done

TOTAL=${#dir_list[@]}
echo "TOTAL = $TOTAL directories to process"

if (( TOTAL == 0 )); then
    echo "All directories finished, nothing to submit."
    exit 0
fi

# -------------------------------------------------------------------
#   Determine SLURM array size
# -------------------------------------------------------------------

PER_TASK=300
ARRAY_SIZE=$(( (TOTAL + PER_TASK - 1) / PER_TASK ))

echo "Submitting array job with $ARRAY_SIZE tasks (up to $PER_TASK dirs per task)"

slurm_script="${bname}_multi.slurm"

# -------------------------------------------------------------------
#   Write the SLURM script
# -------------------------------------------------------------------
cat <<EOF > "$slurm_script"
#!/bin/bash -l

#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_multi
#SBATCH --output=frag_%A_%a.out
#SBATCH --error=frag_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --time=12:00:00
#SBATCH --array=0-$((ARRAY_SIZE - 1))

set -x
export OMP_NUM_THREADS=1

DIR="$DIR"
PER_TASK=$PER_TASK

# Directory list (embedded from parent script)
dir_list=(
EOF

# Append array contents safely
for item in "${dir_list[@]}"; do
    printf "\"%s\"\n" "$item" >> "$slurm_script"
done

# Finish the SLURM script
cat <<'EOF' >> "$slurm_script"
)

START_INDEX=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END_INDEX=$(( START_INDEX + PER_TASK - 1 ))

# Clamp END_INDEX to array size
if (( END_INDEX >= ${#dir_list[@]} )); then
    END_INDEX=$(( ${#dir_list[@]} - 1 ))
fi

for ((i=START_INDEX; i<=END_INDEX; i++)); do
    vz="${dir_list[i]}"
    VZ_DIR="$DIR/$vz"

    if [ ! -d "$VZ_DIR" ]; then
        echo "Directory $VZ_DIR not found, skipping"
        continue
    fi

    cd "$VZ_DIR" || continue

    echo "Running qcxms in $VZ_DIR"
    qcxms --prod > qcxms.out 2>&1
    touch ready
done
EOF

# -------------------------------------------------------------------
#   Submit SLURM array job
# -------------------------------------------------------------------
sbatch "$slurm_script"
echo "Submitted single array job of $ARRAY_SIZE tasks for $TOTAL directories"
