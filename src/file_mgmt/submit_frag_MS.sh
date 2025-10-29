#!/bin/bash

# Set working directory and project name
d=$(pwd)
bname=$(basename "$d")

# Directory where QCxMS inputs are located
DIR="$d/TMPQCXMS"
if [ ! -d "$DIR" ]; then 
    echo "Run mpspred first!"
    exit 1
fi

echo "Starting parallel qcxms run on"

# Clean previous outputs
rm -f qcxms.out T*/*slurm T*/*err T*/*out T*/*err

cd "$DIR" || exit 1

# Collect all directories to process
dir_list=()
for vz in *; do
    [ -d "$vz" ] || continue
    cd "$vz" || continue
    rm -f ready
    cd .. || continue
    dir_list+=("$vz")
done

TOTAL=${#dir_list[@]}
echo "TOTAL = $TOTAL directories"

# Determine max number of array tasks per SLURM policy (e.g., 300)
MAX_TASKS_PER_ARRAY=300
BATCH_SIZE=$(( (TOTAL + MAX_TASKS_PER_ARRAY - 1) / MAX_TASKS_PER_ARRAY ))  # ceil(TOTAL/300)

# Write one job array script per batch
for ((batch=0; batch<BATCH_SIZE; batch++)); do
    batch_start=$((batch * MAX_TASKS_PER_ARRAY))
    batch_end=$(( (batch + 1) * MAX_TASKS_PER_ARRAY - 1 ))
    if (( batch_end >= TOTAL )); then
        batch_end=$((TOTAL - 1))
    fi

    # Subset of directory names for this batch
    sublist=("${dir_list[@]:batch_start:batch_end - batch_start + 1}")
    array_size=${#sublist[@]}

    slurm_script="${bname}_array_batch${batch}.slurm"

    cat <<EOF > "$slurm_script"
#!/bin/bash -l

#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_array_batch${batch}
#SBATCH --output=frag_%A_%a.out
#SBATCH --error=frag_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000
#SBATCH --time=72:00:00
#SBATCH --array=0-$((array_size - 1))

set -x
export OMP_NUM_THREADS=4

cd "$DIR"

dir_list=(${sublist[@]})

vz=\${dir_list[\$SLURM_ARRAY_TASK_ID]}
VZ_DIR=\$DIR/\$vz

if [ ! -d "\$VZ_DIR" ]; then
    echo "Directory \$VZ_DIR not found"
    exit 1
fi

cd "\$VZ_DIR" || exit 1

qcxms --prod > qcxms.out 2>&1
touch ready
EOF

    # Submit each SLURM array job
    sbatch "$slurm_script"
    echo "Submitted batch $batch for directories $batch_start to $batch_end"
done
