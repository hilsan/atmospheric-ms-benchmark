#!/bin/bash

# Initial setup
d=$(pwd)
bname=$(basename "$d")

DIR="$d/TMPQCXMS"
if [ ! -d "$DIR" ]; then 
    echo "Run mpspred first!"
    exit 1
fi

echo "Starting parallel qcxms run on"

# Clean logs
rm -f qcxms.out T*/*slurm T*/*err T*/*out T*/*err
cd "$DIR" || exit

# Collect subdirectories
dir_list=()
for vz in */; do
    cd "$vz" || continue
    rm -f ready
    cd ..
    dir_list+=("${vz%/}")
done

TOTAL=${#dir_list[@]}
echo "TOTAL = $TOTAL directories"

MAX_ARRAY=300
PER_TASK=$(( (TOTAL + MAX_ARRAY - 1) / MAX_ARRAY ))  # ceil division

# Write SLURM script
slurm_script="${bname}_multi.slurm"

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
#SBATCH --array=0-$((MAX_ARRAY - 1))

set -x
export OMP_NUM_THREADS=1

DIR="$DIR"
cd "\$DIR"

# All directories
dir_list=(${dir_list[@]})

PER_TASK=$PER_TASK
START_INDEX=\$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END_INDEX=\$(( START_INDEX + PER_TASK - 1 ))
if (( END_INDEX >= ${#dir_list[@]} )); then
    END_INDEX=\$(( ${#dir_list[@]} - 1 ))
fi

for ((i=START_INDEX; i<=END_INDEX; i++)); do
    vz=\${dir_list[i]}
    VZ_DIR="\$DIR/\$vz"
    if [ -d "\$VZ_DIR" ]; then
        echo "Processing \$VZ_DIR"
        cd "\$VZ_DIR" || continue
        qcxms --prod > qcxms.out 2>&1
        touch ready
        cd "\$DIR"
    fi
done
EOF

# Submit just one array job
sbatch "$slurm_script"
echo "Submitted single array job of 300 tasks for $TOTAL directories"
