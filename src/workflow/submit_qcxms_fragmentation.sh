#!/bin/bash

# Default parameters
CPUS=1
MEM=8G
TIME=12:00:00
MAX_ARRAY=300

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --max-array) MAX_ARRAY="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

DIR=$(pwd)
BNAME=$(basename "$DIR")

# Collect all subdirectories
DIR_LIST=()
for d in */; do
  DIR_LIST+=("${d%/}")
done

TOTAL=${#DIR_LIST[@]}
PER_TASK=$(( (TOTAL + MAX_ARRAY - 1) / MAX_ARRAY ))

SLURM_SCRIPT="${BNAME}_multi.slurm"

# Generate SLURM script
cat <<EOF > "$SLURM_SCRIPT"
#!/bin/bash -l
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_multi
#SBATCH --output=frag_%A_%a.out
#SBATCH --error=frag_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem-per-cpu=$MEM
#SBATCH --time=$TIME
#SBATCH --array=0-$((MAX_ARRAY - 1))

set -x
export OMP_NUM_THREADS=$CPUS

DIR="$DIR"
cd "\$DIR"

DIR_LIST=(${DIR_LIST[@]})

START_INDEX=\$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END_INDEX=\$(( START_INDEX + PER_TASK - 1 ))
if (( END_INDEX >= ${#DIR_LIST[@]} )); then
  END_INDEX=\$(( ${#DIR_LIST[@]} - 1 ))
fi

for ((i=START_INDEX; i<=END_INDEX; i++)); do
  SUBDIR="\${DIR_LIST[i]}"
  echo "Processing \$SUBDIR" | tee -a frag_array_task_\$SLURM_ARRAY_TASK_ID.log
  cd "\$SUBDIR" || continue
  qcxms --prod > qcxms.out 2>&1
  touch ready
  cd "\$DIR"
done
EOF

# Submit the array
sbatch "$SLURM_SCRIPT"
echo "Submitted SLURM array job ($SLURM_SCRIPT) for $TOTAL directories"

