#!/bin/bash
#SBATCH --job-name=mytest
#SBATCH --account=project_2006752
#SBATCH --partition=small #medium #gpusmall
#SBATCH --time=00:15:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=%j_%A_%a.out
##SBATCH --gres=gpu:a100:1,nvme:950

module load python-data/3.10-22.09

# Load Conda environment properly
#eval "$(conda shell.bash hook)"
#conda activate neims

# Define model directory
MODEL_DIR="/users/hsandstr/NEIMS/deep-molecular-massspec/"

# Zero-pad the task ID to create the directory name
DIR_NAME=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")

# Full path to the array job directory
ARRAY_JOB_DIR="$(pwd)/$DIR_NAME"

# Create required directories
mkdir -p "$ARRAY_JOB_DIR/xtb_opt"

# Check if the required input file exists
if [ ! -f "$ARRAY_JOB_DIR/QCxMS/xtbopt.sdf" ]; then
    echo "Error: Input file $ARRAY_JOB_DIR/QCxMS/xtbopt.sdf not found!"
    exit 1
fi

# Copy input file to xtb_opt directory
cp "$ARRAY_JOB_DIR/QCxMS/xtbopt.sdf" "$ARRAY_JOB_DIR/xtb_opt/"

# Change directory to xtb_opt
cd "$ARRAY_JOB_DIR/xtb_opt" || { echo "Error: Failed to enter directory $ARRAY_JOB_DIR/xtb_opt"; exit 1; }

# Run the Python script for spectra prediction
python "$MODEL_DIR/make_spectra_prediction.py" \
    --input_file="$ARRAY_JOB_DIR/xtb_opt/xtbopt.sdf" \
    --output_file="$ARRAY_JOB_DIR/xtb_opt/annotated.sdf" \
    --weights_dir="$MODEL_DIR/massspec_weights"


