#!/bin/bash
#SBATCH --job-name=mytest
#SBATCH --account=project_2006752
#SBATCH --partition=small 
#SBATCH --time=00:15:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=%j_%A_%a.out
##SBATCH --gres=gpu:a100:1,nvme:950

#module load python-data/3.10-22.09
#source activate neims

# Define model directory
MODEL_DIR=/users/hsandstr/NEIMS/deep-molecular-massspec/

# Zero-pad the task ID to create the directory name
DIR_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

# Full path to the array job directory
ARRAY_JOB_DIR=$(pwd)/$DIR_NAME

# Check if the directory exists
if [ -d "$ARRAY_JOB_DIR" ]; then
    # Execute the command
    python "$MODEL_DIR/make_spectra_prediction.py" \
        --input_file="$ARRAY_JOB_DIR/smiles.sdf" \
        --output_file="$ARRAY_JOB_DIR/annotated.sdf" \
        --weights_dir="$MODEL_DIR/massspec_weights"
else
    echo "Directory $ARRAY_JOB_DIR does not exist."
fi
