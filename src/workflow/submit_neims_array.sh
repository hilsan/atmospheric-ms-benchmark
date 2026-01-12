#!/bin/bash
#SBATCH --job-name=NEIMS
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=%j_%A_%a.out
#SBATCH --array=0-68

# 1. Clean the environment inherited from your Jupyter session
module purge
unset PYTHONPATH
unset PYTHONHOME
export MPLBACKEND="Agg"

# 2. Set the path to your Tykky installation
export PATH="/scratch/project_2006752/hsandstr/bin:$PATH" 

MODEL_DIR=/users/hsandstr/NEIMS/deep-molecular-massspec/
DIR_NAME=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
ARRAY_JOB_DIR="$(pwd)/$DIR_NAME"

# 3. Run the prediction
if [ -d "$ARRAY_JOB_DIR" ]; then
    python "$MODEL_DIR/make_spectra_prediction.py" \
        --input_file="$ARRAY_JOB_DIR/structure.sdf" \
        --output_file="$ARRAY_JOB_DIR/annotated.sdf" \
        --weights_dir="$MODEL_DIR/massspec_weights"
else
    echo "Directory $ARRAY_JOB_DIR does not exist."
fi