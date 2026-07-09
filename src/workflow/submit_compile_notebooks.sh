#!/bin/bash
#SBATCH --job-name=compile_notebooks
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=%j_compile_notebooks.out
#SBATCH --error=%j_compile_notebooks.err

set -euo pipefail

NB_DIR="/scratch/project_2006752/hsandstr/Project/atmospheric-ms-benchmark/notebook"
JUPYTER="/users/hsandstr/NEIMS/conda_env_neims/neims/bin/jupyter"

cd "$NB_DIR"

echo "=== Running 2_Franklin_TMS.ipynb ==="
"$JUPYTER" nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=10800 \
    "$NB_DIR/2_Franklin_TMS.ipynb"
echo "Done: 2_Franklin_TMS.ipynb"

echo "=== Running 3_Franklin.ipynb ==="
"$JUPYTER" nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=10800 \
    "$NB_DIR/3_Franklin.ipynb"
echo "Done: 3_Franklin.ipynb"
