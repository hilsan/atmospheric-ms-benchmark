#!/bin/bash
#SBATCH --job-name=qcxms2_scaling_test
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --output=GS_cpu%j.out
#SBATCH --error=GS_cpu%j.err
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00

# Get CPU count from input argument
N_THREADS=$1

# Load required modules
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0

# Set unlimited stack size
ulimit -s unlimited

# Set threading environment
export OMP_NUM_THREADS=$N_THREADS
export MKL_NUM_THREADS=$N_THREADS
export LD_PRELOAD=${MKLROOT}/lib/intel64/libmkl_rt.so

# Clean up PATH to avoid conda conflicts
export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')

# Add QCxMS2 crest path
export PATH="/users/hsandstr/NEIMS/QCxMS2/qcxms2/bin:$PATH"

# Improve thread affinity
export KMP_AFFINITY=granularity=fine,compact,1,0

# Create subdirectory and move there
OUTDIR="cpu_${N_THREADS}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Copy input
cp ../smiles.sdf .

echo "Running with $N_THREADS CPUs in $OUTDIR"
echo "Start time: $(date)"

# Optional debug info
ldd $(which crest) | grep mkl

# Run crest
crest -restart -T $N_THREADS --chrg 1
if [ $? -ne 0 ]; then
    echo "Restart failed, starting fresh..."
    crest smiles.sdf -T $N_THREADS --chrg 1
else
    echo "Restart successful."
fi

echo "End time: $(date)"
