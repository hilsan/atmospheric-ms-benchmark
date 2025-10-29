#!/bin/bash
#SBATCH --job-name=qcxms2_GS_sampling
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=GS_%A_%a.out          # Standard output log (with job and task ID)
##SBATCH --error=GS_%A_%a.err           # Standard error log (with job and task ID)
#SBATCH --ntasks=1                     # Number of tasks per job
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=32G                      # Memory per task
#SBATCH --time=148:00:00                # Maximum runtime per task

module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0

ulimit -s unlimited

# Set MKL and OpenMP threading
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so

export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')
export PATH="/users/hsandstr/NEIMS/QCxMS2/qcxms2/bin:$PATH"

# Improve CPU core affinity for better performance
export KMP_AFFINITY=granularity=fine,compact,1,0

# Debugging
ldd $(which crest) | grep mkl

# Try to restart the CREST calculation
crest -restart  -T 16 --chrg 1
# Check if the restart was successful
if [ $? -ne 0 ]; then
    # If the restart failed, run the original CREST command
    echo "Restart failed, running initial CREST calculation..."
    crest smiles.sdf -T 16 --chrg 1
else
    echo "Restart successful, no need to run the initial calculation."
fi

