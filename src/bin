#!/bin/bash
#SBATCH --job-name=mytest
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_GS_sampling
#SBATCH --output=GS_%A_%a.out          # Standard output log (with job and task ID)
#SBATCH --error=GS_%A_%a.err           # Standard error log (with job and task ID)
#SBATCH --ntasks=1                     # Number of tasks per job
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per task
#SBATCH --time=10:00:00                # Maximum runtime per task

module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0

ulimit -s unlimited

# Set MKL and OpenMP threading
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so

export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')
export PATH="/users/hsandstr/NEIMS/QCxMS2/qcxms2/bin:$PATH"

# Improve CPU core affinity for better performance
export KMP_AFFINITY=granularity=fine,compact,1,0


# Error handling
set -e
trap "echo 'Error occurred; exiting'; exit 1" ERR


DIR_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)


# Logging
exec > "$DIR_NAME/task.log" 2>&1

# Activate environment
#conda activate neims

# Navigate to task directory
cd $DIR_NAME
mkdir -p QCxMS
cd QCxMS

# Copy input files
cp ../*.sdf .

# Run xtb
xtb s*sdf --opt extreme > opt.out

# Setup and run qcxms
mkdir -p MS-run
cd MS-run
cp ../xtbopt.sdf .
echo "tmax 10" >> qcxms.in
echo "iseed 10" >> qcxms.in
qcxms -i xtbopt.sdf
qcxms -i xtbopt.sdf

# Post-processing or cleanup if needed
