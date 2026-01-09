#!/bin/bash
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_GS_sampling
#SBATCH --output=GS_%A_%a.out
#SBATCH --error=GS_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=10:00:00

module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0

ulimit -s unlimited

export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so

export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')

export KMP_AFFINITY=granularity=fine,compact,1,0

set -e
trap "echo 'Error occurred; exiting'; exit 1" ERR

DIR_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)



cd $DIR_NAME
mkdir -p GS-opt
cd GS-opt

# ------------------------------
# Copy input SDF
# ------------------------------
cp ../*.sdf input.sdf

# ------------------------------
# Run xtb optimization
# ------------------------------
xtb input.sdf --opt extreme > opt.out

# Ensure output exists
if [[ ! -f xtbopt.sdf ]]; then
    echo "ERROR: xtbopt.sdf not created!"
    exit 1
fi

# ------------------------------
# Set up QCxMS MS-run
# ------------------------------
mkdir -p MS-run
cd MS-run

cp ../xtbopt.sdf .

# ------------------------------
# Run QCxMS2
# ------------------------------
echo "tmax 25" > qcxms.in
echo "iseed 10" >> qcxms.in


qcxms -i xtbopt.sdf > qcxms.out 2>&1
qcxms -i xtbopt.sdf  > qcxms.out 2>&1




