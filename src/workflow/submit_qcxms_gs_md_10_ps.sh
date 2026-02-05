#!/bin/bash
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_GS_sampling
#SBATCH --output=GS_%A_%a.out
#SBATCH --error=GS_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=20:00:00

module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0

ulimit -s unlimited

# --- THREAD CONTROL (CRITICAL) ---
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export KMP_AFFINITY=granularity=fine,compact,1,0

export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so

# Clean PATH from NEIMS
export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' \
                               -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')

set -e
trap "echo 'Error occurred; exiting'; exit 1" ERR

DIR_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
cd $DIR_NAME

mkdir -p GS-opt
cd GS-opt

cp ../structure.sdf input.sdf

xtb input.sdf --opt extreme > opt.out

if [[ ! -f xtbopt.sdf ]]; then
    echo "ERROR: xtbopt.sdf not created!"
    exit 1
fi

mkdir -p MS-run
cd MS-run

cp ../xtbopt.sdf .

cat > qcxms.in << EOF
tmax 10
iseed 10
method ei
EOF

qcxms -i xtbopt.sdf > qcxms.out 2>&1
