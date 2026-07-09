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

if [ -d "MS-run/TMPQCXMS" ]; then
    echo "TMPQCXMS already exists — GS-MD already done. Skipping to protect in-progress fragmentation."
    exit 0
fi

# Check for previous getieeab failure with Poisson IEE — retry with Gaussian if detected
if [ -f "MS-run/qcxms.out" ] && grep -q "internal error inside getieeab" "MS-run/qcxms.out"; then
    if ! grep -q "^gauss" "MS-run/qcxms.in" 2>/dev/null; then
        echo "Detected getieeab error with Poisson IEE — retrying with Gaussian IEE distribution."
        cd MS-run
        find . -mindepth 1 -not -name "xtbopt.sdf" -delete
        cat > qcxms.in << EOF
tmax 25
iseed 10
method ei
gauss
EOF
        qcxms -i xtbopt.sdf > qcxms.out 2>&1
        qcxms -i xtbopt.sdf > qcxms.out 2>&1
        exit $?
    fi
fi

# Gaussian run completed normally but TMPQCXMS missing — retry
if [ -f "MS-run/qcxms.in" ] && grep -q "^gauss" "MS-run/qcxms.in" 2>/dev/null; then
    if [ -f "MS-run/qcxms.out" ] && grep -q "normal termination" "MS-run/qcxms.out" && ! grep -q "internal error inside getieeab" "MS-run/qcxms.out"; then
        echo "Gaussian run completed but TMPQCXMS missing — rerunning both QCxMS passes."
        cd MS-run
        find . -mindepth 1 -not -name "xtbopt.sdf" -not -name "qcxms.in" -delete
        qcxms -i xtbopt.sdf > qcxms.out 2>&1
        qcxms -i xtbopt.sdf > qcxms.out 2>&1
        exit $?
    fi
fi

if [ ! -f "xtbopt.sdf" ]; then
    cp ../structure.sdf input.sdf
    xtb input.sdf --opt extreme > opt.out
    if [[ ! -f xtbopt.sdf ]]; then
        echo "ERROR: xtbopt.sdf not created!"
        exit 1
    fi
fi

mkdir -p MS-run
cd MS-run

cp ../xtbopt.sdf .

cat > qcxms.in << EOF
tmax 25
iseed 10
method ei
EOF

qcxms -i xtbopt.sdf > qcxms.out 2>&1
qcxms  -i xtbopt.sdf > qcxms.out 2>&1
