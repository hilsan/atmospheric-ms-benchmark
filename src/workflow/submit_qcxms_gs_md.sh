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
export PATH="/users/hsandstr/NEIMS/QCxMS2/qcxms2/bin:$PATH"

export KMP_AFFINITY=granularity=fine,compact,1,0

set -e
trap "echo 'Error occurred; exiting'; exit 1" ERR

DIR_NAME=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

exec > "$DIR_NAME/task.log" 2>&1

cd $DIR_NAME
mkdir -p QCxMS
cd QCxMS

# ------------------------------
# Copy input SDF
# ------------------------------
cp ../*.sdf input.sdf

# ------------------------------
# Detect atom count & check V2000
# ------------------------------
ATOM_LINE=$(sed -n '4p' input.sdf)

if [[ "$ATOM_LINE" != *"V2000"* ]]; then
    echo "Warning: V2000 not on line 4 — searching for it..."
    V2000_LINE=$(grep -n "V2000" input.sdf | head -n1 | cut -d: -f1)

    if [[ -z "$V2000_LINE" ]]; then
        echo "ERROR: No V2000 line found in SDF. Cannot proceed."
        exit 1
    else
        echo "Found V2000 at line $V2000_LINE (non-standard position). Proceeding."
    fi
fi

# Extract atom count from the counts line (first 4 columns)
ATOM_COUNT=$(awk 'NR==4 {print $1}' input.sdf)

echo "Detected $ATOM_COUNT atoms in SDF."


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

if (( ATOM_COUNT > 34 )); then
  qcxms -i xtbopt.sdf  > qcxms.out 2>&1
  qcxms -i xtbopt.sdf --unity -v  > qcxms.out 2>&1
   echo "UNITY not used (atom count $ATOM_COUNT > 34)"

    
else
    qcxms -i xtbopt.sdf  > qcxms.out 2>&1
    qcxms -i xtbopt.sdf  > qcxms.out 2>&1
    echo "UNITY not used (atom count $ATOM_COUNT <= 34)"
fi



