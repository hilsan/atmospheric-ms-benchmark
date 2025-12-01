#!/bin/bash
#SBATCH --job-name=qcxms2_w
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=qcxms2_w_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=336:00:00

# Unset unnecessary environment variables
export PATH=$(echo $PATH | sed -e 's|:/users/hsandstr/NEIMS/conda_env_neims/neims/bin||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin:||' -e 's|/users/hsandstr/NEIMS/conda_env_neims/neims/bin||')

# Load required modules
module purge
module load gcc/11.3.0 openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
module load biopythontools/11.3.0_3.10.6

# Set MKL threading layer to GNU to avoid conflicts
export MKL_THREADING_LAYER=GNU

# Ensure the directory exists
mkdir -p qcxms2_wb97x3c

# Ensure crest_best.xyz is copied into the directory as in.xyz
cp crest_best.xyz qcxms2_wb97x3c/in.xyz

# Change directory and run qcxms2
# Change directory and run qcxms2
cd qcxms2_wb97x3c
qcxms2 in.xyz  -edist gaussian  -T 32 -geolevel gfn2 -iplevel gfn2 -tslevel wb97x3c -ip2level wb97x3c -notsgeo > qcxms2_wb97x3c.log 2>&1

# List of specific files to keep
keep_files=("peaks.csv" "qcxms2_wb97x3c.log" "in.xyz")

# Make a temp dir for safe keeping
mkdir -p keep_temp

# Move specific files to the temp directory
for f in "${keep_files[@]}"; do
    [[ -e $f ]] && mv "$f" keep_temp/
done


# Create a tarball of everything else (excluding keep_temp)
tar --exclude=keep_temp -czf qcxms2_auxiliary.tar.gz ./*

# Delete everything else (except the tarball and keep_temp)
find . -mindepth 1 -maxdepth 1 ! -name 'qcxms2_auxiliary.tar.gz' ! -name 'keep_temp' -exec rm -rf {} +

# Restore the important files and directories
mv keep_temp/* . && rmdir keep_temp


