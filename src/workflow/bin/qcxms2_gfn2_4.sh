#!/bin/bash
#SBATCH --job-name=qcxms2_gfn2
#SBATCH --account=project_2006752
#SBATCH --partition=longrun
#SBATCH --output=qcxms2_gfn2_%j.out
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
mkdir -p qcxms2_gfn2

# Copy the input file into the directory
cp crest_best.xyz qcxms2_gfn2/in.xyz

# Change directory and run qcxms2
cd qcxms2_gfn2
qcxms2 in.xyz -T 32 -geolevel gfn2 -iplevel gfn2 -tslevel gfn2 -ip2level gfn2 -notsgeo > qcxms2_gfn2.log 2>&1

# List of specific files to keep
keep_files=("peaks.csv" "qcxms2_gfn2.log")

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

