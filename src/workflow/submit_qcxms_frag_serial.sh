#!/bin/bash

# -------------------------------------------------------------------
#   Initial setup
# -------------------------------------------------------------------

d=$(pwd)
bname=$(basename "$d")

DIR="$d/TMPQCXMS"
if [ ! -d "$DIR" ]; then 
    echo "Run mpspred first!"
    exit 1
fi

echo "Starting parallel QCxMS run on $DIR"

cd "$DIR" || exit 1

# -------------------------------------------------------------------
#   Determine atom count (once per molecule)
# -------------------------------------------------------------------

# Pick a representative SDF (from any fragment)
REP_FRAG=$(find "$d" -type f -name "*.sdf" | head -n1)

if [ -z "$REP_FRAG" ]; then
    echo "No SDF file found in $DIR"
    exit 1
fi

# Check V2000 line
ATOM_LINE=$(sed -n '4p' "$REP_FRAG")

if [[ "$ATOM_LINE" != *"V2000"* ]]; then
    echo "Warning: V2000 not on line 4 — searching..."
    V2000_LINE=$(grep -n "V2000" "$REP_FRAG" | head -n1 | cut -d: -f1)
    if [[ -z "$V2000_LINE" ]]; then
        echo "ERROR: No V2000 line found in SDF. Cannot proceed."
        exit 1
    else
        echo "Found V2000 at line $V2000_LINE (non-standard position). Proceeding."
    fi
fi

# Extract atom count from counts line
ATOM_COUNT=$(awk 'NR==4 {print $1}' "$REP_FRAG")
echo "Detected atom count: $ATOM_COUNT"

# -------------------------------------------------------------------
#   Collect TMPQCXMS subdirectories
# -------------------------------------------------------------------

dir_list=()

for vz in "$DIR"/*/; do
    [ -d "$vz" ] || continue

    vz_name=$(basename "$vz")

    # Skip already finished jobs
    if [ -f "$vz/qcxms.out" ] && grep -q "normal termination of QCxMS" "$vz/qcxms.out"; then
        echo "Skipping $vz_name (already finished successfully)"
        continue
    fi

    # Clean previous outputs without touching SLURM files
    rm -f "$vz"/*.err "$vz"/*.out "$vz"/*.log

    dir_list+=("$vz_name")
done

TOTAL=${#dir_list[@]}
echo "TOTAL = $TOTAL directories to process"

if (( TOTAL == 0 )); then
    echo "All directories finished, nothing to submit."
    exit 0
fi

# -------------------------------------------------------------------
#   Determine SLURM array size
# -------------------------------------------------------------------

PER_TASK=300
ARRAY_SIZE=$(( (TOTAL + PER_TASK - 1) / PER_TASK ))

echo "Submitting array job with $ARRAY_SIZE tasks (up to $PER_TASK dirs per task)"

slurm_script="${bname}_multi.slurm"

# -------------------------------------------------------------------
#   Write the SLURM script
# -------------------------------------------------------------------
cat <<EOF > "$slurm_script"
#!/bin/bash -l

#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --job-name=qcxms_multi
#SBATCH --output=frag_%A_%a.out
#SBATCH --error=frag_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --time=20:00:00
#SBATCH --array=0-$((ARRAY_SIZE - 1))

set -x
export OMP_NUM_THREADS=1

DIR="$DIR"
PER_TASK=$PER_TASK
ATOM_COUNT=$ATOM_COUNT

# Directory list (embedded from parent script)
dir_list=(
EOF

# Append array contents safely
for item in "${dir_list[@]}"; do
    printf "\"%s\"\n" "$item" >> "$slurm_script"
done

# Finish the SLURM script with atom-count dependent QCxMS
cat <<'EOF' >> "$slurm_script"
)

START_INDEX=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END_INDEX=$(( START_INDEX + PER_TASK - 1 ))

# Clamp END_INDEX to array size
if (( END_INDEX >= ${#dir_list[@]} )); then
    END_INDEX=$(( ${#dir_list[@]} - 1 ))
fi

for ((i=START_INDEX; i<=END_INDEX; i++)); do
    vz="${dir_list[i]}"
    VZ_DIR="$DIR/$vz"

    if [ ! -d "$VZ_DIR" ]; then
        echo "Directory $VZ_DIR not found, skipping"
        continue
    fi

    cd "$VZ_DIR" || continue
    echo "Running QCxMS in $VZ_DIR (ATOM_COUNT=$ATOM_COUNT)"

    # Inside each fragment folder (VZ_DIR) before running QCxMS:
    echo "tmax 10" > qcxms.in
    echo "iseed 10" >> qcxms.in
    echo "method ei" >> qcxms.in
    
    if (( ATOM_COUNT > 34 )); then
        # Large molecule — run UNITY
        pqcxms --prod --unity  > qcxms.out 2>&1
        echo "UNITY used (atom count $ATOM_COUNT > 34)"
    else
        # Small molecule — skip UNITY
        pqcxms  --prod > qcxms.out 2>&1
        echo "UNITY not used (atom count $ATOM_COUNT <= 34)"
    fi


    touch ready
done
EOF

# -------------------------------------------------------------------
#   Submit SLURM array job
# -------------------------------------------------------------------
sbatch "$slurm_script"
echo "Submitted single array job of $ARRAY_SIZE tasks for $TOTAL directories"
