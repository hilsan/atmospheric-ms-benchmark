#!/bin/bash
#SBATCH --job-name=cfmid
#SBATCH --account=project_2006752
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=cfmid_%j.out
#SBATCH --error=cfmid_%j.err

set -euo pipefail

CFMID_DIR=$1
RESULTS_DIR=$CFMID_DIR/ml_results.log
IDX_FILE=$CFMID_DIR/idx_smiles.txt
SIF=/users/hsandstr/cfmid_3.0.2.1.sif

echo "Started: $(date)"

# Suppress TYKKY/Puhti bind paths that cause container creation failures on
# compute nodes where the TYKKY mount point does not exist.
# Strategy: (1) capture any /PUHTI_TYKKY_* path before unsetting so we can
# disable it explicitly via --no-mount /path (a separate code path from the
# global --no-mount bind-paths flag); (2) unset all bind-path env var variants.
TYKKY_EXTRA=""
for _p in $(printf '%s\n' "${APPTAINER_BINDPATH:-}" "${APPTAINER_BIND:-}" \
            "${SINGULARITY_BINDPATH:-}" "${SINGULARITY_BIND:-}" \
            | tr ':,' '\n' | grep '/PUHTI_TYKKY_'); do
    TYKKY_EXTRA="$TYKKY_EXTRA --no-mount $_p"
done
unset APPTAINER_BINDPATH APPTAINER_BIND SINGULARITY_BINDPATH SINGULARITY_BIND

# Run CFM-EI (non-isotope)
{ time apptainer exec \
  --no-mount bind-paths $TYKKY_EXTRA \
  --bind $CFMID_DIR:/cfmid/public \
  --pwd /opt/cfm/bin \
  $SIF \
  ./cfm-predict \
  /cfmid/public/idx_smiles.txt \
  0.001 \
  /trained_models_cfmid2.0/ei_ms_models/ei_nn_iso_new/param_output.log \
  /trained_models_cfmid2.0/ei_ms_models/ei_nn_iso_new/param_config.txt \
  0 \
  /cfmid/public/ml_results.log \
  0 \
  0; } 2>&1

echo "CFM-EI (non-iso) finished: $(date)"

for logfile in "$RESULTS_DIR"/*.log; do
    mol=$(basename "$logfile" .log)
    mkdir -p "$CFMID_DIR/$mol"
    mv "$logfile" "$CFMID_DIR/$mol/${mol}.log"
done
rm -rf "$RESULTS_DIR"

# Run CFM-EI (isotope)
{ time apptainer exec \
  --no-mount bind-paths $TYKKY_EXTRA \
  --bind $CFMID_DIR:/cfmid/public \
  --pwd /opt/cfm/bin \
  $SIF \
  ./cfm-predict \
  /cfmid/public/idx_smiles.txt \
  0.001 \
  /trained_models_cfmid2.0/ei_ms_models/ei_nn_iso_new/param_output.log \
  /trained_models_cfmid2.0/ei_ms_models/ei_nn_iso_new/param_config.txt \
  0 \
  /cfmid/public/ml_results.log \
  1 \
  0; } 2>&1

echo "CFM-EI (iso) finished: $(date)"

for logfile in "$RESULTS_DIR"/*.log; do
    mol=$(basename "$logfile" .log)
    mkdir -p "$CFMID_DIR/$mol"
    mv "$logfile" "$CFMID_DIR/$mol/${mol}_iso.log"
done

rm -rf "$RESULTS_DIR"
rm -f "$IDX_FILE"

echo "Finished: $(date)"