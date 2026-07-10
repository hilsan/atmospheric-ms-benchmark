# Atmospheric MS Benchmark

**By Hilda Sandström, 2025–2026.**

Benchmark of EI mass spectrum prediction methods for atmospheric organic compounds,
comparing QCxMS, QCxMS2, NEIMS, and CFM-ID against experimental NIST spectra.
Two molecular datasets are used: a plain set (`franklin`) and a TMS-derivatized set (`franklin_tms`).

---

## License

The code in this repository is released under the **GNU General Public License v3** (see `LICENSE`).

Data, figures, and documentation are released under **CC BY 4.0**
([Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)).
You are free to share and adapt the material for any purpose, provided appropriate credit is given.

---

## What this project does

1. Prepares two molecular datasets (raw and TMS-derivatized).
2. Submits and manages EI-MS simulations on a SLURM HPC cluster (Puhti, CSC):
   - **QCxMS** — semi-empirical trajectory-based EI fragmentation (GFN1-xTB)
   - **QCxMS2** — second-generation QCxMS with CREST conformer search (GFN2-xTB)
   - **NEIMS** — neural network predictor
   - **CFM-ID** — competitive fragmentation modelling
3. Processes raw simulation outputs into normalised spectra.
4. Compares simulated spectra against experimental NIST spectra using cosine, weighted-dot,
   Tanimoto, precision, and recall metrics.
5. Produces publication-quality plots and summary tables.
6. Tracks computational resource usage (CPU, memory, wall time) per method.

---

## Requirements

### Software (HPC environment — Puhti/CSC)
| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.10+ | All analysis scripts |
| QCxMS | 6.x | EI trajectory simulations |
| PlotMS | 6.2.0 | Compiles QCxMS trajectories into spectra (`~/PlotMS.v.6.2.0/plotms`) |
| getres | — | Checks QCxMS trajectory completion |
| QCxMS2 | 1.0.0 | Second-generation EI simulations |
| CREST | 3.x | Conformer search for QCxMS2 |
| xtb | — | GFN1/GFN2 tight-binding calculations |
| NEIMS | — | Neural network spectrum predictor (conda env at `~/NEIMS/conda_env_neims/neims/`) |
| CFM-ID | — | Competitive fragmentation modelling (at `~/cfmid/`) |
| SLURM | — | Job scheduling on HPC cluster |

### Python packages
Install via the NEIMS conda environment or `module load python-data/3.10-22.09`:
`pandas`, `numpy`, `matplotlib`, `seaborn`, `rdkit`, `scipy`, `ipython`

### Environment variables (set in `~/.bashrc`)
```bash
export DATA_NEIMS="/scratch/.../data/simulation_results/franklin_tms/"
export SCRIPTS_NEIMS="/scratch/.../src"
export PYTHONPATH="/scratch/.../src/:$PYTHONPATH"
```

---

## Repository structure

```
atmospheric-ms-benchmark/
├── data/
│   ├── raw/                        # Raw input data
│   │   └── franklin/goamazon.csv   # Source dataset (SMILES + metadata)
│   ├── processed/                  # Pre-processed datasets
│   │   ├── franklin_tms/           # TMS-derivatized SMILES CSV
│   │   └── franklin/               # Non-derivatized SMILES CSV
│   └── simulation_results/         # All simulation outputs (large — not in git)
│       ├── franklin_tms/           # TMS dataset results
│       │   ├── QCxMS_10_ps/        # QCxMS, 10 ps MD, per-molecule dirs 0000–0060
│       │   ├── QCxMS_25_ps/        # QCxMS, 25 ps MD
│       │   ├── QCxMS_10_ps_iee03/  # QCxMS, 10 ps, custom IEE
│       │   ├── QCxMS2/             # QCxMS2, per-molecule dirs with CREST output
│       │   ├── NEIMS/              # NEIMS, per-molecule annotated.sdf
│       │   ├── CFMID/              # CFM-ID, per-molecule log files
│       │   ├── EXP/                # Experimental NIST spectra (.MSPEC / exp.msp)
│       │   ├── memory_benchmark/   # 9-molecule subset for resource benchmarking
│       │   └── results/            # Per-molecule comparison CSVs (auto-generated)
│       └── franklin/               # Same structure, non-derivatized dataset
│
├── notebook/
│   ├── 2_Franklin_TMS.ipynb        # Main pipeline: TMS dataset (edit CONFIG cell only)
│   ├── 3_Franklin.ipynb            # Main pipeline: plain dataset (edit CONFIG cell only)
│   ├── 4_grid_search_qcxms.ipynb   # QCxMS parameter grid search (MD time, IEE)
│   ├── 5_resource_benchmark.ipynb  # CPU/memory/wall-time analysis across methods
│   ├── benchmark_template.ipynb    # Reusable template for new datasets
│   └── visualization_for_presentation.ipynb  # Publication plots and summary tables
│
├── src/
│   ├── analysis/
│   │   ├── compare_spectra.py      # Similarity metrics: cosine, weighted-dot, Tanimoto, precision, recall, entropy
│   │   ├── diagnose_spectra.py     # Coverage report: which molecules have spectra per method
│   │   ├── patch_entropy.py        # Back-fills Entropy_Similarity in comparison CSVs (requires Python 3.8+)
│   │   └── run_comparison.py       # Runs compare_spectra.py for all peak types via subprocess
│   ├── processing/
│   │   ├── bin_scale_intensities.py          # Bins m/z to integer, normalises to 0–1000, filters peaks
│   │   ├── make_directories_and_sdfs.py      # Creates per-molecule directories and 3D SDF files from SMILES
│   │   ├── make_TMS_derivative_251125_v1.py  # Derivatizes OH/COOH/NH groups with TMS; validates against reference
│   │   ├── make_inchlkey_sdf_for_nist.py     # Generates InChIKey/CAS for NIST lookup; sets up EXP folders
│   │   ├── process_spectra_batch.py          # Orchestrates bin_scale_intensities.py for all molecules/datasets
│   │   ├── remove_duplicate_SMILES_entries.py # Removes duplicate molecules by canonical SMILES
│   │   ├── NIST2MSP.py                       # Converts NIST web-text format to MSP (standalone utility)
│   │   └── copy_unique_folders.py            # Copies molecule folders using duplicate mapping (standalone utility)
│   ├── utils/
│   │   ├── parse_resource_usage.py  # Parses SLURM sacct + log files for CPU/memory/wall-time tables
│   │   ├── plotting.py              # Spectrum loading (CSV/MSP/SDF), binning, normalisation helpers
│   │   ├── processing.py            # SMILES/SDF handling, 3D embedding via RDKit
│   │   ├── check_qcxms_runs.sh      # Bash: runs getres for each molecule, writes analysis_decision.txt
│   │   ├── compare_smiles.py        # Detects duplicate SMILES across files (standalone utility)
│   │   └── select_molecule_by_size.py # Selects small/medium/large molecules by atom count (standalone utility)
│   ├── visualization/
│   │   ├── plot_results.py          # Master plotting: histograms, heatmaps, summary tables from comparison CSVs
│   │   ├── plot_mirrored_spectra.py # Mirror plot: experimental (inverted) vs simulated spectrum for a single molecule
│   │   └── plot_stacked_MS.py       # Stacked spectra plot comparing multiple methods for a single molecule
│   ├── workflow/
│   │   ├── job_submission.py              # Python helpers: submit SLURM arrays, run PlotMS, check CREST status
│   │   ├── submit_qcxms_gs_md_10_ps.sh    # SLURM: QCxMS ground-state MD (10 ps)
│   │   ├── submit_qcxms_gs_md_25_ps.sh    # SLURM: QCxMS ground-state MD (25 ps)
│   │   ├── submit_qcxms_frag_serial_no_unity.sh  # SLURM: QCxMS fragmentation (serial, no unity)
│   │   ├── submit_batch_crest.sh          # SLURM: CREST conformer search for QCxMS2
│   │   ├── submit_qcxms2_job.sh           # SLURM: QCxMS2 run (checks CREST done; cleans up on re-run)
│   │   ├── submit_neims_array.sh          # SLURM: NEIMS array job
│   │   └── submit_compile_notebooks.sh    # SLURM: run notebooks headlessly via jupyter nbconvert
│   ├── features/
│   │   ├── generate_mw.py           # Calculates molecular weight and atom counts from SMILES via RDKit
│   │   └── generate_simpol_groups.py # Counts SIMPOL functional groups using SMARTS patterns
│   ├── aprl_ssp/                    # APRL-SSP substructure search (functional group statistics)
│   │   ├── util.py                  # SMARTS pattern matching engine
│   │   └── userdef.py               # Custom ring/aromatic/nitrophenol detectors
│   └── depricated/                  # Archived older versions (not used)
│
├── reports/                         # Output figures and summary CSVs
└── LICENSE                          # GNU GPL v3
```

---

## Running the benchmark (notebook workflow)

All active workflow is driven by two notebooks: **`2_Franklin_TMS.ipynb`** and **`3_Franklin.ipynb`**.
Edit only the `CONFIG` cell at the top; all other cells are parameter-free.

### Step-by-step

#### 1. Data preparation (sections 1–6, `COMPILE_ONLY = False`)
Set `COMPILE_ONLY = False` and run sections 1–6 to:
- Deduplicate SMILES and apply TMS derivatization (section 1)
- Create molecule directories and submit QCxMS GS-MD jobs (section 2)
- Submit QCxMS2 CREST + fragmentation jobs (section 3)
- Submit NEIMS array job (section 4)
- Submit CFM-ID job (section 5)
- Set up experimental spectra EXP folders for NIST lookup (section 6)

After submitting, wait for all HPC jobs to finish before continuing.

#### 2. QCxMS post-processing (section 2, `COMPILE_ONLY = False`)

Once fragmentation jobs are done, re-run the **check / PlotMS cell** (section 2):
- `check_qcxms_runs.sh` runs `getres` for each molecule and writes `analysis_decision.txt`
  (molecules with < 90% completed trajectories are marked EXCLUDE)
- `run_plotms_for_included` runs PlotMS for all INCLUDE molecules

**PlotMS quality gates:**

| Criterion | Threshold |
|-----------|-----------|
| Successful trajectories | ≥ 90 % of total |
| Base peak counts | ≥ 95 counts |

Console messages:
- `[SKIP] already done` — passes both criteria, skipped
- `[SKIP] X/N (Z%) trajectories complete` — fragmentation incomplete, skipped
- `[SKIP] no TMPQCXMS trajectories found` — fragmentation never ran
- `[REDO] prior run below quality` — re-runs PlotMS if more trajectories completed
- `Done [OK]` / `Done [LOW QUALITY]` — result; LOW QUALITY molecules are excluded from comparison

#### 3. Spectrum processing and comparison (`COMPILE_ONLY = True`, sections 7–9)

Set `COMPILE_ONLY = True` and run sections 7–9:
- **Section 7** — `process_spectra_batch`: reads raw outputs, bins m/z to integers,
  normalises intensities, writes `spectra_all.csv`, `spectra_top20.csv`, `spectra_10pct.csv`
  per molecule per method (including EXP)
- **Section 8** — `diagnose_spectra`: reports which molecules have valid spectra per method
- **Section 9** — `run_comparison`: computes cosine / weighted-dot / Tanimoto / precision / recall
  for each molecule × method × peak-type combination; writes results to `results/`
- **Section 9b** — entropy similarity back-fill (see below)

#### Entropy similarity

`ms_entropy` requires Python ≥ 3.8, but the notebook kernel runs Python 3.6 (NEIMS conda env).
`run_comparison()` handles this automatically: after `compare_spectra.py` finishes it calls
`patch_entropy.py` via the system Python 3.8 at `/appl/opt/python/3.8.14-gnu850/bin/python3.8`,
which has `ms_entropy`, `numpy`, and `pandas` available in the user site-packages.

No manual steps are needed. If you ever need to recompute entropy values directly:

```bash
/appl/opt/python/3.8.14-gnu850/bin/python3.8 src/analysis/patch_entropy.py \
    --base_dir data/simulation_results/franklin_tms
```

Pass `--force` to recompute rows that are already filled.

#### 4. Visualisation (`visualization_for_presentation.ipynb`)

Run all cells. The notebook reads from `data/simulation_results/*/results/` and produces:
- Similarity histograms per metric and peak type
- Heatmaps of mean similarity by method × peak type
- QCxMS variant comparison table (10 ps vs 25 ps vs custom IEE)
- Publication-quality figures saved to `notebook/plots/`

#### 5. Resource benchmark (`5_resource_benchmark.ipynb`)

Reads SLURM accounting data (sacct) or the cached `sacct_cache.json` for the 9-molecule
memory benchmark subset and produces a CPU/memory/wall-time comparison table per method.
Call `save_sacct_cache(job_ids, path)` after new jobs complete to persist the data before
sacct records expire (~1 year).

---

## Experimental spectra

Experimental spectra are stored as `.MSPEC` files (NIST MSP format) in `EXP/{mol}/`.
A symlink `exp.msp → {mol}.MSPEC` must exist for `process_spectra_batch` to find them.
After adding new `.MSPEC` files, run:
```bash
for dir in data/simulation_results/franklin_tms/EXP/*/; do
    mspec=$(ls "$dir"*.MSPEC 2>/dev/null | head -1)
    [ -n "$mspec" ] && [ ! -f "${dir}exp.msp" ] && ln -s "$(basename $mspec)" "${dir}exp.msp"
done
```

---

## Per-molecule directory layout (QCxMS)

```
{dataset}/{variant}/{mol}/
├── GS-opt/
│   ├── xtbopt.sdf          # xtb-optimised geometry
│   └── MS-run/
│       ├── TMPQCXMS/
│       │   └── TMP.{n}/    # One directory per trajectory; contains "ready" when done
│       ├── result.csv       # Raw QCxMS spectrum (m/z, intensity)
│       └── plotms.res       # PlotMS output (quality marker)
└── spectra/
    ├── spectra_all.csv      # Processed spectrum, all peaks
    ├── spectra_top20.csv    # Top 20 peaks
    └── spectra_10pct.csv    # Peaks ≥ 10 % of base peak
```

---

## Legacy / deprecated notes

The `src/depricated/` folder contains earlier versions of the workflow
(pre-notebook, manual bash-driven pipeline). The shell commands in the top
of this README reflect that old workflow and are kept for reference only.
The current workflow is entirely notebook-driven via sections 1–9 above.

---

## TMPQCXMS compression status

QCxMS trajectories are compressed to `TMPQCXMS.tar.gz` per molecule to save disk space.
Compression is managed by `src/workflow/submit_compress_tmpqcxms.sh` with `compress_tmpqcxms_list.txt` as input.

### Current state (as of 2026-07-10)

| Dataset | Compressed (`.tar.gz` only) | Raw dir remaining | Notes |
|---------|----------------------------|-------------------|-------|
| `franklin` | 169 | 28 | Was never in the original list; 14 are dual-state (corrupt .tar.gz + raw dir) |
| `franklin_tms` | 133 | 95 | Job 35378513 hit 2 h limit; 29 dual-state (corrupt .tar.gz + raw dir) |
| `ucb_globes_tracers` | 12 | 14 | `QCxMS_25_ps` (all 13 mols) was never listed; 1 in `QCxMS_10_ps` also missed |

**Dual-state** = both `TMPQCXMS/` directory and `TMPQCXMS.tar.gz` exist.
These archives are **corrupt** (job hit time limit mid-tar).
The updated script now verifies integrity with `tar -tzf` and removes a corrupt archive before re-compressing.

### How to resubmit

`compress_tmpqcxms_list.txt` now covers all 137 remaining raw dirs across all three datasets.
The script time limit has been extended to 12 h.

```bash
cd /scratch/project_2006752/hsandstr/Project/atmospheric-ms-benchmark
N=$(wc -l < src/workflow/compress_tmpqcxms_list.txt)
sbatch --array=1-${N} src/workflow/submit_compress_tmpqcxms.sh
```

After the job completes, clean up `.err`/`.out` files at the project root.

---

## Known archive issues (restrip job 35181425)

During restripping of QCxMS2 auxiliary archives (removing orca.out/geo.out/g98.out to save disk),
several tasks failed. All benchmark data (peaks.csv, qcxms2.log) was unaffected.

| Task | Path | Error | Resolution |
|------|------|-------|------------|
| 3 | franklin/QCxMS2/0003/qcxms2_gfn2 | `mv` to Lustre hit disk quota — archive truncated | Not re-fixed; auxiliary data only, benchmark data intact |
| 6 | franklin/QCxMS2/0006/qcxms2_gfn2 | Same as above | Not re-fixed; benchmark data intact |
| 11 | franklin/QCxMS2/0011/qcxms2_gfn2 | Same as above | Not re-fixed; benchmark data intact |
| 77 | franklin_tms/QCxMS2/0005/qcxms2_gfn2 | NVMe full (archive ~11 GB, 15 GB NVMe) | Resubmitted as job 35182692 with 30 GB NVMe — fixed |
| 78 | franklin_tms/QCxMS2/0006/qcxms2_gfn2 | NVMe full (archive ~15 GB, 15 GB NVMe) | Resubmitted as job 35182692 with 30 GB NVMe — fixed |
| 118 | franklin_tms/QCxMS2/0052/qcxms2_gfn2 | NVMe full (archive ~13 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
| 125 | franklin_tms/QCxMS2/0059/qcxms2_gfn2 | NVMe full (archive ~11 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
| 130 | franklin_tms/QCxMS2_dft/0018/qcxms2_wb97x3c | Corrupt input archive (unexpected EOF) — QCxMS2_dft run crashed mid-archive | Resubmitted as QCxMS2_dft job 35202063 |
| 135 | franklin_tms/memory_benchmark/QCxMS2/0012/qcxms2_gfn2 | NVMe full (archive ~13 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
| 142 | ucb_globes_tracers/QCxMS2/0002/qcxms2_gfn2 | NVMe full (archive ~15 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
| 143 | ucb_globes_tracers/QCxMS2/0004/qcxms2_gfn2 | NVMe full (archive ~13 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
| 145 | ucb_globes_tracers/QCxMS2/0007/qcxms2_gfn2 | NVMe full (archive ~14 GB) | Resubmitted as job 35201234 with 60 GB NVMe — fixed |
