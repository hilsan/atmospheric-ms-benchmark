# Molecular Simulation Workflow

This repository contains scripts to run **QCxMS2**, **QCxMS**, and **NEIMS** simulations, manage restarts, and post-process spectra. The scripts are organized into **SLURM submission scripts**, **Python workflow managers**, and **post-processing utilities**.

---

## Table of Contents

1. [QCxMS2](#qcxms2)  
2. [QCxMS](#qcxms)  
3. [NEIMS](#neims)  
4. [Post-Processing](#post-processing)  
5. [Environment Setup](#environment-setup)  

---

## QCxMS2

**Purpose:** Full QCxMS2 runs including CREST sampling, geometry optimization, and QCxMS2 spectra calculations.

### Submission Scripts

- **`submit_batch_crest.sh`**  
  - Submits CREST sampling jobs for multiple molecules.  
  - Handles restart automatically.  
  - **Usage:**  
    ```bash
    sbatch submit_batch_crest.sh
    ```
  - **Notes:** Adjust `-T` (threads) and `--chrg` (charge) as needed.

- **`submit_qcxms2_job.sh`**  
  - Runs QCxMS2 on a given molecule.  
  - Accepts flags for levels and methods: `-i INPUT_FILE -g GEOLEVEL -t TSLEVEL -p IPLEVEL -d DIST_METHOD`.  
  - **Usage example:**  
    ```bash
    ./submit_qcxms2_job.sh -i in.xyz -g 6 -t 6 -p 6 -d gfn2
    ```

- **`submit_qcxms_gs_md.sh`**  
  - Submits GS-level molecular dynamics sampling jobs.  
  - Handles array jobs and creates necessary QCxMS input files.  
  - **Usage:**  
    ```bash
    sbatch submit_qcxms_gs_md.sh
    ```

### Python Managers

- **`restart_manager.py`**  
  - Detects failed or incomplete QCxMS2 / CREST runs and resubmits them.  
  - **Usage example:**  
    ```bash
    python restart_manager.py --base-dir /path/to/molecules \
                              --scripts-dir /path/to/slurm_scripts \
                              --job-type gfn2 \
                              --dirs 0000 0010 \
                              --dry-run
    ```
  - **Job Types:** `crest`, `gfn2`, `gfn2_gaussian`, `wb97x3c_gaussian`.

- **`run_manager.py`**  
  - Submits and monitors QCxMS fragmentation jobs across multiple directories.  
  - Works with **`submit_qcxms_fragmentation.sh`** (previously `submit_frag_parallell.sh`).  
  - **Usage example:**  
    ```bash
    python run_manager.py --start 16 --end 68 \
                          --base-dir /path/to/base \
                          --submit-script /path/to/submit_qcxms_fragmentation.sh \
                          --poll 600
    ```
  - Use `--dry-run` to simulate without submission.

---

## QCxMS

**Purpose:** Standalone QCxMS fragmentation jobs, usually after MD/CREST sampling.

- **`submit_qcxms_fragmentation.sh`**  
  - Handles array submission for multiple molecule directories.  
  - Generates a temporary SLURM script and submits jobs in parallel.  
  - **Usage example:**  
    ```bash
    ./submit_qcxms_fragmentation.sh --cpus 4 --mem 16G --time 12:00:00 --max-array 200
    ```
  - **Notes:** Adjust `--cpus`, `--mem`, `--time`, and `--max-array` according to cluster limits.  
  - **Logs:** Each array task creates `frag_array_task_<ID>.log`. Essential output is kept; auxiliary files are archived.

---

## NEIMS

**Purpose:** Predict molecular spectra using neural network models.

- **`submit_neims_array.sh`**  
  - Submits NEIMS jobs as a SLURM array.  
  - Reads `smiles.sdf` in each directory and generates `annotated.sdf`.  
  - **Usage:**  
    ```bash
    sbatch submit_neims_array.sh
    ```
  - **Notes:** Array indices are controlled via `#SBATCH --array=0-68`. Adjust to match number of molecules.

---

## Post-Processing

**Purpose:** Automate scaling, comparison, and plotting of simulated vs. experimental spectra.

- **`postprocessing_manager.py`**  
  - Performs the following steps:
    1. Intensity scaling (`scale_intensities.py`)  
    2. Spectra comparison (`compare_spectra.py` or `_nomax.py`)  
    3. Similarity histogram plotting (`plot_similarity_histograms.py`)  
  - **Usage example:**  
    ```bash
    python postprocessing_manager.py --dir /path/to/base --nomax
    ```
  - Optional flags:
    - `--skip-scaling`
    - `--skip-comparison`
    - `--skip-plots`

---

## Environment Setup

- Python 3.10 required for NEIMS and post-processing.  
- Modules required for QCxMS2 / CREST:
  - `gcc`, `openmpi`, `intel-oneapi-mkl`, `biopythontools` (check versions in scripts)  
- Recommended to use a conda environment for NEIMS:
  ```bash
  conda activate neims

