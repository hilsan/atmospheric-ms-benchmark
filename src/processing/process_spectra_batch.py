"""
Batch spectrum processing: bin, scale, and filter spectra for all
molecules and datasets in a simulation results directory.
"""

import os
import subprocess

EXPECTED_SPECTRA = ["spectra_all.csv", "spectra_10pct.csv", "spectra_top20.csv"]


def get_input_file(dataset, subdir, subdir_path):
    """Return the expected raw spectrum file path for a given dataset/molecule."""
    dl = dataset.lower()
    if dl == "neims":
        return os.path.join(subdir_path, "annotated.sdf")
    elif dl == "exp":
        return os.path.join(subdir_path, "exp.msp")
    elif dl in ("cfm-id", "cfmid"):
        return os.path.join(subdir_path, f"{subdir}.log")
    elif dl == "qcxms2":
        gfn2_dir = os.path.join(subdir_path, "qcxms2_gfn2")
        peaks = os.path.join(gfn2_dir, "peaks.csv")
        aux   = os.path.join(gfn2_dir, "qcxms2_auxiliary.tar.gz")
        return peaks if os.path.exists(aux) else None
    else:  # QCxMS variants
        return os.path.join(subdir_path, "GS-opt", "MS-run", "result.csv")


def spectra_complete(spectra_dir):
    """Return True if all expected spectrum output files exist."""
    return all(os.path.exists(os.path.join(spectra_dir, f)) for f in EXPECTED_SPECTRA)


def process_spectra_batch(sim_base_dir, datasets, override=False):
    """Bin and scale spectra for all molecules in each dataset.

    Calls bin_scale_intensities.py for each molecule, producing
    spectra_all.csv, spectra_top20.csv, and spectra_10pct.csv.

    Parameters
    ----------
    sim_base_dir : str
        Root directory containing one subfolder per method/dataset.
    datasets : list[str]
        Dataset names to process (e.g. ["CFMID", "EXP", "NEIMS", "QCxMS_25_ps"]).
    override : bool
        If True, reprocess even if output spectra already exist.
    """
    sim_base_dir = os.path.abspath(sim_base_dir)
    script_path = os.path.abspath(
        os.path.join(sim_base_dir, "../../../src/processing/bin_scale_intensities.py")
    )

    print(f"Simulation base directory: {sim_base_dir}")
    print(f"Processing script:         {script_path}")
    print(f"Override existing:         {override}\n")

    for dataset in datasets:
        dataset_path = os.path.join(sim_base_dir, dataset)
        if not os.path.exists(dataset_path):
            print(f"[{dataset}] Folder missing, skipping.")
            continue

        skipped = processed = failed_input = 0

        for i in range(61):
            subdir = f"{i:04d}"
            subdir_path = os.path.join(dataset_path, subdir)
            if not os.path.exists(subdir_path):
                continue

            spectra_dir = os.path.join(subdir_path, "spectra")

            if not override and spectra_complete(spectra_dir):
                skipped += 1
                continue

            input_file = get_input_file(dataset, subdir, subdir_path)
            if input_file is None or not os.path.exists(input_file):
                print(f"  [{dataset}/{subdir}] Input not found: {input_file}")
                failed_input += 1
                continue

            os.makedirs(spectra_dir, exist_ok=True)
            result = subprocess.run(
                [
                    "python", script_path,
                    "-i", input_file,
                    "-o", os.path.join(spectra_dir, "spectra"),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if result.returncode != 0:
                print("  [{}/{}] ERROR:\n{}".format(dataset, subdir, result.stderr.decode()))
            else:
                print("  [{}/{}] Done".format(dataset, subdir))
            processed += 1

        print(
            f"[{dataset}] processed={processed}  "
            f"skipped={skipped}  missing_input={failed_input}"
        )
