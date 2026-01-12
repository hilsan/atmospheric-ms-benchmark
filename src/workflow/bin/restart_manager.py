#!/usr/bin/env python3
"""
RestartManager for QCxMS2 workflows.

This script checks directories for failed or incomplete runs (CREST or QCxMS2),
cleans necessary directories, and submits the appropriate SLURM scripts.

Usage:
    python restart_manager.py --job-type gfn2 --dirs 0000 0010 --dry-run
"""

import os
import subprocess
import argparse
from pathlib import Path

class RestartManager:
    def __init__(self, base_dir, scripts_dir, job_type, dir_range=None, dry_run=False):
        """
        Initialize the RestartManager.

        Parameters:
            base_dir (str): Base path containing the calculation directories.
            scripts_dir (str): Path to SLURM scripts.
            job_type (str): One of ["crest", "gfn2", "gfn2_gaussian", "wb97x3c_gaussian"].
            dir_range (list of str): Specific directories to check (optional).
            dry_run (bool): If True, do not submit jobs, just print actions.
        """
        self.base_dir = Path(base_dir)
        self.scripts_dir = Path(scripts_dir)
        self.job_type = job_type
        self.dir_range = dir_range
        self.dry_run = dry_run

        # Map job types to SLURM script names
        self.job_scripts = {
            "crest": "batch_crest.sh",
            "gfn2": "qcxms2_gfn2_4.sh",
            "gfn2_gaussian": "qcxms2_gfn2_gaussian.sh",
            "wb97x3c_gaussian": "qcxms2_wB97x3c_gaussian.sh",
        }

        if job_type not in self.job_scripts:
            raise ValueError(f"Unknown job type: {job_type}")

    def check_crest_completed(self, dir_path):
        """Check if CREST finished successfully."""
        crest_dir = dir_path / "QCxMS2"
        if not crest_dir.exists():
            return False
        for out_file in crest_dir.glob("G*out"):
            with open(out_file) as f:
                if "CREST terminated normally" in f.read():
                    return True
        return False

    def check_qcxms2_completed(self, dir_path, sub_dir="qcxms2_gfn2"):
        """Check if QCxMS2 finished successfully in the given subdirectory."""
        qcx_dir = dir_path / "QCxMS2" / sub_dir
        if not qcx_dir.exists():
            return False
        for log_file in qcx_dir.glob("*.log"):
            with open(log_file) as f:
                content = f.read()
                if "QCxMS2 terminated normally" in content:
                    return True
        return False

    def submit_job(self, dir_path, script_name):
        """Submit a SLURM job using sbatch."""
        script_path = self.scripts_dir / script_name
        if not script_path.exists():
            print(f"❌ SLURM script not found: {script_path}")
            return

        if self.dry_run:
            print(f"[Dry-run] Would submit {script_path} in {dir_path}")
            return

        print(f"Submitting {script_path} in {dir_path}")
        subprocess.run(["sbatch", str(script_path)], cwd=dir_path)

    def clean_directory(self, dir_path, sub_dir="qcxms2_gfn2"):
        """Remove files in the QCxMS2 subdirectory before restart."""
        target_dir = dir_path / "QCxMS2" / sub_dir
        if target_dir.exists():
            print(f"Cleaning {target_dir}")
            for item in target_dir.iterdir():
                if item.is_dir():
                    subprocess.run(["rm", "-rf", str(item)])
                else:
                    item.unlink(missing_ok=True)

    def run(self):
        """Run the restart manager over the selected directories."""
        dirs_to_check = self.dir_range or [f"{i:04d}" for i in range(69)]
        for dir_name in dirs_to_check:
            dir_path = self.base_dir / dir_name
            if not dir_path.exists():
                print(f"Directory {dir_path} does not exist. Skipping.")
                continue

            # CREST check
            if self.job_type == "crest":
                if not self.check_crest_completed(dir_path):
                    print(f"CREST failed or incomplete in {dir_path}")
                    self.submit_job(dir_path / "QCxMS2", self.job_scripts["crest"])
                else:
                    print(f"CREST completed normally in {dir_path}")
                continue

            # QCxMS2 checks
            sub_dir_map = {
                "gfn2": "qcxms2_gfn2",
                "gfn2_gaussian": "qcxms2_gfn2",
                "wb97x3c_gaussian": "qcxms2_wb97x3c",
            }
            sub_dir = sub_dir_map.get(self.job_type, "")
            if self.check_qcxms2_completed(dir_path, sub_dir=sub_dir):
                print(f"QCxMS2 {self.job_type} completed normally in {dir_path}. Skipping.")
                continue

            # Otherwise, clean and submit
            print(f"QCxMS2 {self.job_type} failed or incomplete in {dir_path}")
            self.clean_directory(dir_path, sub_dir=sub_dir)
            self.submit_job(dir_path / "QCxMS2", self.job_scripts[self.job_type])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manage restarts of CREST/QCxMS2 runs.")
    parser.add_argument("--base-dir", required=True, help="Base directory with calculation subfolders")
    parser.add_argument("--scripts-dir", required=True, help="Directory containing SLURM scripts")
    parser.add_argument("--job-type", required=True,
                        choices=["crest", "gfn2", "gfn2_gaussian", "wb97x3c_gaussian"],
                        help="Type of job to manage")
    parser.add_argument("--dirs", nargs="*", help="Specific directories to check (e.g., 0000 0001)")
    parser.add_argument("--dry-run", action="store_true", help="Do not actually submit jobs, just print actions")
    args = parser.parse_args()

    manager = RestartManager(
        base_dir=args.base_dir,
        scripts_dir=args.scripts_dir,
        job_type=args.job_type,
        dir_range=args.dirs,
        dry_run=args.dry_run
    )
    manager.run()

