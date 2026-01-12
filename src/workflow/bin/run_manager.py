#!/usr/bin/env python3
"""
run_manager.py

A Python replacement for run_all_qcxms.sh:
- Iterates over molecule directories
- Submits QCxMS fragmentation jobs
- Extracts SLURM job ID from output
- Polls SLURM until jobs complete
- Prints clean status messages

Run:
    python workflow/run_manager.py --start 16 --end 68 \
        --base-dir /path/to/base \
        --submit-script $SCRIPTS_NEIMS/submit_frag_parallell.sh

You should have already done:
    export PYTHONPATH=$SRC_PATH:$PYTHONPATH
"""

import os
import re
import time
import argparse
import subprocess
from pathlib import Path


# ---------------------------
# Helper utilities
# ---------------------------

def pad_dirnum(n: int) -> str:
    """Convert an integer like 16 → '0016' to match directory names."""
    return f"{n:04d}"


def run_subprocess(cmd, cwd=None):
    """Run command and return stdout as string. Raise on failure."""
    result = subprocess.run(
        cmd,
        cwd=cwd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    return result.stdout


def extract_job_id(submission_output: str):
    """
    Extract SLURM job ID from submission output.
    Matches typical:
        'Submitted batch job 123456'
    """
    match = re.search(r"Submitted batch job\s+(\d+)", submission_output)
    return match.group(1) if match else None


def slurm_job_running(job_id: str) -> bool:
    """Return True if job is still present in squeue."""
    result = subprocess.run(
        ["squeue", "-j", job_id],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return job_id in result.stdout


# ---------------------------
# Main workflow
# ---------------------------

def process_directory(dirnum: str, base_dir: Path, submit_script: Path,
                      poll_interval: int, dry_run: bool):
    """
    Handle one directory:
    - validate path
    - submit job
    - poll for completion
    """
    run_dir = base_dir / dirnum / "QCxMS" / "MS-run"
    print(f"\n📁 Processing: {run_dir}")

    if not run_dir.is_dir():
        print(f"❌ Directory missing, skipping.")
        return

    if dry_run:
        print(f"DRY RUN: Would run {submit_script} in {run_dir}")
        return

    # Submit job
    try:
        print(f"🚀 Submitting job via {submit_script}")
        output = run_subprocess([str(submit_script)], cwd=run_dir)
        print(output)
    except Exception as e:
        print(f"❌ Failed to submit job: {e}")
        return

    # Extract job ID
    job_id = extract_job_id(output)
    if not job_id:
        print("❌ Could not parse SLURM job ID. Skipping.")
        return

    print(f"🆔 SLURM Job ID: {job_id}")
    print(f"🕒 Polling every {poll_interval} seconds...")

    # Poll
    while slurm_job_running(job_id):
        print(f"⏳ Job {job_id} still running...")
        time.sleep(poll_interval)

    print(f"✅ Job {job_id} finished.")


# ---------------------------
# CLI
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit and monitor QCxMS fragmentation jobs over multiple directories."
    )
    parser.add_argument("--start", type=str, required=True,
                        help="Start directory number (e.g., 16 or 0016)")
    parser.add_argument("--end", type=str, required=True,
                        help="End directory number (e.g., 68 or 0068)")
    parser.add_argument("--base-dir", type=str, required=True,
                        help="Base directory containing numbered dirs.")
    parser.add_argument("--submit-script", type=str, required=True,
                        help="SLURM submission script (full path).")
    parser.add_argument("--poll", type=int, default=600,
                        help="Polling interval in seconds (default 600).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Do everything except submit jobs.")

    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    submit_script = Path(args.submit_script).resolve()

    if not base_dir.is_dir():
        raise NotADirectoryError(f"Invalid base directory: {base_dir}")

    if not submit_script.is_file():
        raise FileNotFoundError(f"Submit script not found: {submit_script}")

    # Convert start/end to ints
    start_n = int(args.start)
    end_n = int(args.end)

    for n in range(start_n, end_n + 1):
        dirnum = pad_dirnum(n)
        process_directory(
            dirnum=dirnum,
            base_dir=base_dir,
            submit_script=submit_script,
            poll_interval=args.poll,
            dry_run=args.dry_run
        )

    print("\n🎉 All directories processed.\n")


if __name__ == "__main__":
    main()

