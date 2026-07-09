"""
Helper functions for submitting and monitoring HPC/SLURM jobs
in the atmospheric-ms-benchmark workflow.
"""

import os
import re
import subprocess


def _queued_dirs():
    """Return set of working directories for all jobs currently in the queue."""
    result = subprocess.run(
        ["squeue", "-u", os.environ.get("USER", ""),
         "--format=%Z", "--noheader"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return set(result.stdout.strip().splitlines())


def submit_slurm_array(sim_dir, script_path, array_spec,
                        done_marker=None, n_mols=None):
    """Submit a SLURM array job, skipping already-complete molecules.

    Parameters
    ----------
    sim_dir : str
        Working directory for sbatch.
    script_path : str
        Absolute path to the SLURM script.
    array_spec : str
        Default value for --array, e.g. "0-60". Overridden when done_marker
        is provided and some molecules are already complete.
    done_marker : str, optional
        Path relative to each molecule folder that signals completion,
        e.g. "GS-opt/xtbopt.sdf" (GS-MD) or "annotated.sdf" (NEIMS).
        When given, only pending molecules are included in the array.
    n_mols : int, optional
        Total number of molecules to check (required when done_marker is set).
    """
    os.makedirs(sim_dir, exist_ok=True)
    sim_dir = os.path.abspath(sim_dir)

    if done_marker is not None and n_mols is not None:
        pending, done_count = [], 0
        for i in range(n_mols):
            mol = f"{i:04d}"
            marker_path = os.path.join(sim_dir, mol, done_marker)
            tar_path = marker_path + ".tar.gz"
            if os.path.exists(marker_path) or os.path.exists(tar_path):
                if not os.path.exists(marker_path) and os.path.exists(tar_path):
                    print(f"  [ARCHIVED] {mol}: {done_marker} compressed as tar.gz")
                done_count += 1
            else:
                pending.append(str(i))
        print(f"  Done: {done_count}/{n_mols}  |  Pending: {len(pending)}")
        if not pending:
            print(f"[SKIP] All {n_mols} molecules already done ({done_marker})")
            return
        array_spec = ",".join(pending)
        print(f"  Submitting for: {', '.join(f'{int(p):04d}' for p in pending)}")

    if sim_dir in _queued_dirs():
        print(f"[SKIP] {sim_dir} already in queue")
        return
    print(f"Using simulation directory: {sim_dir}")
    result = subprocess.run(
        ["sbatch", f"--array={array_spec}", script_path],
        cwd=sim_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    print(result.stdout)
    if result.stderr:
        print(result.stderr)


def _plotms_passes(msrun_dir, min_frag_frac=0.90, min_base_peak=95):
    """Return True if plotms.res exists and passes quality gates.

    When TMPQCXMS is archived (tar.gz), the trajectory threshold is skipped
    and only plotms.res quality is checked.
    """
    plotms_res = os.path.join(msrun_dir, "plotms.res")
    if not os.path.exists(plotms_res):
        return False
    n_ready, n_total = _frag_completion(msrun_dir)
    if (n_ready, n_total) == _ARCHIVED:
        n_runs, base = _parse_plotms_quality(plotms_res)
        return (n_runs is not None and base is not None and base >= min_base_peak)
    threshold = int(n_total * min_frag_frac) if n_total else 0
    n_runs, base = _parse_plotms_quality(plotms_res)
    return (n_runs is not None and n_runs >= threshold
            and base is not None and base >= min_base_peak)


def submit_qcxms_frag_jobs(sim_dir, script_path, n_mols=61, skip_done_check=False):
    """Submit QCxMS fragmentation jobs for molecules that still need work.

    Skips molecules where plotms.res already passes quality (filesystem check,
    no squeue needed). For the remainder, does a single squeue call and skips
    anything already queued. Only then submits the rest.

    Parameters
    ----------
    skip_done_check : bool
        If True, bypass the plotms.res quality check entirely and call the
        bash script for every mol that has an MS-run directory (subject only
        to the queue check). The bash script itself skips individual TMP.xxx
        dirs that already have normal termination, so only unfinished
        trajectories are re-submitted. plotms.res files are left untouched.
    """
    done, no_dir, to_check = [], [], []
    for i in range(n_mols):
        mol = f"{i:04d}"
        mol_path = os.path.abspath(os.path.join(sim_dir, mol, "GS-opt", "MS-run"))
        if not os.path.isdir(mol_path):
            no_dir.append(mol)
            continue
        if not skip_done_check and _plotms_passes(mol_path):
            done.append(mol)
            continue
        to_check.append((mol, mol_path))

    label = "plotms OK" if not skip_done_check else "skipped (FORCE_FRAG)"
    print(f"  Done ({label}): {len(done)}  |  No dir: {len(no_dir)}  |  To check: {len(to_check)}")

    if not to_check:
        print("  Nothing to submit.")
        return

    queued = _queued_dirs()
    already_queued, to_submit = [], []
    for mol, mol_path in to_check:
        tmpqcxms = os.path.join(mol_path, "TMPQCXMS")
        tmptar = os.path.join(mol_path, "TMPQCXMS.tar.gz")
        if os.path.isfile(tmptar):
            done.append(mol)
            print(f"  [ARCHIVED] {mol}: TMPQCXMS.tar.gz found — treating as done")
            continue
        if tmpqcxms in queued or mol_path in queued:
            already_queued.append(mol)
        else:
            to_submit.append((mol, mol_path))

    print(f"  Already queued:   {len(already_queued)}  |  Submitting: {len(to_submit)}")
    if already_queued:
        print(f"  Queued:  {' '.join(already_queued)}")

    for mol, mol_path in to_submit:
        result = subprocess.run(
            ["bash", script_path],
            cwd=mol_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        print(f"  Submitted {mol}: {result.stdout.strip()}")
        if result.stderr:
            print(f"    stderr: {result.stderr.strip()}")


_ARCHIVED = (-1, -1)  # sentinel: TMPQCXMS compressed, treat as fully complete


def _frag_completion(msrun_dir):
    """Return (n_ready, n_total) trajectory counts for a MS-run directory.

    Returns _ARCHIVED (-1, -1) when TMPQCXMS has been compressed to
    TMPQCXMS.tar.gz — callers should treat this as fully complete.
    """
    tmpdir = os.path.join(msrun_dir, "TMPQCXMS")
    tarfile = os.path.join(msrun_dir, "TMPQCXMS.tar.gz")
    if os.path.isfile(tarfile):
        return _ARCHIVED
    if not os.path.isdir(tmpdir):
        return 0, 0
    subdirs = [d for d in os.listdir(tmpdir)
               if os.path.isdir(os.path.join(tmpdir, d)) and d.startswith("TMP.")]
    n_total = len(subdirs)
    n_ready = sum(
        1 for d in subdirs
        if os.path.exists(os.path.join(tmpdir, d, "ready"))
    )
    return n_ready, n_total


def _parse_plotms_quality(plotms_res_path):
    """Parse plotms.res and return (n_successful, base_peak_counts).

    Looks for:
      ``X  successfull runs!!``
      ``Theoretical counts in 100 % signal:   X``

    Returns (None, None) if either value cannot be parsed.
    """
    try:
        text = open(plotms_res_path).read()
    except OSError:
        return None, None
    m_runs = re.search(r"(\d+)\s+successfull runs!!", text)
    m_base = re.search(r"Theoretical counts in 100 %\s+signal:\s+(\d+)", text)
    n_runs = int(m_runs.group(1)) if m_runs else None
    base   = int(m_base.group(1)) if m_base else None
    return n_runs, base


def run_plotms_for_included(wrkdir, plotms_bin=None,
                             min_frag_frac=0.90, min_base_peak=95):
    """Parse analysis_decision.txt and run PlotMS for all INCLUDE folders.

    Skips a molecule if ``plotms.res`` already exists **and** passes quality:
      - successful runs ≥ ``min_frag_frac`` × total trajectories
      - base peak (theoretical counts at 100 % signal) ≥ ``min_base_peak``

    Only attempts PlotMS when ≥ ``min_frag_frac`` of TMPQCXMS trajectories
    have a ``ready`` flag.

    Parameters
    ----------
    wrkdir : str
        Root directory containing molecule folders and analysis_decision.txt.
    plotms_bin : str, optional
        Path to the PlotMS binary. Defaults to ~/PlotMS.v.6.2.0/plotms.
    min_frag_frac : float
        Minimum fraction of completed trajectories to attempt PlotMS (0.90).
    min_base_peak : int
        Minimum base-peak counts required to mark a run as good (95).
    """
    if plotms_bin is None:
        plotms_bin = os.path.expanduser("~/PlotMS.v.6.2.0/plotms")

    include_file = os.path.join(wrkdir, "analysis_decision.txt")
    if not os.path.exists(include_file):
        print(f"  analysis_decision.txt not found in {wrkdir} — skipping PlotMS (run fragmentation first)")
        return
    pattern = re.compile(r"^(\d{4})\s*:.*->\s*INCLUDE")
    include_folders = set()
    with open(include_file) as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                include_folders.add(match.group(1))

    print(f"Folders to process ({len(include_folders)}): {sorted(include_folders)}")

    for folder in sorted(os.listdir(wrkdir)):
        if folder not in include_folders:
            continue
        msrun_dir = os.path.join(wrkdir, folder, "GS-opt", "MS-run")
        if not os.path.isdir(msrun_dir):
            print(f"[SKIP] {folder}: MS-run not found")
            continue

        output_file = os.path.join(msrun_dir, "plotms.res")
        n_ready, n_total = _frag_completion(msrun_dir)
        archived = (n_ready, n_total) == _ARCHIVED

        # Check if a prior PlotMS run already produced good-quality output
        if os.path.exists(output_file):
            n_runs, base = _parse_plotms_quality(output_file)
            threshold = int(n_total * min_frag_frac) if not archived and n_total else 0
            if (n_runs is not None and (archived or n_runs >= threshold)
                    and base is not None and base >= min_base_peak):
                suffix = " [TMPQCXMS archived]" if archived else ""
                print(f"[SKIP] {folder}: already done "
                      f"({n_runs} runs, base peak {base}){suffix}")
                continue
            else:
                print(f"[REDO] {folder}: prior run below quality "
                      f"({n_runs} runs, base peak {base}) — re-running")

        # Pre-flight: enough trajectories completed?
        if archived:
            print(f"[SKIP] {folder}: TMPQCXMS archived but no valid plotms.res — skipping PlotMS")
            continue
        if n_total == 0:
            print(f"[SKIP] {folder}: no TMPQCXMS trajectories found")
            continue
        frac = n_ready / n_total
        if frac < min_frag_frac:
            print(f"[SKIP] {folder}: {n_ready}/{n_total} ({frac:.0%}) trajectories "
                  f"complete (need ≥{min_frag_frac:.0%})")
            continue

        print(f"Running PlotMS for {folder} ({n_ready}/{n_total} trajectories)...")
        try:
            with open(output_file, "w") as out_f:
                subprocess.run(
                    [plotms_bin],
                    cwd=msrun_dir,
                    stdout=out_f,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            n_runs, base = _parse_plotms_quality(output_file)
            if n_runs is not None and base is not None:
                ok = (n_runs >= int(n_total * min_frag_frac) and base >= min_base_peak)
                status = "OK" if ok else "LOW QUALITY"
                print(f"  Done [{status}] — {n_runs} runs, base peak {base}")
            else:
                print(f"  Done → {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"  [ERROR] PlotMS failed for {folder}: {e}")


def submit_crest_jobs(sim_dir, script_path, folders):
    """Submit CREST conformer search jobs for a list of molecule folders.

    Skips molecules where CREST already terminated normally.

    Parameters
    ----------
    sim_dir : str
        Base directory containing molecule folders.
    script_path : str
        Absolute path to the CREST SLURM script.
    folders : list[str]
        Folder names to process (e.g. ["0005", "0013"]).
    """
    queued = _queued_dirs()
    for mol in folders:
        mol_path = os.path.join(sim_dir, mol)
        if not os.path.isdir(mol_path):
            print(f"Skipping {mol}: directory not found")
            continue
        if not os.path.isfile(os.path.join(mol_path, "structure.sdf")):
            print(f"Skipping {mol}: structure.sdf not found")
            continue

        if mol_path in queued:
            print(f"{mol}: already in queue, skipping.")
            continue

        finished = False
        for log_name in ("crest.log", "crest_restart.log"):
            log_path = os.path.join(mol_path, log_name)
            if os.path.exists(log_path):
                with open(log_path) as f:
                    if "CREST terminated normally." in f.read():
                        finished = True
                        break
        if finished:
            print(f"{mol}: CREST already finished, skipping.")
            continue

        print(f"{mol}: Submitting CREST job...")
        result = subprocess.run(
            ["sbatch", script_path],
            cwd=mol_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode == 0:
            print(f"  {result.stdout.strip()}")
        else:
            print(f"  Submission failed: {result.stderr.strip()}")


def submit_qcxms2_jobs(base_dir, bash_script, method="gfn2", num_folders=61, max_concurrent=None):
    """Submit QCxMS2 jobs for molecules where CREST has finished.

    Skips molecules where QCxMS2 already completed, CREST is not done, or
    the job is already in the SLURM queue.

    Parameters
    ----------
    max_concurrent : int or None
        If set, use SLURM afterany dependencies to keep at most this many
        jobs running simultaneously (rolling window over submitted job IDs).
    """
    base_dir = os.path.abspath(base_dir)
    queued = _queued_dirs()

    # Seed the rolling window with job IDs already in the queue for this base_dir,
    # so new submissions correctly chain onto existing ones.
    submitted_ids = []
    if max_concurrent:
        result = subprocess.run(
            ["squeue", "-u", os.environ.get("USER", ""),
             "--format=%i %Z", "--noheader"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        for line in result.stdout.strip().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].startswith(base_dir):
                submitted_ids.append(parts[0])
        if submitted_ids:
            print(f"  Found {len(submitted_ids)} already-queued jobs for this dataset")

    for i in range(num_folders):
        folder = f"{i:04d}"
        workdir = os.path.join(base_dir, folder)
        if not os.path.isdir(workdir):
            print(f"Skipping {folder}: directory not found")
            continue

        if workdir in queued:
            print(f"Skipping {folder}: already in queue")
            continue

        out_dir = os.path.join(workdir, f"qcxms2_{method}")

        perm_fail = os.path.join(out_dir, "QCXMS2_PERMANENT_FAIL")
        if os.path.exists(perm_fail):
            print(f"Skipping {folder}: permanently failed (msreact/unrecoverable error)")
            continue

        # tar archive is only written after a successful run completes — reliable done marker
        tar_done = os.path.join(out_dir, "qcxms2_auxiliary.tar.gz")
        if os.path.exists(tar_done):
            print(f"Skipping {folder}: QCxMS2 already completed (tar exists)")
            continue
        # fallback: log-based check for runs completed before NVMe workflow
        qcxms2_log = os.path.join(out_dir, "qcxms2.log")
        if os.path.exists(qcxms2_log):
            with open(qcxms2_log) as f:
                if "QCxMS2 terminated normally" in f.read():
                    print(f"Skipping {folder}: QCxMS2 already completed (log)")
                    continue

        crest_file = os.path.join(workdir, "crest_best.xyz")
        if not os.path.exists(crest_file):
            print(f"Skipping {folder}: CREST not finished (no crest_best.xyz)")
            continue
        crest_done = False
        for log_name in ("crest.log", "crest_restart.log"):
            log_path = os.path.join(workdir, log_name)
            if os.path.exists(log_path):
                with open(log_path) as f:
                    if "CREST terminated normally" in f.read():
                        crest_done = True
                        break
        if not crest_done:
            print(f"Skipping {folder}: CREST not finished (log incomplete)")
            continue

        cmd = ["sbatch"]
        if max_concurrent and len(submitted_ids) >= max_concurrent:
            # depend on the job that was submitted max_concurrent slots ago
            dep_id = submitted_ids[-max_concurrent]
            cmd += [f"--dependency=afterany:{dep_id}"]
        cmd += [bash_script, base_dir, method]

        print(f"Submitting {folder} ({method})")
        result = subprocess.run(
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode == 0:
            print(f"  {result.stdout.strip()}")
            # extract job ID from "Submitted batch job 12345"
            parts = result.stdout.strip().split()
            if parts and parts[-1].isdigit():
                submitted_ids.append(parts[-1])
        else:
            print(f"  Submission failed: {result.stderr.strip()}")


def submit_cfmid_job(cfmid_dir, script_path, n_mols=None):
    """Submit a CFMID batch job via sbatch. Skips if all molecules already have output."""
    cfmid_dir = os.path.abspath(cfmid_dir)

    if n_mols is not None:
        done = [
            f"{i:04d}" for i in range(n_mols)
            if os.path.exists(os.path.join(cfmid_dir, f"{i:04d}", f"{i:04d}.log"))
        ]
        print(f"  CFMID done: {len(done)}/{n_mols}")
        if len(done) == n_mols:
            print(f"[SKIP] CFMID already completed for all {n_mols} molecules")
            return
        if done:
            print(f"  Completed: {' '.join(done)}")

    result = subprocess.run(
        [
            "sbatch",
            f"-o{cfmid_dir}/cfmid_%j.out",
            f"-e{cfmid_dir}/cfmid_%j.err",
            script_path,
            cfmid_dir,
        ],
        cwd=cfmid_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if result.returncode != 0:
        print("Submission failed:")
        print(result.stderr)
    else:
        print("Submitted:", result.stdout.strip())


def write_cfmid_idx_smiles(base_dir):
    """Write idx_smiles.txt from per-molecule smiles.smi files.

    Each output line: <folder_name> <smiles>
    """
    from pathlib import Path

    base_dir = Path(base_dir)
    out_file = base_dir / "idx_smiles.txt"
    lines = []
    for d in sorted(p for p in base_dir.iterdir() if p.is_dir()):
        smi_file = d / "smiles.smi"
        if not smi_file.exists():
            continue
        smi = smi_file.read_text().strip()
        if smi:
            lines.append(f"{d.name} {smi}")
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Wrote {len(lines)} molecules to {out_file}")


def check_crest_status(datasets):
    """Check CREST completion status across multiple datasets.

    Parameters
    ----------
    datasets : dict[str, str]
        Mapping of dataset name → absolute path to the QCxMS2 directory.
    """
    squeue_result = subprocess.run(
        ["squeue", "-u", os.environ.get("USER", ""), "--format=%Z", "--noheader"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    running_dirs = set(squeue_result.stdout.strip().splitlines())

    for dataset_name, dataset in datasets.items():
        print(f"\n{'='*60}")
        print(f"Dataset: {dataset_name}")
        print(f"{'='*60}")
        finished, still_running, not_started, missing = [], [], [], []
        for i in range(61):
            mol = f"{i:04d}"
            mol_path = os.path.join(dataset, mol)
            if not os.path.isdir(mol_path):
                missing.append(mol)
                continue
            done = False
            for log_name in ("crest.log", "crest_restart.log"):
                log_path = os.path.join(mol_path, log_name)
                if os.path.exists(log_path):
                    with open(log_path) as f:
                        if "CREST terminated normally." in f.read():
                            done = True
                            break
            if done:
                finished.append(mol)
            elif mol_path in running_dirs:
                still_running.append(mol)
            else:
                not_started.append(mol)

        print(f"  Finished:       {len(finished):3d} — {finished}")
        print(f"  Still running:  {len(still_running):3d} — {still_running}")
        print(f"  Not done/stuck: {len(not_started):3d} — {not_started}")
        print(f"  Missing dir:    {len(missing):3d} — {missing}")


_QCXMS_GSMD_SCRIPTS = {
    "QCxMS_10_ps":       "submit_qcxms_gs_md_10_ps.sh",
    "QCxMS_25_ps":       "submit_qcxms_gs_md_25_ps.sh",
    "QCxMS_10_ps_iee03": "submit_qcxms_gs_md_10_ps_iee03.sh",
}


def setup_memory_benchmark(bench_dir, sim_base, bench_mols, qcxms_variants, src_root):
    """Set up the memory benchmark directory and submit fast jobs.

    Creates per-molecule subdirs, symlinks structure.sdf from the production
    run, then submits NEIMS, CFMID, and QCxMS GS-MD array jobs.
    CREST/QCxMS2 are slow — submit separately via submit_crest_for_benchmark.

    Returns a list of SLURM job ID strings.
    """
    from pathlib import Path
    bench_dir = Path(bench_dir)
    sim_base  = Path(sim_base)
    src_root  = Path(src_root)
    task_ids  = ",".join(str(int(rep["mol_id"])) for rep in bench_mols)
    job_ids   = []

    def _symlink_sdf(src_sub, dst_sub, mol_id):
        dst_dir = bench_dir / dst_sub / mol_id
        dst_dir.mkdir(parents=True, exist_ok=True)
        src = sim_base / src_sub / mol_id / "structure.sdf"
        dst = dst_dir / "structure.sdf"
        if not dst.exists() and src.exists():
            dst.symlink_to(src)

    def _submit_array(subdir, script_name):
        result = subprocess.run(
            ["sbatch", f"--array={task_ids}",
             str(src_root / "workflow" / script_name)],
            cwd=str(bench_dir / subdir),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        print(result.stdout.strip() or result.stderr.strip())
        m = re.search(r"\b(\d{5,})\b", result.stdout)
        if m:
            job_ids.append(m.group(1))

    # NEIMS
    for rep in bench_mols:
        _symlink_sdf("NEIMS", "NEIMS", rep["mol_id"])
    print(f"Submitting NEIMS array ({task_ids}) ...")
    _submit_array("NEIMS", "submit_neims_array.sh")

    # CFMID
    cfmid_bench = bench_dir / "CFMID"
    cfmid_bench.mkdir(parents=True, exist_ok=True)
    idx_lines = [f"{rep['mol_id']} {rep['smiles']}" for rep in bench_mols]
    (cfmid_bench / "idx_smiles.txt").write_text("\n".join(idx_lines) + "\n")
    print(f"\nSubmitting CFMID ({len(idx_lines)} molecules) ...")
    result = subprocess.run(
        ["sbatch", str(src_root / "workflow" / "run_cfmid.sh"), str(cfmid_bench)],
        cwd=str(cfmid_bench),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    print(result.stdout.strip() or result.stderr.strip())
    m = re.search(r"\b(\d{5,})\b", result.stdout)
    if m:
        job_ids.append(m.group(1))

    # QCxMS GS-MD (all variants)
    for variant in qcxms_variants:
        for rep in bench_mols:
            _symlink_sdf(variant, variant, rep["mol_id"])
        script = _QCXMS_GSMD_SCRIPTS.get(variant)
        if not script:
            print(f"Warning: no GS-MD script configured for {variant}, skipping.")
            continue
        print(f"\nSubmitting QCxMS GS-MD ({variant}): {task_ids}")
        _submit_array(variant, script)

    # QCxMS2: symlink structure.sdf only (CREST submitted separately)
    for rep in bench_mols:
        _symlink_sdf("QCxMS2", "QCxMS2", rep["mol_id"])

    print(f"\nBenchmark directory: {bench_dir}")
    print(f"Submitted job IDs: {job_ids}")
    print("After GS-MD finishes, submit fragmentation manually for each molecule.")
    return job_ids


def submit_crest_for_benchmark(bench_dir, bench_mols, crest_script):
    """Submit CREST for benchmark molecules that don't have crest_best.xyz yet.

    Returns a list of submitted SLURM job ID strings.
    """
    from pathlib import Path
    bench_dir = Path(bench_dir)
    job_ids   = []

    for rep in bench_mols:
        mol_id   = rep["mol_id"]
        mol_path = bench_dir / "QCxMS2" / mol_id
        if not mol_path.is_dir():
            print(f"{mol_id}: directory not found, skipping.")
            continue
        if (mol_path / "crest_best.xyz").exists():
            print(f"{mol_id}: CREST already done, skipping.")
            continue
        print(f"{mol_id}: Submitting CREST ...")
        result = subprocess.run(
            ["sbatch", crest_script],
            cwd=str(mol_path),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        print(f"  {result.stdout.strip() or result.stderr.strip()}")
        m = re.search(r"\b(\d{5,})\b", result.stdout)
        if m:
            job_ids.append(m.group(1))

    return job_ids


def check_and_submit_qcxms2(bench_dir, reps, sizes, crest_script, qcxms2_script,
                              submit=False):
    """Check QCxMS2 status for benchmark molecules and optionally submit missing jobs.

    For each molecule:
      - No crest_best.xyz           → submit CREST
      - crest_best.xyz, no peaks.csv → submit QCxMS2
        (script auto-detects getieeab error and retries with -edist gaussian)

    Parameters
    ----------
    bench_dir : str
        Root memory_benchmark directory.
    reps : dict[str, list[dict]]
        Molecule representatives grouped by size label.
    sizes : tuple[str]
        Size labels in display order, e.g. ("small", "medium", "large").
    crest_script : str
        Absolute path to submit_batch_crest.sh.
    qcxms2_script : str
        Absolute path to submit_qcxms2_job.sh.
    submit : bool
        If False, only print status (dry run). If True, submit missing jobs.
    """
    from pathlib import Path
    bench_dir = Path(bench_dir)

    print(f"{'Mol':<6} {'Size':<8} {'crest_best':>10} {'peaks.csv':>10}  Action")
    print("-" * 55)

    for size in sizes:
        for rep in reps[size]:
            mol_id    = rep["mol_id"]
            mol_path  = bench_dir / "QCxMS2" / mol_id
            crest_done = False
            for log_name in ("crest.log", "crest_restart.log"):
                log_path = mol_path / log_name
                if log_path.exists() and "CREST terminated normally." in log_path.read_text(errors="replace"):
                    crest_done = True
                    break
            has_crest = crest_done
            _gfn2 = mol_path / "qcxms2_gfn2"
            has_peaks = (
                (_gfn2 / "qcxms2_auxiliary.tar.gz").exists()
                or (
                    (_gfn2 / "qcxms2.log").exists()
                    and "QCxMS2 terminated normally" in (_gfn2 / "qcxms2.log").read_text(errors="replace")
                )
            )
            is_perm_fail = (_gfn2 / "QCXMS2_PERMANENT_FAIL").exists()

            if is_perm_fail:
                action = "PERM_FAIL"
            elif has_peaks:
                action = "done"
            elif has_crest:
                action = "submit QCxMS2"
            else:
                action = "submit CREST"

            print(f"{mol_id:<6} {size:<8} {str(has_crest):>10} {str(has_peaks):>10}  {action}")

            if not submit or has_peaks or is_perm_fail:
                continue

            if has_crest:
                result = subprocess.run(
                    ["sbatch", qcxms2_script, str(mol_path), "gfn2"],
                    cwd=str(mol_path),
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True,
                )
                print(f"  → {result.stdout.strip() or result.stderr.strip()}")
            else:
                result = subprocess.run(
                    ["sbatch", crest_script],
                    cwd=str(mol_path),
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True,
                )
                print(f"  → {result.stdout.strip() or result.stderr.strip()}")


def check_stuck_molecules(datasets_stuck, min_size_kb=10):
    """Report file presence and restart file sizes for stuck molecules.

    Parameters
    ----------
    datasets_stuck : list[tuple[str, str, list[str]]]
        Each entry: (dataset_name, qcxms2_dir, list_of_stuck_folders).
    min_size_kb : int
        Minimum expected restart file size in KB to be considered valid.
    """
    min_size = min_size_kb * 1024
    for dataset_name, dataset, stuck in datasets_stuck:
        print(f"\n{'='*40}")
        print(f"{dataset_name} — stuck molecules")
        print(f"{'='*40}")
        for mol in stuck:
            mol_path = os.path.join(dataset, mol)
            files = set(os.listdir(mol_path))
            err_tail = ""
            err_path = os.path.join(mol_path, "crest.err")
            if os.path.exists(err_path):
                with open(err_path) as f:
                    lines = f.readlines()
                    err_tail = lines[-1].strip() if lines else ""
            print(
                f"  {mol}: restart={'crest.restart' in files} "
                f"backup={'crest.restart.bak' in files} "
                f"xtbopt={'xtbopt.sdf' in files} "
                f"log={'crest.log' in files} | err: {err_tail}"
            )
            for fname in ("crest.restart", "crest.restart.bak"):
                fpath = os.path.join(mol_path, fname)
                if os.path.exists(fpath):
                    size = os.path.getsize(fpath)
                    status = "OK" if size >= min_size else "CORRUPTED"
                    print(f"    {fname}: {size/1024:.1f} KB — {status}")
                else:
                    print(f"    {fname}: MISSING")
