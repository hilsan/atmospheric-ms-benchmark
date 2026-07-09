"""
Parse CPU/memory resource usage from SLURM accounting (sacct) and
method-specific log files for the atmospheric-MS-benchmark methods.

Strategy per method
-------------------
QCxMS GS-MD   : sacct matched by JobName + WorkDir; log in {variant}/GS_*.out
QCxMS frag    : sacct matched by JobName; log in {mol}/GS-opt/MS-run/TMPQCXMS/frag_*.out
CREST         : log in {mol}/crest.log or crest_restart.log  +  sacct by WorkDir
QCxMS2        : log in {mol}/qcxms2_gfn2/qcxms2.log         +  sacct by WorkDir
NEIMS         : SLURM .out in NEIMS base dir (has `time` output); sacct by ArrayJobID
CFMID         : SLURM .out in CFMID dir; sacct by JobID

All timing functions return values in seconds; memory in MB.
"""

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional


# ─────────────────────────────────────────────────────────────────────────────
# sacct helpers
# ─────────────────────────────────────────────────────────────────────────────

SACCT_FMT = (
    "JobID,JobName,WorkDir,Elapsed,TotalCPU,MaxRSS,AllocCPUS,AllocTRES,State"
)

# Module-level sacct cache — load once via load_sacct_cache().
_SACCT_CACHE = None  # {"allocation_records": [...], "batch_step_records": [...]}


def load_sacct_cache(cache_path):
    # type: (str) -> None
    """Load a sacct JSON cache file so queries fall back to it when sacct
    returns no results (e.g. after accounting records have been purged)."""
    global _SACCT_CACHE
    with open(cache_path) as f:
        _SACCT_CACHE = json.load(f)
    print(f"Loaded sacct cache: {len(_SACCT_CACHE.get('allocation_records', []))} "
          f"allocation + {len(_SACCT_CACHE.get('batch_step_records', []))} batch-step records")


def save_sacct_cache(job_ids_str, cache_path):
    # type: (str, str) -> None
    """Query sacct for the given comma-separated job IDs and save results to a
    JSON cache file.  Call this right after jobs complete, before records expire."""
    fmt = "JobID,JobName,WorkDir,Partition,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ExitCode"
    r1 = subprocess.run(
        ["sacct", "--format", fmt, "--parsable2", "--noheader", "--allocations", "-j", job_ids_str],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alloc_rows = []
    for line in r1.stdout.decode().splitlines():
        p = line.split("|")
        if len(p) >= 10:
            alloc_rows.append(dict(zip(fmt.split(","), p)))
    r2 = subprocess.run(
        ["sacct", "--format", "JobID,MaxRSS,Elapsed,State", "--parsable2", "--noheader", "-j", job_ids_str],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    step_rows = []
    for line in r2.stdout.decode().splitlines():
        p = line.split("|")
        if len(p) >= 4 and p[0].strip().endswith(".batch"):
            step_rows.append({"JobID": p[0], "MaxRSS": p[1], "Elapsed": p[2], "State": p[3]})
    import datetime
    cache = {"allocation_records": alloc_rows, "batch_step_records": step_rows,
             "saved": datetime.datetime.now().isoformat(), "job_ids": job_ids_str}
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2)
    print(f"Saved {len(alloc_rows)} allocation + {len(step_rows)} batch-step records → {cache_path}")


def _elapsed_to_seconds(elapsed):
    # type: (str) -> float
    """Convert sacct Elapsed string (D-HH:MM:SS or HH:MM:SS) to seconds."""
    elapsed = elapsed.strip()
    if not elapsed or elapsed in ("-", "00:00:00"):
        return 0.0
    try:
        if "-" in elapsed:
            days, rest = elapsed.split("-", 1)
            d = int(days)
        else:
            d, rest = 0, elapsed
        parts = rest.split(":")
        h, m, s = int(parts[0]), int(parts[1]), float(parts[2])
        return d * 86400 + h * 3600 + m * 60 + s
    except Exception:
        return 0.0


def _rss_to_mb(rss):
    # type: (str) -> float
    """Convert sacct MaxRSS string (e.g. '12345K', '1.2G') to MB."""
    rss = rss.strip()
    if not rss or rss in ("-", ""):
        return 0.0
    try:
        if rss.endswith("K"):
            return float(rss[:-1]) / 1024
        elif rss.endswith("M"):
            return float(rss[:-1])
        elif rss.endswith("G"):
            return float(rss[:-1]) * 1024
        else:
            return float(rss) / 1024  # assume KB
    except Exception:
        return 0.0


def _parse_tres_gpu(tres):
    # type: (str) -> int
    """Extract GPU count from AllocTRES string (e.g. 'gpu=1')."""
    if not tres:
        return 0
    match = re.search(r"gres/gpu=(\d+)", tres)
    return int(match.group(1)) if match else 0


def sacct_query(filters, user=None):
    # type: (List[str], Optional[str]) -> List[Dict]
    """Run sacct with given filter flags and return parsed rows.

    Falls back to the in-memory cache loaded by load_sacct_cache() when
    sacct returns no results (accounting records purged after ~1 year).

    Parameters
    ----------
    filters : list of str
        Extra sacct flags, e.g. ["-j", "12345"] or ["--name=crest"].
    user : str, optional
        Limit to this username.
    """
    cmd = ["sacct", "--format", SACCT_FMT, "--parsable2", "--noheader",
           "--allocations"]
    if user:
        cmd += ["-u", user]
    cmd += filters
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    rows = []
    if result.returncode == 0:
        for line in result.stdout.decode().splitlines():
            parts = line.split("|")
            if len(parts) < 9:
                continue
            rows.append({
                "JobID":     parts[0],
                "JobName":   parts[1],
                "WorkDir":   parts[2],
                "Elapsed":   parts[3],
                "TotalCPU":  parts[4],
                "MaxRSS":    parts[5],
                "AllocCPUS": parts[6],
                "AllocTRES": parts[7],
                "State":     parts[8],
            })

    if rows:
        return rows

    # Fallback: search the in-memory cache
    if _SACCT_CACHE is None:
        return []

    # Extract job IDs from -j filter
    job_filter = set()
    for i, f in enumerate(filters):
        if f == "-j" and i + 1 < len(filters):
            for jid in filters[i + 1].split(","):
                job_filter.add(jid.strip())

    cache_rows = []
    for rec in _SACCT_CACHE.get("allocation_records", []):
        jid = rec.get("JobID", "")
        # Match exact job ID or array subjob (e.g. "12345_3" matches filter "12345_3")
        if not job_filter or jid in job_filter or jid.split(".")[0] in job_filter:
            cache_rows.append({
                "JobID":     rec.get("JobID", ""),
                "JobName":   rec.get("JobName", ""),
                "WorkDir":   rec.get("WorkDir", ""),
                "Elapsed":   rec.get("Elapsed", ""),
                "TotalCPU":  rec.get("TotalCPU", ""),
                "MaxRSS":    rec.get("MaxRSS", ""),
                "AllocCPUS": rec.get("AllocCPUS", ""),
                "AllocTRES": rec.get("AllocTRES", ""),
                "State":     rec.get("State", ""),
            })
    return cache_rows


def sacct_by_workdir(workdir, user=None, start="2024-01-01"):
    # type: (str, Optional[str], str) -> List[Dict]
    """Return sacct rows whose WorkDir matches workdir (exact prefix match)."""
    workdir = str(Path(workdir).resolve())
    rows = sacct_query(["--starttime", start], user=user)
    return [r for r in rows if r["WorkDir"].startswith(workdir)]


def sacct_by_jobid(job_id):
    # type: (str) -> Optional[Dict]
    """Return sacct row for a specific job ID, or None."""
    rows = sacct_query(["-j", str(job_id)])
    return rows[0] if rows else None


def sacct_maxrss(job_id):
    # type: (str) -> float
    """Get peak RSS (MB) for a job from the .batch step record.

    Unlike sacct_query (which uses --allocations), this queries step records
    where MaxRSS is actually populated.  Falls back to the cache loaded by
    load_sacct_cache() when sacct returns no results.
    """
    cmd = ["sacct", "--format", "JobID,MaxRSS", "--parsable2", "--noheader",
           "-j", str(job_id)]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        for line in result.stdout.decode().splitlines():
            parts = line.split("|")
            if len(parts) >= 2 and parts[0].strip().endswith(".batch"):
                val = _rss_to_mb(parts[1])
                if val:
                    return val

    # Fallback: search batch-step records in cache
    if _SACCT_CACHE is not None:
        batch_id = str(job_id) + ".batch"
        for rec in _SACCT_CACHE.get("batch_step_records", []):
            if rec.get("JobID", "").strip() == batch_id:
                return _rss_to_mb(rec.get("MaxRSS", ""))
    return 0.0


def sacct_summary(rows):
    # type: (List[Dict]) -> Dict
    """Aggregate a list of sacct rows (e.g. array subjobs) into one summary."""
    if not rows:
        return {}
    wall_s   = max(_elapsed_to_seconds(r["Elapsed"]) for r in rows)
    cpu_s    = sum(_elapsed_to_seconds(r["TotalCPU"]) for r in rows)
    peak_mb  = max(_rss_to_mb(r["MaxRSS"]) for r in rows)
    n_cpus   = int(rows[0].get("AllocCPUS", 0) or 0)
    n_gpus   = _parse_tres_gpu(rows[0].get("AllocTRES", ""))
    state    = rows[0].get("State", "")
    return {
        "wall_s":  wall_s,
        "cpu_s":   cpu_s,
        "peak_mb": peak_mb,
        "n_cpus":  n_cpus,
        "n_gpus":  n_gpus,
        "state":   state,
        "n_jobs":  len(rows),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Log file parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_crest_timing(mol_dir):
    # type: (str) -> Dict
    """Parse wall time from CREST log in mol_dir.

    CREST 3.x prints:
      "CREST runtime (total)               0 d,  2 h, 15 min, 21.833 sec"
    followed by:
      " * wall-time:     0 d,  2 h, 15 min, 21.833 sec"
    """
    mol_dir = Path(mol_dir)
    time_re = re.compile(
        r"(\d+)\s*d,\s*(\d+)\s*h,\s*(\d+)\s*min,\s*([\d.]+)\s*sec"
    )
    for log_name in ("crest.log", "crest_restart.log"):
        log_path = mol_dir / log_name
        if not log_path.exists():
            continue
        text = log_path.read_text(errors="replace")
        # CREST 3.x: "CREST runtime (total)   0 d,  2 h, 15 min, 21.833 sec"
        match = re.search(
            r"CREST runtime \(total\)\s+" + time_re.pattern, text
        )
        if match:
            d, h, m, s = int(match.group(1)), int(match.group(2)), int(match.group(3)), float(match.group(4))
            wall_s = d * 86400 + h * 3600 + m * 60 + s
            return {"wall_s": wall_s, "source": log_name}
        # Older CREST: "Total run time:  0 d,  0 h, 12 m, 34.2 s"
        match = re.search(
            r"Total run time:\s*(\d+)\s*d,\s*(\d+)\s*h,\s*(\d+)\s*m,\s*([\d.]+)\s*s",
            text,
        )
        if match:
            d, h, m, s = int(match.group(1)), int(match.group(2)), int(match.group(3)), float(match.group(4))
            wall_s = d * 86400 + h * 3600 + m * 60 + s
            return {"wall_s": wall_s, "source": log_name}
    return {}


def parse_qcxms2_timing(mol_dir, method="gfn2"):
    # type: (str, str) -> Dict
    """Parse wall time from qcxms2.log.

    QCxMS2 typically prints timing near the end of the log.
    Falls back to the SLURM output file qcxms2_*.out in mol_dir.
    """
    mol_dir = Path(mol_dir)
    log_path = mol_dir / ("qcxms2_" + method) / "qcxms2.log"
    if log_path.exists():
        text = log_path.read_text(errors="replace")
        # QCxMS2 prints "Overall wall time  : 19h :33m :28s"
        match = re.search(
            r"Overall wall time\s*:\s*(\d+)h\s*:(\d+)m\s*:(\d+)s",
            text,
        )
        if match:
            wall_s = int(match.group(1)) * 3600 + int(match.group(2)) * 60 + int(match.group(3))
            return {"wall_s": wall_s, "source": "qcxms2.log"}
        match = re.search(
            r"[Tt]otal\s+run\s+time[:\s]+(\d+)\s*h[,\s]*(\d+)\s*min[,\s]*([\d.]+)\s*sec",
            text,
        )
        if match:
            wall_s = int(match.group(1)) * 3600 + int(match.group(2)) * 60 + float(match.group(3))
            return {"wall_s": wall_s, "source": "qcxms2.log"}
        match = re.search(
            r"(\d+)\s*d[,\s]*(\d+)\s*h[,\s]*(\d+)\s*min[,\s]*([\d.]+)\s*sec",
            text,
        )
        if match:
            wall_s = (int(match.group(1)) * 86400 + int(match.group(2)) * 3600
                      + int(match.group(3)) * 60 + float(match.group(4)))
            return {"wall_s": wall_s, "source": "qcxms2.log"}

    return parse_slurm_time_output(mol_dir, prefix="qcxms2_")


def parse_slurm_time_output(search_dir, prefix=""):
    # type: (str, str) -> Dict
    """Parse wall/user/sys time from a SLURM .out file containing `time` output.

    The NEIMS script wraps the call in `{ time python ... } 2>&1`, producing:
        real    0m12.345s
        user    0m10.123s
        sys     0m1.234s
    """
    search_dir = Path(search_dir)
    pattern = prefix + "*.out" if prefix else "*.out"
    candidates = sorted(search_dir.glob(pattern))
    for out_file in candidates:
        text = out_file.read_text(errors="replace")
        match = re.search(r"^real\s+(\d+)m([\d.]+)s", text, re.MULTILINE)
        if match:
            wall_s = int(match.group(1)) * 60 + float(match.group(2))
            return {"wall_s": wall_s, "source": out_file.name}
    return {}


def parse_qcxms_gs_slurm_output(variant_dir, mol_index):
    # type: (str, int) -> Dict
    """Find and parse the GS-MD SLURM output file for a molecule.

    Output file pattern: GS_{arrayJobId}_{mol_index}.out in variant_dir.
    Returns dict with 'slurm_out' path and 'array_job_id' if found.
    """
    variant_dir = Path(variant_dir)
    candidates = sorted(variant_dir.glob("GS_*_{}.out".format(mol_index)))
    for out_file in candidates:
        match = re.search(r"GS_(\d+)_\d+\.out", out_file.name)
        job_id = match.group(1) if match else None
        return {"slurm_out": str(out_file), "array_job_id": job_id}
    return {}


def find_slurm_jobid_in_dir(directory, prefix=""):
    # type: (str, str) -> Optional[str]
    """Find the first SLURM output file in directory and extract the job ID."""
    directory = Path(directory)
    pattern = prefix + "*.out" if prefix else "*.out"
    for out_file in sorted(directory.glob(pattern)):
        match = re.search(r"(\d{5,})", out_file.stem)
        if match:
            return match.group(1)
    return None


def find_qcxms2_jobid(mol_dir):
    # type: (str) -> Optional[str]
    """Extract SLURM job ID from qcxms2_{jobid}_*.out files in mol_dir."""
    mol_dir = Path(mol_dir)
    for out_file in sorted(mol_dir.glob("qcxms2_*.out")):
        match = re.search(r"qcxms2_(\d+)", out_file.name)
        if match:
            return match.group(1)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# High-level per-molecule collectors
# ─────────────────────────────────────────────────────────────────────────────

def collect_qcxms_resources(sim_base_dir, mol_id, variant, user=None):
    # type: (str, str, str, Optional[str]) -> Dict
    """Collect GS-MD + fragmentation resource usage for one QCxMS molecule.

    Wall time and CPU time come from sacct (total job elapsed, not just xTB
    opt).  MaxRSS is fetched separately from the .batch step record.
    """
    mol_dir     = Path(sim_base_dir) / variant / mol_id
    gsopt_dir   = mol_dir / "GS-opt"
    msrun_dir   = gsopt_dir / "MS-run"
    variant_dir = Path(sim_base_dir) / variant

    result = {"variant": variant, "mol_id": mol_id, "method": "qcxms"}

    # ── GS-MD ────────────────────────────────────────────────────────────────
    # Use sacct for total job time (includes xTB opt + full QCxMS GS-MD run).
    # opt.out only has the xTB geometry optimisation — it misses the GS-MD.
    gs_result = {}
    mol_index = int(mol_id)
    gs_info = parse_qcxms_gs_slurm_output(str(variant_dir), mol_index)
    if gs_info.get("array_job_id"):
        job_spec = "{}_{}".format(gs_info["array_job_id"], mol_index)
        rows = sacct_query(["-j", job_spec], user=user)
        if rows:
            gs_result = sacct_summary(rows)
            gs_result["source"] = "sacct"
        peak_mb = sacct_maxrss(job_spec)
        if peak_mb:
            gs_result["peak_mb"] = peak_mb

    # Fallback: sacct by WorkDir
    if not gs_result and gsopt_dir.exists():
        rows = sacct_by_workdir(str(gsopt_dir), user=user)
        if rows:
            gs_result = sacct_summary(rows)
            gs_result["source"] = "sacct_workdir"
            job_id = rows[0].get("JobID", "").split("_")[0]
            if job_id:
                peak_mb = sacct_maxrss(job_id)
                if peak_mb:
                    gs_result["peak_mb"] = peak_mb

    result["gsmd"] = gs_result

    # ── Fragmentation ────────────────────────────────────────────────────────
    frag_result = {}
    tmpqcxms_dir = msrun_dir / "TMPQCXMS"
    if tmpqcxms_dir.exists():
        job_id = find_slurm_jobid_in_dir(str(tmpqcxms_dir), prefix="frag_")
        if job_id:
            rows = sacct_query(["-j", job_id], user=user)
            if rows:
                frag_result = sacct_summary(rows)
                frag_result["source"] = "sacct"
            peak_mb = sacct_maxrss(job_id)
            if peak_mb:
                frag_result["peak_mb"] = peak_mb

    # Fallback: parse tmpqcxms.out for wall time
    if not frag_result:
        tmpqcxms_out = msrun_dir / "tmpqcxms.out"
        if tmpqcxms_out.exists():
            text = tmpqcxms_out.read_text(errors="replace")
            matches = re.findall(r"wall time \(min\)\s+([\d.]+)", text)
            if matches:
                frag_result["wall_s"] = float(matches[-1]) * 60
                frag_result["source"] = "tmpqcxms.out"

    result["frag"] = frag_result
    return result


def collect_crest_resources(sim_base_dir, mol_id, user=None):
    # type: (str, str, Optional[str]) -> Dict
    """Collect CREST resource usage for one molecule."""
    mol_dir = Path(sim_base_dir) / "QCxMS2" / mol_id
    result  = {"mol_id": mol_id, "method": "crest"}

    # Log file timing (most reliable for wall time)
    timing = parse_crest_timing(str(mol_dir))
    result.update(timing)

    # sacct for CPU time and memory
    rows = sacct_by_workdir(str(mol_dir), user=user)
    crest_rows = [r for r in rows if "crest" in r["JobName"].lower()]
    if crest_rows:
        summ = sacct_summary(crest_rows)
        # Keep log-derived wall_s if available (more accurate than sacct elapsed)
        if "wall_s" not in result:
            result["wall_s"] = summ.get("wall_s", 0)
        result.update({k: v for k, v in summ.items() if k not in result})
        # Get MaxRSS from step record
        job_id = crest_rows[0].get("JobID", "").split(".")[0]
        if job_id:
            peak_mb = sacct_maxrss(job_id)
            if peak_mb:
                result["peak_mb"] = peak_mb

    return result


def collect_qcxms2_resources(sim_base_dir, mol_id, method="gfn2", user=None):
    # type: (str, str, str, Optional[str]) -> Dict
    """Collect QCxMS2 resource usage for one molecule."""
    mol_dir = Path(sim_base_dir) / "QCxMS2" / mol_id
    result  = {"mol_id": mol_id, "method": "qcxms2"}

    # Log file timing
    timing = parse_qcxms2_timing(str(mol_dir), method=method)
    result.update(timing)

    # Try to get job ID directly from output filename
    job_id = find_qcxms2_jobid(str(mol_dir))
    if job_id:
        rows = sacct_query(["-j", job_id], user=user)
        if rows:
            summ = sacct_summary(rows)
            if "wall_s" not in result:
                result["wall_s"] = summ.get("wall_s", 0)
            result.update({k: v for k, v in summ.items() if k not in result})
        peak_mb = sacct_maxrss(job_id)
        if peak_mb:
            result["peak_mb"] = peak_mb
    else:
        # Fallback: sacct by WorkDir
        rows = sacct_by_workdir(str(mol_dir), user=user)
        qcxms2_rows = [r for r in rows if "qcxms2" in r["JobName"].lower()]
        if qcxms2_rows:
            summ = sacct_summary(qcxms2_rows)
            if "wall_s" not in result:
                result["wall_s"] = summ.get("wall_s", 0)
            result.update({k: v for k, v in summ.items() if k not in result})
            job_id = qcxms2_rows[0].get("JobID", "").split(".")[0]
            if job_id:
                peak_mb = sacct_maxrss(job_id)
                if peak_mb:
                    result["peak_mb"] = peak_mb

    return result


def collect_neims_resources(sim_base_dir, mol_id, user=None):
    # type: (str, str, Optional[str]) -> Dict
    """Collect NEIMS resource usage for one molecule."""
    neims_dir = Path(sim_base_dir) / "NEIMS"
    mol_index = int(mol_id)
    result    = {"mol_id": mol_id, "method": "neims"}

    # NEIMS output files are in the base NEIMS dir and named {job}_{array}_{task}.out
    candidates = sorted(neims_dir.glob("*_*_{}.out".format(mol_index)))
    for out_file in candidates:
        text = out_file.read_text(errors="replace")
        match = re.search(r"^real\s+(\d+)m([\d.]+)s", text, re.MULTILINE)
        if match:
            result["wall_s"] = int(match.group(1)) * 60 + float(match.group(2))
            result["source"] = out_file.name

            id_match = re.match(r"(\d+)_(\d+)_\d+\.out", out_file.name)
            if id_match:
                job_id  = "{}_{}".format(id_match.group(1), mol_index)
                rows    = sacct_query(["-j", job_id], user=user)
                if rows:
                    summ = sacct_summary(rows)
                    result.update({k: v for k, v in summ.items()
                                   if k not in ("wall_s",)})
                peak_mb = sacct_maxrss(job_id)
                if peak_mb:
                    result["peak_mb"] = peak_mb
            break

    if "wall_s" not in result:
        rows = sacct_by_workdir(str(neims_dir), user=user)
        neims_rows = [r for r in rows if "neims" in r["JobName"].lower()]
        if neims_rows:
            result.update(sacct_summary(neims_rows))

    return result


def collect_cfmid_resources(sim_base_dir, user=None):
    # type: (str, Optional[str]) -> Dict
    """Collect CFMID resource usage (runs all molecules as one job).

    The script runs two timed apptainer calls (non-iso + iso) whose `real`
    timing lines both appear in the same SLURM output file.  We sum them.
    """
    cfmid_dir = Path(sim_base_dir) / "CFMID"
    result    = {"method": "cfmid"}

    candidates = sorted(cfmid_dir.glob("cfmid_*.out"))
    for out_file in candidates:
        text = out_file.read_text(errors="replace")
        matches = re.findall(r"^real\s+(\d+)m([\d.]+)s", text, re.MULTILINE)
        if matches:
            total = sum(int(m[0]) * 60 + float(m[1]) for m in matches)
            result.update({"wall_s": total, "source": out_file.name})
            break

    job_id = find_slurm_jobid_in_dir(str(cfmid_dir), prefix="cfmid_")
    if job_id:
        rows = sacct_query(["-j", job_id], user=user)
        if rows:
            summ = sacct_summary(rows)
            if "wall_s" not in result:
                result["wall_s"] = summ.get("wall_s", 0)
            result.update({k: v for k, v in summ.items() if k not in result})
        peak_mb = sacct_maxrss(job_id)
        if peak_mb:
            result["peak_mb"] = peak_mb

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Molecule size helpers
# ─────────────────────────────────────────────────────────────────────────────

def heavy_atom_count(smiles):
    # type: (str) -> int
    """Count non-hydrogen atoms in a SMILES string (rdkit if available, else regex)."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return mol.GetNumHeavyAtoms()
    except ImportError:
        pass
    smiles_noh = re.sub(r"\[H\]", "", smiles)
    return len(re.findall(r"[BCNOFPSIKCL]", smiles_noh, re.IGNORECASE))


def select_size_representatives(csv_path, smiles_col="SMILES", id_col=None):
    # type: (str, str, Optional[str]) -> Dict
    """Pick small/medium/large molecules from a processed CSV."""
    import pandas as pd

    df = pd.read_csv(csv_path)

    if smiles_col not in df.columns:
        for col in ["SMILES", "Modified_SMILES", "smiles", "Smiles"]:
            if col in df.columns:
                smiles_col = col
                break

    df["_ha"] = df[smiles_col].apply(
        lambda s: heavy_atom_count(str(s)) if pd.notna(s) else 0
    )
    df = df[df["_ha"] > 0].reset_index()
    df["_mol_id"] = df["index"].apply(lambda i: "{:04d}".format(i))

    df_sorted = df.sort_values("_ha")
    n         = len(df_sorted)

    small  = df_sorted.iloc[n // 6]
    medium = df_sorted.iloc[n // 2]
    large  = df_sorted.iloc[5 * n // 6]

    result = {}
    for label, row in [("small", small), ("medium", medium), ("large", large)]:
        result[label] = {
            "mol_id":      row["_mol_id"],
            "smiles":      row[smiles_col],
            "heavy_atoms": int(row["_ha"]),
        }

    return result


def select_size_representatives_multi(csv_path, n=3, smiles_col="SMILES",
                                       id_col=None):
    # type: (str, int, str, Optional[str]) -> Dict
    """Pick n small/medium/large molecules from a processed CSV."""
    import pandas as pd

    df = pd.read_csv(csv_path)

    if smiles_col not in df.columns:
        for col in ["SMILES", "Modified_SMILES", "smiles", "Smiles"]:
            if col in df.columns:
                smiles_col = col
                break

    df["_ha"] = df[smiles_col].apply(
        lambda s: heavy_atom_count(str(s)) if pd.notna(s) else 0
    )
    df = df[df["_ha"] > 0].reset_index()
    df["_mol_id"] = df["index"].apply(lambda i: "{:04d}".format(i))

    df_sorted = df.sort_values("_ha").reset_index(drop=True)
    total = len(df_sorted)
    t = total // 3

    terciles = {
        "small":  df_sorted.iloc[:t].reset_index(drop=True),
        "medium": df_sorted.iloc[t:2 * t].reset_index(drop=True),
        "large":  df_sorted.iloc[2 * t:].reset_index(drop=True),
    }

    result = {}
    for label, terc in terciles.items():
        tn = len(terc)
        indices = [int((i + 0.5) * tn / n) for i in range(n)]
        result[label] = [
            {
                "mol_id":      terc.iloc[i]["_mol_id"],
                "smiles":      terc.iloc[i][smiles_col],
                "heavy_atoms": int(terc.iloc[i]["_ha"]),
            }
            for i in indices
        ]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Formatting helpers
# ─────────────────────────────────────────────────────────────────────────────

def fmt_wall(seconds):
    # type: (float) -> str
    if not seconds:
        return u"\u2014"
    h, rem = divmod(int(seconds), 3600)
    m, s   = divmod(rem, 60)
    if h:
        return "{}h {:02d}m".format(h, m)
    return "{}m {:02d}s".format(m, s)


def fmt_mem(mb):
    # type: (float) -> str
    if not mb:
        return u"\u2014"
    return "{:.1f} GB".format(mb / 1024) if mb >= 1024 else "{:.0f} MB".format(mb)
