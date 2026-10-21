#!/usr/bin/env python3
"""genomeforge — command-line front door for the GenomeForge pipeline.

  genomeforge run           launch (or resume) the pipeline on the cluster
  genomeforge samplesheet   build the lane-level sample sheet from FASTQ directories
  genomeforge check         verify a site's references, containers and tools before a run
  genomeforge resources     stage the VQSR training sets and GATK calling regions
  genomeforge status        show the state of a run directory
  genomeforge qc-summary    rebuild the per-sample QC table from a results directory
  genomeforge test          dry-run the whole DAG on a synthetic layout (install check)

The CLI is stdlib-only. Every analytical tool runs inside the CallForge containers that
Nextflow starts per task; nothing is installed on the host.
"""
import argparse
import glob
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

try:
    from . import __version__
except ImportError:                                   # run as a plain script (share/bin shim)
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from genomeforge_cli import __version__

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.environ.get("GENOMEFORGE_HOME") or os.path.dirname(PKG_DIR)
BIN = os.path.join(REPO_ROOT, "bin")
PROG = "genomeforge"

ENV_HELP = """\
environment (all optional; command-line flags override them):
  GENOMEFORGE_HOME        pipeline install directory (default: where this CLI lives)
  GENOMEFORGE_PROFILE     site profile, e.g. noguchi or slurm (default: slurm)
  GENOMEFORGE_PARAMS      -params-file with the site paths (needed with the generic slurm profile)
  GENOMEFORGE_RUNDIR      run directory holding work/, tmp/, logs/, results/ (default: ./genomeforge_run)
  GENOMEFORGE_PARTITION   SLURM partition for the head job (default: the cluster default)
  GENOMEFORGE_EXCLUDE_NODES  nodes to keep free, e.g. node4
"""


# ----------------------------------------------------------------------------- helpers
def _say(msg):
    sys.stderr.write(f"[{PROG}] {msg}\n")


def _die(msg, code=1):
    _say(f"ERROR: {msg}")
    sys.exit(code)


def _load_bin(script):
    path = os.path.join(BIN, script)
    if not os.path.isfile(path):
        _die(f"cannot find {path} (is GENOMEFORGE_HOME set to the pipeline directory?)", 2)
    spec = importlib.util.spec_from_file_location(script.replace(".py", ""), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_bin(script, prog, rest):
    mod = _load_bin(script)
    old = sys.argv
    sys.argv = [prog, *rest]
    try:
        return mod.main() or 0
    finally:
        sys.argv = old


def _have(tool):
    return shutil.which(tool) is not None


def _site_defaults():
    return {
        "profile":   os.environ.get("GENOMEFORGE_PROFILE") or "slurm",
        "params":    os.environ.get("GENOMEFORGE_PARAMS") or None,
        "rundir":    os.environ.get("GENOMEFORGE_RUNDIR") or os.path.join(os.getcwd(), "genomeforge_run"),
        "partition": os.environ.get("GENOMEFORGE_PARTITION") or None,
    }


def _nextflow_params(profile, params_file, overrides):
    """Resolve params as Nextflow would: profile config -> params file -> --key value."""
    params = {}
    if _have("nextflow"):
        try:
            out = subprocess.run(["nextflow", "config", REPO_ROOT, "-flat", "-profile", profile],
                                 capture_output=True, text=True, timeout=120,
                                 env={**os.environ, "NXF_OFFLINE": "true", "NXF_ANSI_LOG": "false"}).stdout
            for line in out.splitlines():
                if line.startswith("params."):
                    k, _, v = line[7:].partition(" = ")
                    v = v.strip()
                    if v.startswith(("'", '"')) and v.endswith(("'", '"')):
                        v = v[1:-1]
                    params[k.strip()] = None if v == "null" else v
        except (subprocess.SubprocessError, OSError) as e:
            _say(f"WARN: could not resolve the {profile} profile with `nextflow config` ({e})")
    if params_file:
        params.update(_read_params_file(params_file))
    i = 0
    while i < len(overrides):
        a = overrides[i]
        if a.startswith("--") and i + 1 < len(overrides) and not overrides[i + 1].startswith("-"):
            params[a[2:]] = overrides[i + 1]; i += 2
        else:
            i += 1
    return params


def _read_params_file(path):
    if not os.path.isfile(path):
        _die(f"params file not found: {path}")
    try:
        import yaml                                    # present on most systems; optional
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        return {k: (None if v is None else str(v)) for k, v in data.items()}
    except ImportError:
        data = {}
        for line in open(path):
            line = line.split("#", 1)[0].rstrip()
            if ":" in line and not line.startswith((" ", "\t")):
                k, _, v = line.partition(":")
                v = v.strip().strip("'\"")
                data[k.strip()] = None if v in ("", "null", "~") else v
        return data


# ----------------------------------------------------------------------------- run
RUN_EPILOG = """\
examples:
  # Noguchi HPC (bundled site profile; nothing else to configure)
  genomeforge run --profile noguchi --input samplesheet.csv --rundir /hpc/data/<lab>/<project>/genomeforge_run

  # any SLURM cluster: generic profile + a params file with the site paths
  genomeforge run --profile slurm --params-file site.params.yaml --input samplesheet.csv

  # dry-run the whole DAG on the cluster first (seconds; no containers executed)
  genomeforge run --profile noguchi --input samplesheet.csv --stub

  # pilot: two or three samples, fewer shards, then the full cohort with the same command (-resume is implied)
  genomeforge run --profile noguchi --input pilot.csv --scatter_count 8 --outdir results_pilot

  # override any pipeline parameter inline (they pass straight through to Nextflow)
  genomeforge run --profile noguchi --input samplesheet.csv --filtering hard --archive_cram true --trim false

  # keep one node free for other users; the head job and every task stay off it
  genomeforge run --profile noguchi --input samplesheet.csv --exclude-nodes node4

  # show the command that would be submitted, without submitting it
  genomeforge run --profile noguchi --input samplesheet.csv --print-cmd

pipeline parameters (pass-through; defaults in nextflow.config, site paths from the profile/params file):
  --outdir DIR              results directory inside the run dir (default results)
  --scatter_count N         HaplotypeCaller / joint-genotyping shards (48)
  --calling_regions FILE    GATK wgs_calling_regions interval_list or BED (site profile)
  --filtering vqsr|hard     VQSR (default) or GATK hard filters; auto-fallback if training sets are missing
  --snp_tranche F / --indel_tranche F     VQSR truth-sensitivity levels (99.7 / 99.7)
  --trim true|false         fastp streamed into bwa-mem2 (true)
  --sequencing_center STR   read-group CN tag (unset)
  --run_vep / --run_pca / --run_contamination true|false     switch stages off
  --archive_cram true|false CRAM next to each final BAM (false)
  --qc_min_depth N  --qc_max_dup F  --qc_max_contam F  --qc_min_mapped F     QC flag thresholds
  --fasta / --dbsnp / --mills / --known_indels / --hapmap / --omni / --g1000 / --clinvar
  --gnomad_dir / --gnomad_pattern / --gnomad_fields / --somalier_sites / --vep_cache / --vep_release
  --callforge_sif / --annotate_sif / --cohortqc_sif / --burden_sif / --bind_paths / --slurm_partition
  any other Nextflow option (-resume is added for you; -w/-work-dir is set to <rundir>/work)

""" + ENV_HELP


def cmd_run(rest):
    d = _site_defaults()
    ap = argparse.ArgumentParser(
        prog=f"{PROG} run", allow_abbrev=False, formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Launch or resume the GenomeForge pipeline. The Nextflow head job is submitted with "
                     "sbatch (so it survives logout) and work/, tmp/, logs/ and results/ live under the run "
                     "directory on disk. Re-running the same command after a failure resumes it."),
        epilog=RUN_EPILOG)
    ap.add_argument("--input", metavar="CSV", help="sample sheet (one row per lane-pair; see `genomeforge samplesheet --help`)")
    ap.add_argument("--profile", default=d["profile"], help=f"site profile: noguchi, slurm, local, stub (default {d['profile']})")
    ap.add_argument("--params-file", default=d["params"], metavar="YAML", help="site parameter file (copy params.example.yaml)")
    ap.add_argument("--rundir", default=d["rundir"], metavar="DIR", help=f"run directory (default {d['rundir']})")
    ap.add_argument("--partition", default=d["partition"], help="SLURM partition for the head job (default: cluster default)")
    ap.add_argument("--exclude-nodes", default=os.environ.get("GENOMEFORGE_EXCLUDE_NODES"), metavar="LIST",
                    help="SLURM nodes to keep free for other users, e.g. node4 (head job and every task avoid them)")
    ap.add_argument("--max-restarts", type=int, default=20, help="self-resubmissions of the head job after an unfinished exit (default 20)")
    ap.add_argument("--restart-delay", type=int, default=600, help="seconds to wait before a self-resubmission (default 600)")
    ap.add_argument("--stub", action="store_true", help="dry-run the DAG with stub tasks (no containers executed)")
    ap.add_argument("--foreground", action="store_true", help="run Nextflow in this shell instead of an sbatch head job")
    ap.add_argument("--no-resume", action="store_true", help="start a fresh run instead of resuming")
    ap.add_argument("--print-cmd", action="store_true", help="print the nextflow command and exit")
    known, passthrough = ap.parse_known_args(rest)

    if not known.input and not any(a == "--input" for a in passthrough) and not known.print_cmd:
        ap.error("--input samplesheet.csv is required (build it with `genomeforge samplesheet`)")
    if known.profile == "slurm" and not known.params_file and not known.stub:
        _say("WARN: the generic slurm profile carries no site paths; pass --params-file site.params.yaml "
             "(or GENOMEFORGE_PARAMS) unless every path is given on the command line")

    rundir = os.path.abspath(known.rundir)
    for sub in ("work", "tmp", "logs"):
        os.makedirs(os.path.join(rundir, sub), exist_ok=True)
    if not os.access(rundir, os.W_OK):
        _die(f"run directory is not writable: {rundir}")

    main_nf = os.path.join(REPO_ROOT, "main.nf")
    if not os.path.isfile(main_nf):
        _die(f"main.nf not found under {REPO_ROOT}; set GENOMEFORGE_HOME to the pipeline directory", 2)
    cmd = ["nextflow", "-log", os.path.join(rundir, "nextflow.log"), "run", main_nf,
           "-profile", known.profile, "-work-dir", os.path.join(rundir, "work")]
    if known.params_file:
        cmd += ["-params-file", os.path.abspath(known.params_file)]
    if known.exclude_nodes and "--slurm_exclude" not in passthrough:
        cmd += ["--slurm_exclude", known.exclude_nodes]
    if known.input:
        cmd += ["--input", os.path.abspath(known.input)]
    if not known.no_resume:
        cmd.append("-resume")
    if known.stub:
        cmd.append("-stub-run")
    cmd += passthrough
    shown = " ".join(shlex.quote(c) for c in cmd)
    _say("+ " + shown)
    if known.print_cmd:
        print(shown)
        return 0

    env = {**os.environ,
           "NXF_WORK": os.path.join(rundir, "work"), "TMPDIR": os.path.join(rundir, "tmp"),
           "NXF_TEMP": os.path.join(rundir, "tmp"), "NXF_ANSI_LOG": "false", "NXF_OFFLINE": "true",
           "NXF_OPTS": f"-Xms1g -Xmx6g -Djava.io.tmpdir={os.path.join(rundir, 'tmp')}"}
    env.setdefault("NXF_SINGULARITY_CACHEDIR", os.path.expanduser("~/.genomeforge/container_cache"))

    if known.foreground or known.profile in ("local", "stub") or not _have("sbatch"):
        if not _have("nextflow"):
            _die("`nextflow` not found on PATH (needs Nextflow >= 24.10 and Java 17+)", 127)
        if not known.foreground and known.profile not in ("local", "stub"):
            _say("WARN: sbatch not found; running the head job in the foreground")
        return subprocess.call(cmd, cwd=rundir, env=env)

    stamp = time.strftime("%Y%m%d_%H%M%S")
    import hashlib
    jobname = "genomeforge-" + hashlib.sha1(rundir.encode()).hexdigest()[:8]     # one head chain per run directory
    head = os.path.join(rundir, "logs", f"head_{stamp}.sbatch")
    done = os.path.join(rundir, "logs", f"head_{stamp}.complete")
    exports = "\n".join(f"export {k}={shlex.quote(v)}" for k, v in env.items()
                        if k in ("NXF_WORK", "TMPDIR", "NXF_TEMP", "NXF_ANSI_LOG", "NXF_OFFLINE", "NXF_OPTS", "NXF_SINGULARITY_CACHEDIR"))
    sbatch_opts = [f"--job-name={jobname}", "--cpus-per-task=2", "--mem=8G", "--requeue", "--open-mode=append",
                   f"--output={rundir}/logs/head_{stamp}.out"]
    if known.partition:
        sbatch_opts.append(f"--partition={known.partition}")
    if known.exclude_nodes:
        sbatch_opts.append(f"--exclude={known.exclude_nodes}")
    with open(head, "w") as fh:
        fh.write("#!/bin/bash\n" + "".join(f"#SBATCH {o}\n" for o in sbatch_opts))
        fh.write(f"""# GenomeForge head job — self-resuming.
#  * SLURM --requeue restarts it after a node failure; the queue state survives a controller restart.
#  * On any unfinished exit (Nextflow crash, power trip while the node was up) it sleeps and
#    resubmits itself with -resume, up to GF_MAX attempts. A marker file ends the chain.
#  * Only one head job per run directory ever runs: a second instance exits immediately.
set -u
cd {shlex.quote(rundir)}
{exports}
HEAD={shlex.quote(head)}
DONE={shlex.quote(done)}
GF_ATTEMPT=${{GF_ATTEMPT:-1}}
GF_MAX={known.max_restarts}
GF_DELAY={known.restart_delay}
ulimit -n 65536 2>/dev/null || true
[ -f "$DONE" ] && {{ echo "run already completed ($DONE); nothing to do"; exit 0; }}
others=$(squeue -h -u "$USER" -n {jobname} -t RUNNING -o %i | grep -vx "${{SLURM_JOB_ID:-x}}" || true)
[ -z "$others" ] || {{ echo "another head job ($others) is already running for this run directory; exiting"; exit 0; }}
echo "=== GenomeForge head job $SLURM_JOB_ID on $(hostname) at $(date)  attempt $GF_ATTEMPT/$GF_MAX"
{shown}
rc=$?
echo "=== nextflow exited with status $rc at $(date)"
if [ $rc -eq 0 ]; then
    touch "$DONE"; echo "run complete"; exit 0
fi
if [ "$GF_ATTEMPT" -ge "$GF_MAX" ]; then
    echo "giving up after $GF_ATTEMPT attempts; fix the cause and run the same genomeforge command again (it resumes)"; exit $rc
fi
echo "unfinished: resubmitting in $GF_DELAY s as attempt $((GF_ATTEMPT + 1))/$GF_MAX (resume)"
sleep "$GF_DELAY"
sbatch --export=ALL,GF_ATTEMPT=$((GF_ATTEMPT + 1)) "$HEAD"
exit $rc
""")
    out = subprocess.run(["sbatch", "--parsable", head], capture_output=True, text=True)
    if out.returncode != 0:
        _die(f"sbatch failed: {out.stderr.strip()}")
    job = out.stdout.strip()
    print(f"submitted head job {job} ({jobname})\n  run dir   : {rundir}\n  head      : {head}\n"
          f"  nodes     : {'all except ' + known.exclude_nodes if known.exclude_nodes else 'all in the partition'}\n"
          f"  resilience: SLURM requeue + up to {known.max_restarts} self-resubmissions ({known.restart_delay}s apart), always -resume\n"
          f"  monitor   : genomeforge status --rundir {rundir}    |    tail -f {rundir}/logs/head_{stamp}.out")
    return 0


# ----------------------------------------------------------------------------- samplesheet
def cmd_samplesheet(rest):
    if not rest or rest in (["-h"], ["--help"]):
        mod = _load_bin("make_samplesheet.py")
        print(mod.__doc__.replace("make_samplesheet.py", f"{PROG} samplesheet"))
        print("examples:\n"
              f"  {PROG} samplesheet --roots /data/delivery-part1 /data/delivery-part2 --out samplesheet.csv\n"
              f"  {PROG} samplesheet --roots /data/delivery --metadata phenotypes.csv --exclude S1,S2 --out samplesheet.csv\n"
              f"  {PROG} samplesheet --fastq-csv pairs.csv --out samplesheet.csv\n"
              f"  {PROG} samplesheet --roots /data/delivery --sample-from filename --min-bytes 500000000 --out samplesheet.csv\n")
        return 0
    return _run_bin("make_samplesheet.py", f"{PROG} samplesheet", rest)


# ----------------------------------------------------------------------------- check
def cmd_check(rest):
    d = _site_defaults()
    ap = argparse.ArgumentParser(
        prog=f"{PROG} check", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Verify that a site is ready to run: tools on PATH, container images, reference FASTA "
                     "with bwa-mem2 index / .fai / .dict, every VCF with its index, gnomAD files, VEP cache, "
                     "and the contig style of the annotation files versus the reference. Nothing is run."),
        epilog=f"examples:\n  {PROG} check --profile noguchi\n  {PROG} check --profile slurm --params-file site.params.yaml\n"
               f"  {PROG} check --profile slurm --params-file site.params.yaml --fasta /refs/other.fa\n\n" + ENV_HELP)
    ap.add_argument("--profile", default=d["profile"])
    ap.add_argument("--params-file", default=d["params"], metavar="YAML")
    known, overrides = ap.parse_known_args(rest)
    p = _nextflow_params(known.profile, known.params_file, overrides)
    ok = True

    def row(status, label, detail=""):
        nonlocal ok
        if status == "MISSING":
            ok = False
        print(f"  {status:<8} {label:<34} {detail}")

    print(f"{PROG} check  (profile {known.profile}{', params ' + known.params_file if known.params_file else ''})\n\ntools:")
    for t, note in (("nextflow", "Nextflow >= 24.10"), ("java", "Java 17+"), ("sbatch", "SLURM"),
                    ("apptainer", "container engine"), ("singularity", "container engine (alias on Apptainer sites)")):
        path = shutil.which(t)
        row("OK" if path else ("MISSING" if t in ("nextflow", "java") else "absent"), t, path or note)
    if not (shutil.which("apptainer") or shutil.which("singularity")):
        row("MISSING", "container engine", "apptainer or singularity is required")

    print("\ncontainers:")
    for k in ("callforge_sif", "annotate_sif", "cohortqc_sif", "burden_sif"):
        v = p.get(k)
        row("OK" if v and os.path.isfile(v) else "MISSING", k, v or "(unset)")

    print("\nreference:")
    fa = p.get("fasta")
    if fa and os.path.isfile(fa):
        row("OK", "fasta", fa)
        for ext in (".fai", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
            row("OK" if os.path.isfile(fa + ext) else "MISSING", f"  fasta{ext}", "")
        dct = fa.rsplit(".", 1)[0] + ".dict"
        row("OK" if os.path.isfile(dct) else "MISSING", "  .dict", dct)
        style = None
        if os.path.isfile(fa + ".fai"):
            with open(fa + ".fai") as fh:
                first = fh.readline().split("\t")[0]
            style = "chr" if first.startswith("chr") else "nochr"
            n = sum(1 for _ in open(fa + ".fai"))
            alts = sum(1 for l in open(fa + ".fai") if any(t in l.split("\t")[0] for t in ("_alt", "_decoy", "EBV", "HLA-")))
            row("OK" if alts == 0 else "WARN", "  contigs", f"{n} contigs, {alts} alt/decoy, {style}-style names")
    else:
        row("MISSING", "fasta", fa or "(unset)"); style = None
    cr = p.get("calling_regions")
    if cr and os.path.isfile(cr):
        row("OK", "calling_regions", cr)
    else:
        row("WARN", "calling_regions", (cr + "  (not found: whole chromosomes will be used; run `genomeforge resources fetch`)") if cr else "(unset: whole chromosomes will be used)")

    print("\nannotation / known sites:")
    for k in ("dbsnp", "mills", "known_indels", "hapmap", "omni", "g1000", "clinvar", "somalier_sites"):
        v = p.get(k)
        if not v or not os.path.isfile(v):
            row("WARN" if k in ("hapmap", "omni", "g1000") else "MISSING", k,
                (v or "(unset)") + ("  -> VQSR falls back to hard filters; run `genomeforge resources fetch`" if k in ("hapmap", "omni", "g1000") else ""))
            continue
        idx = os.path.isfile(v + ".tbi") or os.path.isfile(v + ".csi")
        needs_idx = k != "somalier_sites"                      # somalier reads its sites file sequentially
        detail = v + ("" if idx or not needs_idx else "  (index .tbi MISSING)")
        st = "OK" if idx or not needs_idx else "MISSING"
        if style and k != "clinvar" and _have("tabix") and idx:
            try:
                first = subprocess.run(["tabix", "-l", v], capture_output=True, text=True, timeout=60).stdout.split("\n")[0]
                vs = "chr" if first.startswith("chr") else "nochr"
                if first and vs != style:
                    st, detail = "MISSING", f"{v}  (contig style {vs} != reference {style})"
            except (subprocess.SubprocessError, OSError):
                pass
        row(st, k, detail)
    gd, gp = p.get("gnomad_dir"), p.get("gnomad_pattern") or "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
    if gd and os.path.isdir(gd):
        have = [c for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"] if os.path.isfile(os.path.join(gd, gp.replace("{chrom}", c)))]
        row("OK" if len(have) >= 24 else ("WARN" if have else "MISSING"), "gnomad", f"{len(have)}/24 chromosome files in {gd}")
    else:
        row("MISSING", "gnomad_dir", gd or "(unset)")
    vc, vr = p.get("vep_cache"), p.get("vep_release") or "112"
    cache = os.path.join(vc or "", "homo_sapiens", f"{vr}_GRCh38")
    row("OK" if vc and os.path.isdir(cache) else ("WARN" if (p.get("run_vep") or "true") == "false" else "MISSING"), "vep_cache", cache)

    print("\nscheduler:")
    row("OK" if p.get("slurm_partition") else "INFO", "slurm_partition", p.get("slurm_partition") or "(cluster default)")
    print("\n" + ("READY: no missing items." if ok else "NOT READY: fix the MISSING items above."))
    return 0 if ok else 1


# ----------------------------------------------------------------------------- resources
def cmd_resources(rest):
    ap = argparse.ArgumentParser(
        prog=f"{PROG} resources", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Stage public hg38 resources that a typical annotation bundle lacks: the VQSR training "
                     "sets (hapmap 3.3, 1000G omni 2.5, 1000G phase-1 high-confidence SNPs) and the GATK "
                     "wgs_calling_regions interval list. Idempotent and resumable; ~2.5 GB."),
        epilog=f"examples:\n  {PROG} resources fetch /hpc/opt/callforge/refs/annotation_db\n  {PROG} resources fetch /refs/GRCh38\n")
    ap.add_argument("action", choices=["fetch"], help="fetch: download into DEST")
    ap.add_argument("dest", metavar="DEST", help="destination directory (the one your site profile points at)")
    a = ap.parse_args(rest)
    script = os.path.join(BIN, "fetch_resources.sh")
    return subprocess.call(["bash", script, a.dest])


# ----------------------------------------------------------------------------- status
def cmd_status(rest):
    d = _site_defaults()
    ap = argparse.ArgumentParser(
        prog=f"{PROG} status", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Show the state of a run directory: head job, task counts per process, recent log lines, results present.",
        epilog=f"examples:\n  {PROG} status\n  {PROG} status --rundir /hpc/data/<lab>/<project>/genomeforge_run --lines 40\n")
    ap.add_argument("--rundir", default=d["rundir"])
    ap.add_argument("--lines", type=int, default=15, help="log lines to show (default 15)")
    a = ap.parse_args(rest)
    rundir = os.path.abspath(a.rundir)
    if not os.path.isdir(rundir):
        _die(f"run directory not found: {rundir}")
    print(f"run dir: {rundir}")
    if _have("squeue"):
        q = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%i %j %T %M %R"], capture_output=True, text=True).stdout
        mine = [l for l in q.splitlines() if "genomeforge" in l or "nf-" in l]
        print(f"SLURM jobs ({len(mine)} nextflow/genomeforge jobs in your queue):")
        for l in mine[:12]:
            print("  " + l)
        if len(mine) > 12:
            print(f"  ... {len(mine) - 12} more")
    heads = sorted(glob.glob(os.path.join(rundir, "logs", "head_*.out")))
    if heads:
        last = heads[-1]; stamp = os.path.basename(last)[5:-4]
        state = "COMPLETE" if os.path.isfile(os.path.join(rundir, "logs", f"head_{stamp}.complete")) else "in progress / stopped"
        attempts = sum(1 for l in open(last, errors="replace") if l.startswith("=== GenomeForge head job"))
        print(f"head chain: {os.path.basename(last)}  attempts={attempts}  state={state}")
        with open(last, errors="replace") as fh:
            tail = fh.readlines()[-6:]
        for l in tail:
            print("  " + l.rstrip())
    log = os.path.join(rundir, "nextflow.log")
    trace = None
    for root in (os.path.join(rundir, "results", "pipeline_info"),):
        if os.path.isfile(os.path.join(root, "trace.txt")):
            trace = os.path.join(root, "trace.txt")
    if trace:
        counts = {}
        with open(trace) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ip, istat = hdr.index("process"), hdr.index("status")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(ip, istat):
                    continue
                proc = f[ip].split(" ")[0]
                counts.setdefault(proc, {}).setdefault(f[istat], 0)
                counts[proc][f[istat]] += 1
        print("tasks per process (from trace.txt):")
        for proc in sorted(counts):
            print(f"  {proc:<22} " + "  ".join(f"{s}={n}" for s, n in sorted(counts[proc].items())))
    if os.path.isfile(log):
        with open(log, errors="replace") as fh:
            lines = fh.readlines()
        print(f"last {a.lines} lines of nextflow.log:")
        for l in lines[-a.lines:]:
            print("  " + l.rstrip())
    else:
        print("no nextflow.log yet")
    res = os.path.join(rundir, "results")
    if os.path.isdir(res):
        print("results:", ", ".join(sorted(os.listdir(res))))
    return 0


# ----------------------------------------------------------------------------- qc-summary
def cmd_qc_summary(rest):
    ap = argparse.ArgumentParser(
        prog=f"{PROG} qc-summary", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Rebuild cohort_qc_summary.tsv/.md from a results directory, e.g. with different "
                     "thresholds or after editing the sample sheet's sex column."),
        epilog=f"examples:\n  {PROG} qc-summary --results genomeforge_run/results --samplesheet samplesheet.csv\n"
               f"  {PROG} qc-summary --results results --samplesheet samplesheet.csv --min-depth 25 --out-dir qc_strict\n")
    ap.add_argument("--results", required=True, metavar="DIR")
    ap.add_argument("--samplesheet", required=True, metavar="CSV")
    ap.add_argument("--out-dir", default=".", metavar="DIR")
    ap.add_argument("--min-depth", default="20"); ap.add_argument("--max-dup", default="0.25")
    ap.add_argument("--max-contam", default="0.03"); ap.add_argument("--min-mapped", default="0.95")
    a = ap.parse_args(rest)
    r = os.path.abspath(a.results)
    qc = os.path.join(r, "03_sample_qc"); co = os.path.join(r, "07_cohort_qc")
    os.makedirs(a.out_dir, exist_ok=True)
    args = ["--flagstat", os.path.join(qc, "samtools"), "--markdup", os.path.join(qc, "markdup"),
            "--mosdepth", os.path.join(qc, "mosdepth"), "--wgs-metrics", os.path.join(qc, "wgs_metrics"),
            "--verifybamid", os.path.join(qc, "verifybamid"),
            "--somalier-samples", os.path.join(co, "somalier", "cohort.samples.tsv"),
            "--somalier-pairs", os.path.join(co, "somalier", "cohort.pairs.tsv"),
            "--bcftools-stats", os.path.join(co, "cohort.bcftools_stats.txt"),
            "--samplesheet", a.samplesheet, "--min-depth", a.min_depth, "--max-dup", a.max_dup,
            "--max-contam", a.max_contam, "--min-mapped", a.min_mapped,
            "--out-tsv", os.path.join(a.out_dir, "cohort_qc_summary.tsv"),
            "--out-md", os.path.join(a.out_dir, "cohort_qc_summary.md")]
    return _run_bin("summarize_qc.py", f"{PROG} qc-summary", args)


# ----------------------------------------------------------------------------- test
def cmd_test(rest):
    ap = argparse.ArgumentParser(
        prog=f"{PROG} test", formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Installation check: builds a synthetic FASTQ/reference layout in a temporary directory "
                     "(Novogene, bcl2fastq and generic file names), builds the sample sheet, and dry-runs the "
                     "complete DAG with stub tasks on the local executor. Needs only Nextflow and Java."),
        epilog=f"examples:\n  {PROG} test\n  {PROG} test --keep --dir /tmp/gf_test\n")
    ap.add_argument("--dir", help="directory for the synthetic layout (default: a temporary one)")
    ap.add_argument("--keep", action="store_true", help="keep the directory afterwards")
    a = ap.parse_args(rest)
    if not _have("nextflow"):
        _die("`nextflow` not found on PATH", 127)
    tdir = a.dir or tempfile.mkdtemp(prefix="genomeforge_test_")
    harness = os.path.join(REPO_ROOT, "test", "make_stub_env.sh")
    rc = subprocess.call(["bash", harness, tdir], stdout=subprocess.DEVNULL)
    if rc != 0:
        _die("could not build the synthetic layout")
    T = tdir
    cmd = ["nextflow", "run", os.path.join(REPO_ROOT, "main.nf"), "-profile", "stub", "-stub-run",
           "--input", "samplesheet.csv", "--outdir", "results", "--callforge_sif", "x", "--cohortqc_sif", "x", "--burden_sif", "x",
           "--fasta", f"{T}/ref/ref.fa", "--dbsnp", f"{T}/db/dbsnp.chr.vcf.gz", "--mills", f"{T}/db/mills.vcf.gz",
           "--known_indels", f"{T}/db/known_indels.vcf.gz", "--hapmap", f"{T}/db/hapmap.vcf.gz", "--omni", f"{T}/db/omni.vcf.gz",
           "--g1000", f"{T}/db/g1000.vcf.gz", "--clinvar", f"{T}/db/clinvar.vcf.gz", "--gnomad_dir", f"{T}/db",
           "--somalier_sites", f"{T}/db/sites.hg38.vcf.gz", "--vep_cache", f"{T}/vep", "--chromosomes", "chr1,chr22,chrM"]
    _say(f"synthetic layout in {T}; running the stub DAG ...")
    rc = subprocess.call(cmd, cwd=T, env={**os.environ, "NXF_ANSI_LOG": "false"})
    if rc == 0:
        print(f"\nOK: GenomeForge {__version__} compiles and the full DAG completes in stub mode.")
    else:
        print(f"\nFAILED (exit {rc}); see {T}/.nextflow.log")
    if not a.keep and not a.dir:
        shutil.rmtree(T, ignore_errors=True)
    return rc


# ----------------------------------------------------------------------------- dispatch
SUBCOMMANDS = {
    "run":         (cmd_run,        "launch or resume the pipeline (sbatch head job; resumable)"),
    "samplesheet": (cmd_samplesheet, "build the lane-level sample sheet from FASTQ directories or a pairs CSV"),
    "check":       (cmd_check,      "verify references, containers and tools for a site profile"),
    "resources":   (cmd_resources,  "fetch the VQSR training sets and GATK calling regions"),
    "status":      (cmd_status,     "show the state of a run directory"),
    "qc-summary":  (cmd_qc_summary, "rebuild the per-sample QC table from a results directory"),
    "test":        (cmd_test,       "dry-run the whole DAG on a synthetic layout (installation check)"),
}
ALIASES = {"make-samplesheet": "samplesheet", "sheet": "samplesheet", "fetch": "resources"}


def _top_help():
    lines = [f"genomeforge {__version__} — whole-genome germline joint calling (FASTQ → annotated, QC'd cohort callset)", "",
             "usage: genomeforge <command> [options]", "", "commands:"]
    for name, (_, desc) in SUBCOMMANDS.items():
        lines.append(f"  {name:<13} {desc}")
    lines += ["", "  --version      print version", "  --help         show this help", "",
              "typical session:",
              "  genomeforge check --profile noguchi                                    # site ready?",
              "  genomeforge resources fetch /hpc/opt/callforge/refs/annotation_db       # once per site",
              "  genomeforge samplesheet --roots /data/delivery --metadata phenotypes.csv --out samplesheet.csv",
              "  genomeforge run --profile noguchi --input samplesheet.csv --rundir /hpc/data/<lab>/<project>/genomeforge_run --stub",
              "  genomeforge run --profile noguchi --input samplesheet.csv --rundir /hpc/data/<lab>/<project>/genomeforge_run",
              "  genomeforge status --rundir /hpc/data/<lab>/<project>/genomeforge_run",
              "", "Run `genomeforge <command> --help` for every option and more examples.", "", ENV_HELP.rstrip()]
    return "\n".join(lines)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(_top_help())
        return 0
    if argv[0] in ("-V", "--version", "version"):
        print(f"genomeforge {__version__}")
        return 0
    cmd, rest = argv[0], argv[1:]
    cmd = ALIASES.get(cmd, cmd)
    if cmd not in SUBCOMMANDS:
        sys.stderr.write(f"[{PROG}] unknown command '{cmd}'\n\n{_top_help()}\n")
        return 2
    return SUBCOMMANDS[cmd][0](rest)


if __name__ == "__main__":
    sys.exit(main())
