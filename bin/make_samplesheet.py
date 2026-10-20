#!/usr/bin/env python3
"""make_samplesheet.py — build the GenomeForge sample sheet (one row per lane-pair) from
FASTQ directories of any Illumina naming convention. READ-ONLY on the data.

  make_samplesheet.py --roots DIR [DIR ...] --out samplesheet.csv \\
      [--metadata clinical.csv] [--exclude S1,S2] [--min-bytes 1000000000] [--sample-from auto|dir|filename]

  make_samplesheet.py --fastq-csv pairs.csv --out samplesheet.csv       # from an existing sample_id,fastq_1,fastq_2 list

Recognised file names (R1/R2 mates must differ only in the read token):
  Novogene      <sample>_<library>_<flowcell>_L<lane>_<1|2>.fq.gz
  bcl2fastq     <sample>_S<n>_L00<lane>_R<1|2>_001.fastq.gz
  generic       <sample>[_<flowcell>][_L<lane>][_.]R?<1|2>[_001].(fastq|fq)[.gz]
Flowcell and lane are taken from the file name when present, otherwise from the first read
header (Illumina `@instrument:run:flowcell:lane:...`). A library token is used when the name
carries one, else the sample id.

Columns written: sample_id,lane_id,library,flowcell,lane,fastq_1,fastq_2,sex,phenotype
  lane_id = <flowcell>.L<lane>  (read-group ID; must be unique within a sample)
Metadata CSV (optional) needs a sample_id column; sex and phenotype are copied if present.

Fails loudly on unpaired FASTQs, duplicate lane IDs, unreadable files, or files smaller than
--min-bytes (a few hundred MB "lane" is an interrupted transfer) unless --allow-small.
"""
import argparse, csv, gzip, os, re, sys
from collections import defaultdict

PATTERNS = [
    ("novogene",  re.compile(r"^(?P<sample>.+?)_(?P<lib>[A-Za-z0-9]+-[A-Za-z0-9]+)_(?P<fc>[A-Za-z0-9]{5,})_L(?P<lane>\d+)_(?P<read>[12])\.(fastq|fq)(\.gz)?$")),
    ("bcl2fastq", re.compile(r"^(?P<sample>.+?)_S(?P<snum>\d+)(?:_L(?P<lane>\d{3}))?_R(?P<read>[12])(?:_\d{3})?\.(fastq|fq)(\.gz)?$")),
    ("generic",   re.compile(r"^(?P<sample>.+?)(?:_(?P<fc>[A-Z0-9]{5,}))?(?:_L(?P<lane>\d+))?[._]R?(?P<read>[12])(?:_\d{3})?\.(fastq|fq)(\.gz)?$")),
]
HEADER = re.compile(r"^@(?P<inst>[^:\s]+):(?P<run>\d+):(?P<fc>[^:\s]+):(?P<lane>\d+):")
SKIP_DIRS = {"02.bam", "work", "results", "undetermined"}


def first_header(path):
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt", errors="replace") as fh:
            line = fh.readline()
        m = HEADER.match(line)
        return (m.group("fc"), m.group("lane")) if m else (None, None)
    except Exception:
        return (None, None)


def parse_name(fn):
    for scheme, rx in PATTERNS:
        m = rx.match(fn)
        if m:
            d = m.groupdict()
            return scheme, d["sample"], d.get("lib"), d.get("fc"), (d.get("lane") or "").lstrip("0") or d.get("lane"), d["read"]
    return None


def scan(roots, sample_from):
    files = []
    for root in roots:
        for dp, dns, fns in os.walk(root):
            dns[:] = [d for d in dns if d.lower() not in SKIP_DIRS and not d.lower().startswith("03.result")]
            for fn in sorted(fns):
                if not re.search(r"\.(fastq|fq)(\.gz)?$", fn):
                    continue
                parsed = parse_name(fn)
                if not parsed:
                    sys.stderr.write(f"[samplesheet] WARNING: unrecognised FASTQ name skipped: {fn}\n"); continue
                scheme, sample, lib, fc, lane, read = parsed
                parent = os.path.basename(dp)
                if sample_from == "dir" or (sample_from == "auto" and fn.startswith(parent) and parent.lower() not in ("fastq", "fastqs", "raw", "rawdata", "01.rawdata")):
                    sample = parent
                path = os.path.abspath(os.path.join(dp, fn))
                if not fc or not lane:
                    hfc, hlane = first_header(path)
                    fc = fc or hfc or "NA"; lane = lane or hlane or "1"
                files.append({"sample": sample, "lib": lib or sample, "fc": fc, "lane": lane, "read": read,
                              "path": path, "size": os.path.getsize(path), "scheme": scheme})
    return files


def read_pairs_csv(path):
    files = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            sid = r["sample_id"].strip()
            for read, col in (("1", "fastq_1"), ("2", "fastq_2")):
                p = os.path.abspath(r[col].strip())
                fc, lane = (r.get("flowcell") or "").strip(), (r.get("lane") or "").strip()
                if not fc or not lane:
                    hfc, hlane = first_header(p); fc = fc or hfc or "NA"; lane = lane or hlane or "1"
                files.append({"sample": sid, "lib": (r.get("library") or sid).strip(), "fc": fc, "lane": lane, "read": read,
                              "path": p, "size": os.path.getsize(p) if os.path.exists(p) else 0, "scheme": "csv"})
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--roots", nargs="+", help="directories to scan recursively")
    src.add_argument("--fastq-csv", help="existing CSV with sample_id,fastq_1,fastq_2[,flowcell,lane,library]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--metadata", help="CSV with sample_id[,sex,phenotype,...]")
    ap.add_argument("--exclude", default="", help="comma list of sample_ids to leave out")
    ap.add_argument("--sample-from", choices=["auto", "dir", "filename"], default="auto",
                    help="sample id from the parent directory (Novogene layout) or the file name; auto = directory when the name starts with it")
    ap.add_argument("--min-bytes", type=int, default=1_000_000_000, help="smallest acceptable lane FASTQ (default 1 GB)")
    ap.add_argument("--allow-small", action="store_true", help="keep files below --min-bytes instead of failing")
    ap.add_argument("--drop-unpaired", action="store_true", help="skip unpaired files instead of failing")
    a = ap.parse_args()

    excl = {s.strip() for s in a.exclude.split(",") if s.strip()}
    files = read_pairs_csv(a.fastq_csv) if a.fastq_csv else scan(a.roots, a.sample_from)
    files = [f for f in files if f["sample"] not in excl]
    units = defaultdict(dict)
    for f in files:
        key = (f["sample"], f["fc"], f["lane"], f["lib"])
        if f["read"] in units[key]:
            sys.exit(f"[samplesheet] ERROR: two R{f['read']} files for sample={key[0]} flowcell={key[1]} lane={key[2]}:\n  {units[key][f['read']]['path']}\n  {f['path']}")
        units[key][f["read"]] = f

    meta = {}
    if a.metadata:
        with open(a.metadata, newline="") as fh:
            for r in csv.DictReader(fh):
                sid = (r.get("sample_id") or "").strip()
                if sid: meta[sid] = r

    rows, problems = [], []
    for (sample, fc, lane, lib), d in sorted(units.items()):
        if "1" not in d or "2" not in d:
            msg = f"unpaired: {sample} {fc} L{lane} has only R{list(d)[0]} ({list(d.values())[0]['path']})"
            if a.drop_unpaired: sys.stderr.write(f"[samplesheet] WARNING {msg} (dropped)\n"); continue
            problems.append(msg); continue
        for r in ("1", "2"):
            f = d[r]
            if not os.access(f["path"], os.R_OK): problems.append(f"unreadable: {f['path']}")
            if f["size"] < a.min_bytes:
                msg = f"too small ({f['size']/1e9:.2f} GB): {f['path']}"
                if a.allow_small: sys.stderr.write(f"[samplesheet] WARNING {msg}\n")
                else: problems.append(msg)
        m = meta.get(sample, {})
        rows.append({"sample_id": sample, "lane_id": f"{fc}.L{lane}", "library": lib, "flowcell": fc, "lane": lane,
                     "fastq_1": d["1"]["path"], "fastq_2": d["2"]["path"],
                     "sex": (m.get("sex") or "").strip(), "phenotype": (m.get("phenotype") or "").strip()})

    seen = defaultdict(set)
    for r in rows:
        if r["lane_id"] in seen[r["sample_id"]]:
            problems.append(f"duplicate lane_id {r['lane_id']} for {r['sample_id']} (two libraries on one lane: exclude one or split the sample)")
        seen[r["sample_id"]].add(r["lane_id"])
    bad_ids = sorted({r["sample_id"] for r in rows if not re.fullmatch(r"[A-Za-z0-9_.-]+", r["sample_id"])})
    if bad_ids: problems.append(f"sample ids with characters outside [A-Za-z0-9_.-]: {', '.join(bad_ids)}")

    if problems:
        sys.stderr.write("[samplesheet] ERROR — sheet NOT written:\n" + "".join(f"  - {p}\n" for p in problems))
        sys.exit(1)
    if not rows:
        sys.exit("[samplesheet] ERROR: no lane-pairs found")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample_id", "lane_id", "library", "flowcell", "lane", "fastq_1", "fastq_2", "sex", "phenotype"])
        w.writeheader(); w.writerows(rows)
    samples = {r["sample_id"] for r in rows}
    multi = sum(1 for s in samples if len(seen[s]) > 1)
    schemes = sorted({f["scheme"] for f in files})
    print(f"[samplesheet] wrote {a.out}: {len(rows)} lane-pairs, {len(samples)} samples ({multi} multi-lane), naming: {','.join(schemes)}, "
          f"{sum(1 for s in samples if s in meta)} with metadata, "
          f"{sum(os.path.getsize(r['fastq_1']) + os.path.getsize(r['fastq_2']) for r in rows)/1e12:.2f} TB")
    if a.metadata:
        missing = sorted(samples - set(meta))
        if missing: print(f"[samplesheet] samples without metadata: {', '.join(missing)}")
        extra = sorted(set(meta) - samples)
        if extra: print(f"[samplesheet] metadata rows without sequence data: {', '.join(extra)}")


if __name__ == "__main__":
    main()
