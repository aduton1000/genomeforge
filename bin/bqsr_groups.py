#!/usr/bin/env python3
"""Split the contigs of a sequence dictionary into N contiguous, length-balanced groups.

Writes bqsr_groups/<0000..>.intervals (one contig name per line, GATK -L format); the
last file is named <NNNN>.last.intervals so the pipeline can add the unmapped reads to
it. Whole contigs only, in dictionary order: no read is ever split across shards and the
gathered BAM is in reference order.
"""
import argparse
import os

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("dict", help="reference .dict file")
ap.add_argument("-n", "--groups", type=int, default=24, help="number of groups (default 24)")
ap.add_argument("-o", "--outdir", default="bqsr_groups")
a = ap.parse_args()

ctgs = []
for line in open(a.dict):
    if line.startswith("@SQ"):
        f = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:])
        ctgs.append((f["SN"], int(f["LN"])))
if not ctgs:
    raise SystemExit(f"no @SQ lines in {a.dict}")

n = max(1, a.groups)
target = sum(l for _, l in ctgs) / n
groups, cur, size = [], [], 0
for name, ln in ctgs:
    if cur and size + ln > target and len(groups) < n - 1:
        groups.append(cur)
        cur, size = [], 0
    cur.append(name)
    size += ln
if cur:
    groups.append(cur)

os.makedirs(a.outdir, exist_ok=True)
for i, g in enumerate(groups):
    suffix = ".last" if i == len(groups) - 1 else ""
    with open(os.path.join(a.outdir, "%04d%s.intervals" % (i, suffix)), "w") as fh:
        fh.write("\n".join(g) + "\n")
print(f"bqsr groups: {len(groups)} (from {len(ctgs)} contigs, target {target/1e6:.0f} Mb each)")
