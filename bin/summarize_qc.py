#!/usr/bin/env python3
"""summarize_qc.py — one row per sample from flagstat, MarkDuplicates, mosdepth,
CollectWgsMetrics, verifyBamID2, somalier and bcftools stats; flags (never drops) samples
against WGS thresholds and checks sex against the sample sheet."""
import argparse, csv, glob, os, re, statistics, sys
from collections import defaultdict

AUTOS = {f"chr{i}" for i in range(1, 23)}


def read_flagstat(p):
    d = {}
    for line in open(p):
        m = re.match(r"(\d+) \+ \d+ (.*)", line.strip())
        if not m: continue
        n, k = int(m.group(1)), m.group(2)
        if k.startswith("in total"): d["total_reads"] = n
        elif k.startswith("mapped (") and "primary" not in k: d["mapped"] = n
        elif k.startswith("properly paired"): d["properly_paired"] = n
        elif k.startswith("duplicates") or k.startswith("primary duplicates"): d.setdefault("flag_dups", n)
    if d.get("total_reads"):
        d["pct_mapped"] = round(d.get("mapped", 0) / d["total_reads"], 4)
        d["pct_properly_paired"] = round(d.get("properly_paired", 0) / d["total_reads"], 4)
    return d


def read_picard_metrics(p):
    """Return list of dict rows of the first METRICS CLASS table."""
    rows, header, in_tab = [], None, False
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith("## METRICS CLASS"): in_tab = True; continue
        if in_tab:
            if not line.strip(): break
            if header is None: header = line.split("\t"); continue
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def read_markdup(p):
    rows = read_picard_metrics(p)
    ex = sum(float(r.get("READ_PAIRS_EXAMINED") or 0) for r in rows)
    du = sum(float(r.get("READ_PAIR_DUPLICATES") or 0) for r in rows)
    od = sum(float(r.get("READ_PAIR_OPTICAL_DUPLICATES") or 0) for r in rows)
    return {"dup_rate": round(du / ex, 4) if ex else None, "optical_dup_rate": round(od / ex, 4) if ex else None,
            "n_libraries": len(rows)}


def read_mosdepth(p):
    means, lens = {}, {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "chrom" or len(f) < 4: continue
        try: means[f[0]] = float(f[3]); lens[f[0]] = int(f[1])
        except ValueError: pass
    auto = [(means[c], lens[c]) for c in AUTOS if c in means]
    auto_mean = sum(m * l for m, l in auto) / sum(l for _, l in auto) if auto else None
    d = {"mosdepth_mean": means.get("total"), "auto_mean": round(auto_mean, 2) if auto_mean else None,
         "chrX_mean": means.get("chrX"), "chrY_mean": means.get("chrY")}
    if auto_mean:
        xr = (means.get("chrX") or 0) / auto_mean; yr = (means.get("chrY") or 0) / auto_mean
        d["X_ratio"] = round(xr, 3); d["Y_ratio"] = round(yr, 3)
        d["sex_coverage"] = "M" if yr > 0.15 and xr < 0.8 else ("F" if yr < 0.05 and xr > 0.8 else "ambiguous")
    return d


def read_wgs_metrics(p):
    rows = read_picard_metrics(p)
    if not rows: return {}
    r = rows[0]
    return {k.lower(): r.get(k) for k in ("MEAN_COVERAGE", "MEDIAN_COVERAGE", "SD_COVERAGE", "PCT_10X", "PCT_20X", "PCT_30X", "PCT_EXC_DUPE", "PCT_EXC_MAPQ", "HET_SNP_SENSITIVITY")}


def read_selfsm(p):
    with open(p) as fh:
        r = list(csv.DictReader(fh, delimiter="\t"))
    if not r: return {}
    v = r[0].get("FREEMIX", "NA")
    try: return {"freemix": round(float(v), 4)}
    except ValueError: return {"freemix": None}


def read_somalier_samples(p):
    out = {}
    if not p or not os.path.isfile(p): return out
    with open(p) as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"): hdr = [h.lstrip("#") for h in f]; continue
            if hdr is None or len(f) < len(hdr): continue
            r = dict(zip(hdr, f))
            # Without a pedigree somalier leaves `sex` at -9 even with --infer; infer it here from the
            # X heterozygosity ratio and the Y depth, the same signals somalier's HTML uses.
            sex = {"1": "M", "2": "F"}.get(r.get("sex", ""))
            if sex is None:
                try:
                    xhet, xalt, xn = float(r.get("X_het") or 0), float(r.get("X_hom_alt") or 0), float(r.get("X_n") or 0)
                    depth, ydepth = float(r.get("gt_depth_mean") or r.get("depth_mean") or 0), float(r.get("Y_depth_mean") or 0)
                except ValueError:
                    xhet = xalt = xn = depth = ydepth = 0.0
                ratio = xhet / (xhet + xalt) if (xhet + xalt) > 0 else None
                yfrac = ydepth / depth if depth > 0 else None
                if xn < 50 or ratio is None: sex = "unknown"
                elif ratio >= 0.2 and (yfrac is None or yfrac < 0.15): sex = "F"
                elif ratio < 0.1 and (yfrac is None or yfrac >= 0.15): sex = "M"
                else: sex = "ambiguous"
            out[r["sample_id"]] = {"somalier_sex": sex, "somalier_X_het": r.get("X_het"), "somalier_Y_depth": r.get("Y_depth_mean")}
    return out


def read_somalier_pairs(p, thr=0.2):
    rel = defaultdict(list)
    if not p or not os.path.isfile(p): return rel
    with open(p) as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"): hdr = [h.lstrip("#") for h in f]; continue
            if hdr is None: continue
            r = dict(zip(hdr, f))
            try: v = float(r.get("relatedness", "nan"))
            except ValueError: continue
            if v >= thr:
                rel[r["sample_a"]].append(f"{r['sample_b']}({v:.2f})"); rel[r["sample_b"]].append(f"{r['sample_a']}({v:.2f})")
    return rel


def read_bcftools_psc(p):
    out = {}
    if not p or not os.path.isfile(p): return out
    for line in open(p):
        if not line.startswith("PSC"): continue
        f = line.rstrip("\n").split("\t")
        # 0:PSC 1:id 2:sample 3:nRefHom 4:nNonRefHom 5:nHets 6:nTransitions 7:nTransversions 8:nIndels 9:avgDP 10:nSingletons
        try:
            hom, het, ti, tv, ind, sing = int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[8]), int(f[10])
        except (ValueError, IndexError): continue
        out[f[2]] = {"n_snv_nonref": hom + het, "het_hom_ratio": round(het / hom, 3) if hom else None,
                     "titv": round(ti / tv, 3) if tv else None, "n_indels": ind, "n_singletons": sing}
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("flagstat", "markdup", "mosdepth", "wgs-metrics", "verifybamid"): ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--somalier-samples"); ap.add_argument("--somalier-pairs"); ap.add_argument("--bcftools-stats")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--min-depth", type=float, default=20); ap.add_argument("--max-dup", type=float, default=0.25)
    ap.add_argument("--max-contam", type=float, default=0.03); ap.add_argument("--min-mapped", type=float, default=0.95)
    ap.add_argument("--out-tsv", required=True); ap.add_argument("--out-md", required=True)
    a = ap.parse_args()

    sheet = {}
    with open(a.samplesheet, newline="") as fh:
        for r in csv.DictReader(fh):
            sheet.setdefault(r["sample_id"], {"sex": (r.get("sex") or "").strip().upper()[:1], "phenotype": (r.get("phenotype") or "").strip(), "lanes": 0})
            sheet[r["sample_id"]]["lanes"] += 1

    data = defaultdict(dict)
    def collect(dirn, pattern, parser, suffix):
        for p in glob.glob(os.path.join(dirn, pattern)):
            sid = os.path.basename(p)[:-len(suffix)]
            try: data[sid].update(parser(p))
            except Exception as e: data[sid]["parse_error"] = f"{suffix}: {e}"
    collect(a.flagstat, "*.flagstat", read_flagstat, ".flagstat")
    collect(a.markdup, "*.markdup.metrics.txt", read_markdup, ".markdup.metrics.txt")
    collect(a.mosdepth, "*.mosdepth.summary.txt", read_mosdepth, ".mosdepth.summary.txt")
    collect(a.wgs_metrics, "*.wgs_metrics.txt", read_wgs_metrics, ".wgs_metrics.txt")
    collect(a.verifybamid, "*.selfSM", read_selfsm, ".selfSM")
    som = read_somalier_samples(a.somalier_samples); rel = read_somalier_pairs(a.somalier_pairs); psc = read_bcftools_psc(a.bcftools_stats)
    for sid, d in som.items(): data[sid].update(d)
    for sid, d in psc.items(): data[sid].update(d)

    cols = ["sample_id", "phenotype", "n_lanes", "total_reads", "pct_mapped", "pct_properly_paired", "dup_rate", "optical_dup_rate",
            "mean_coverage", "median_coverage", "pct_10x", "pct_20x", "pct_30x", "mosdepth_mean", "X_ratio", "Y_ratio",
            "sex_sheet", "sex_coverage", "somalier_sex", "freemix", "n_snv_nonref", "het_hom_ratio", "titv", "n_indels", "n_singletons",
            "related_to", "qc_flags", "qc_status"]
    rows = []
    for sid in sorted(set(data) | set(sheet)):
        d = data.get(sid, {}); s = sheet.get(sid, {})
        flags = []
        depth = None
        for k in ("mean_coverage", "mosdepth_mean"):
            try: depth = float(d.get(k)); break
            except (TypeError, ValueError): pass
        if depth is not None and depth < a.min_depth: flags.append(f"low_depth({depth:.1f}x)")
        if d.get("dup_rate") is not None and d["dup_rate"] > a.max_dup: flags.append(f"high_dup({d['dup_rate']:.2f})")
        if d.get("freemix") is not None and d["freemix"] > a.max_contam: flags.append(f"contamination({d['freemix']:.3f})")
        if d.get("pct_mapped") is not None and d["pct_mapped"] < a.min_mapped: flags.append(f"low_mapping({d['pct_mapped']:.3f})")
        sx = s.get("sex", "")
        inferred = d.get("somalier_sex") if d.get("somalier_sex") in ("M", "F") else d.get("sex_coverage")
        if sx in ("M", "F") and inferred in ("M", "F") and sx != inferred: flags.append(f"sex_mismatch(sheet={sx},inferred={inferred})")
        if d.get("sex_coverage") == "ambiguous": flags.append("sex_ambiguous")
        if d.get("somalier_sex") and d.get("sex_coverage") in ("M", "F") and d["somalier_sex"] in ("M", "F") and d["somalier_sex"] != d["sex_coverage"]:
            flags.append("sex_inference_disagree")
        if rel.get(sid): flags.append("related:" + ";".join(rel[sid]))
        if d.get("titv") is not None and not (1.9 <= d["titv"] <= 2.2): flags.append(f"titv({d['titv']})")
        if not d: flags.append("no_qc_data")
        row = {"sample_id": sid, "phenotype": s.get("phenotype", ""), "n_lanes": s.get("lanes", ""), "sex_sheet": sx,
               "related_to": ";".join(rel.get(sid, [])), "qc_flags": ";".join(flags), "qc_status": "FLAG" if flags else "PASS"}
        for c in cols:
            if c not in row: row[c] = d.get(c, "")
        rows.append(row)

    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

    def med(key):
        v = []
        for r in rows:
            try: v.append(float(r[key]))
            except (TypeError, ValueError): pass
        return (f"{statistics.median(v):.3g}" if v else "n/a"), (f"{min(v):.3g}" if v else "n/a"), (f"{max(v):.3g}" if v else "n/a")
    flagged = [r for r in rows if r["qc_status"] != "PASS"]
    with open(a.out_md, "w") as fh:
        fh.write(f"# Cohort QC summary\n\n{len(rows)} samples, {len(flagged)} flagged (thresholds: depth >= {a.min_depth}x, dup <= {a.max_dup}, FREEMIX <= {a.max_contam}, mapped >= {a.min_mapped}).\n\n")
        fh.write("| metric | median | min | max |\n|---|---|---|---|\n")
        for k, lab in (("mean_coverage", "mean coverage (CollectWgsMetrics)"), ("pct_20x", "fraction >= 20x"), ("dup_rate", "duplicate rate"),
                       ("pct_mapped", "mapped fraction"), ("freemix", "FREEMIX"), ("titv", "Ti/Tv"), ("het_hom_ratio", "het/hom"), ("n_snv_nonref", "non-ref SNVs")):
            m, lo, hi = med(k); fh.write(f"| {lab} | {m} | {lo} | {hi} |\n")
        if flagged:
            fh.write("\n## Flagged samples\n\n| sample | flags |\n|---|---|\n")
            for r in flagged: fh.write(f"| {r['sample_id']} | {r['qc_flags']} |\n")
        sexes = defaultdict(int)
        for r in rows: sexes[r.get("somalier_sex") if r.get("somalier_sex") in ("M", "F") else (r.get("sex_coverage") or "unknown")] += 1
        fh.write(f"\nInferred sex counts: {dict(sexes)}\n")
    print(f"[summarize_qc] {len(rows)} samples, {len(flagged)} flagged -> {a.out_tsv}")


if __name__ == "__main__":
    main()
