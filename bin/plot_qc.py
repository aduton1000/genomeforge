#!/usr/bin/env python3
"""plot_qc.py — cohort QC figures + a self-contained HTML dashboard for a GenomeForge run.

Every figure is optional: a missing or unreadable input skips that figure with a note, so
the stage never fails the pipeline. Figures are written as PNG (and SVG) under --out-dir and
embedded into --dashboard.
"""
import argparse, base64, csv, glob, gzip, io, json, math, os, re, sys, traceback
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

PAL = ["#2FA38E", "#3B6EA8", "#B8862B", "#9E2B25", "#6F4E9C", "#5A8F3A", "#C75D8A", "#4F7C8A"]
AUTOS = [f"chr{i}" for i in range(1, 23)]
CHROMS = AUTOS + ["chrX", "chrY", "chrM"]
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.grid": True, "grid.alpha": 0.25})
FIGS = []        # (section, title, filename, caption)
NOTES = []


def save(fig, name, section, title, caption=""):
    os.makedirs(OUT, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name + ".png"))
    try: fig.savefig(os.path.join(OUT, name + ".svg"))
    except Exception: pass
    plt.close(fig)
    FIGS.append((section, title, name + ".png", caption))


def figure(section, title):
    """Decorator: run a plotting function, never let it break the run."""
    def deco(fn):
        def wrapped(*a, **k):
            try:
                r = fn(*a, **k)
                if r is False: NOTES.append(f"{title}: skipped (no input)")
            except Exception as e:
                NOTES.append(f"{title}: skipped ({type(e).__name__}: {e})")
                plt.close("all")
        wrapped.section, wrapped.title = section, title
        return wrapped
    return deco


def read_picard(p):
    rows, header, on = [], None, False
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith("## METRICS CLASS"): on = True; continue
        if on:
            if not line.strip(): break
            if header is None: header = line.split("\t"); continue
            rows.append(dict(zip(header, line.split("\t"))))
    return rows


def sample_of(path, suffix):
    return os.path.basename(path)[:-len(suffix)]


def stats_section(path, tag):
    out = []
    for line in open(path):
        if line.startswith(tag + "\t"):
            out.append(line.rstrip("\n").split("\t")[1:])
    return out


def status_colors(qc):
    return [PAL[3] if qc.get(s, "PASS") != "PASS" else PAL[0] for s in qc.get("_order", [])]


# ----------------------------------------------------------------------------- figures
@figure("Alignment", "Mean coverage per sample")
def fig_coverage(qc):
    d = qc[qc["mean_coverage"].notna()].sort_values("mean_coverage")
    if d.empty: return False
    fig, ax = plt.subplots(figsize=(max(7, 0.16 * len(d)), 4.2))
    cols = [PAL[3] if s != "PASS" else PAL[0] for s in d["qc_status"]]
    ax.bar(d["sample_id"], d["mean_coverage"], color=cols)
    ax.axhline(ARGS.min_depth, color=PAL[2], ls="--", lw=1, label=f"threshold {ARGS.min_depth}x")
    ax.set_ylabel("mean coverage (x)"); ax.set_title("Mean genome coverage (CollectWgsMetrics)"); ax.legend(loc="upper left")
    ax.tick_params(axis="x", rotation=90, labelsize=6)
    save(fig, "coverage_per_sample", "Alignment", "Mean coverage per sample", "Red bars are samples with at least one QC flag.")


@figure("Alignment", "Cumulative coverage")
def fig_cum_coverage(mosdepth_dir):
    files = sorted(glob.glob(os.path.join(mosdepth_dir, "*.mosdepth.global.dist.txt")))
    if not files: return False
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, f in enumerate(files):
        d = pd.read_csv(f, sep="\t", header=None, names=["chrom", "depth", "frac"])
        d = d[d["chrom"] == "total"]
        ax.plot(d["depth"], d["frac"], lw=0.8, alpha=0.6, color=PAL[i % len(PAL)])
    ax.set_xlim(0, 80); ax.set_ylim(0, 1.02); ax.set_xlabel("depth (x)"); ax.set_ylabel("fraction of genome ≥ depth")
    ax.axvline(20, color=PAL[2], ls="--", lw=1); ax.set_title(f"Cumulative coverage distribution ({len(files)} samples)")
    save(fig, "coverage_cumulative", "Alignment", "Cumulative coverage", "Each line is one sample (mosdepth, MAPQ ≥ 20).")


@figure("Alignment", "Coverage by chromosome")
def fig_chrom_heatmap(mosdepth_dir):
    files = sorted(glob.glob(os.path.join(mosdepth_dir, "*.mosdepth.summary.txt")))
    if not files: return False
    rows, names = [], []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d = d[d["chrom"].isin(CHROMS)].set_index("chrom")["mean"]
        auto = d.reindex(AUTOS).mean()
        if not auto: continue
        rows.append((d.reindex(CHROMS) / auto).values); names.append(sample_of(f, ".mosdepth.summary.txt"))
    if not rows: return False
    m = np.array(rows, dtype=float)
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.13 * len(names))))
    im = ax.imshow(m, aspect="auto", cmap="RdBu_r", vmin=0, vmax=2)
    ax.set_xticks(range(len(CHROMS))); ax.set_xticklabels([c[3:] for c in CHROMS]); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=5)
    ax.grid(False); fig.colorbar(im, ax=ax, label="mean depth / autosomal mean"); ax.set_title("Relative coverage by chromosome")
    save(fig, "coverage_by_chromosome", "Alignment", "Coverage by chromosome",
         "Normalised to each sample's autosomal mean: chrX ≈ 1 and chrY ≈ 0 for females, ≈ 0.5 each for males; a deviating autosome suggests aneuploidy or CNV; chrM reflects mitochondrial copy number.")


@figure("Alignment", "Mapping and duplication")
def fig_mapping(qc):
    d = qc.sort_values("sample_id")
    if d["pct_mapped"].isna().all(): return False
    fig, axes = plt.subplots(3, 1, figsize=(max(7, 0.16 * len(d)), 7), sharex=True)
    for ax, col, lab, thr in ((axes[0], "pct_mapped", "mapped fraction", ARGS.min_mapped), (axes[1], "pct_properly_paired", "properly paired fraction", None),
                              (axes[2], "dup_rate", "duplicate rate", ARGS.max_dup)):
        ax.bar(d["sample_id"], pd.to_numeric(d[col], errors="coerce"), color=PAL[1]); ax.set_ylabel(lab)
        if thr is not None: ax.axhline(thr, color=PAL[2], ls="--", lw=1)
    axes[0].set_title("Mapping, pairing and duplication per sample"); axes[2].tick_params(axis="x", rotation=90, labelsize=6)
    save(fig, "mapping_duplication", "Alignment", "Mapping and duplication", "Dashed lines are the QC thresholds.")


@figure("Alignment", "Insert size")
def fig_insert(samtools_dir):
    files = sorted(glob.glob(os.path.join(samtools_dir, "*.stats")))
    if not files: return False
    fig, ax = plt.subplots(figsize=(7, 4.2)); modes = []
    for i, f in enumerate(files):
        rows = stats_section(f, "IS")
        if not rows: continue
        x = np.array([int(r[0]) for r in rows]); y = np.array([float(r[1]) for r in rows]); y = y / y.sum() if y.sum() else y
        ax.plot(x, y, lw=0.8, alpha=0.6, color=PAL[i % len(PAL)]); modes.append(x[y.argmax()])
    ax.set_xlim(0, 1000); ax.set_xlabel("insert size (bp)"); ax.set_ylabel("fraction of pairs")
    ax.set_title(f"Insert-size distributions (mode median {np.median(modes):.0f} bp)" if modes else "Insert-size distributions")
    save(fig, "insert_size", "Alignment", "Insert size", "samtools stats IS section, one line per sample.")


@figure("Alignment", "GC content")
def fig_gc(samtools_dir):
    files = sorted(glob.glob(os.path.join(samtools_dir, "*.stats")))
    if not files: return False
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, f in enumerate(files):
        rows = stats_section(f, "GCF")
        if not rows: continue
        x = np.array([float(r[0]) for r in rows]); y = np.array([float(r[1]) for r in rows]); y = y / y.sum() if y.sum() else y
        ax.plot(x, y, lw=0.8, alpha=0.6, color=PAL[i % len(PAL)])
    ax.set_xlabel("GC (%)"); ax.set_ylabel("fraction of reads"); ax.set_title("Read GC content (first reads)")
    save(fig, "gc_content", "Alignment", "GC content", "A shifted or bimodal curve points to contamination or library bias.")


@figure("Reads", "Lane yield and Q30")
def fig_lanes(fastp_dir, sheet):
    files = sorted(glob.glob(os.path.join(fastp_dir, "*.fastp.json")))
    if not files: return False
    recs = []
    for f in files:
        j = json.load(open(f)); s = j.get("summary", {})
        lane = os.path.basename(f)[:-len(".fastp.json")]; sample = lane.split(".")[0]
        recs.append({"sample": sample, "lane": lane, "reads": s.get("after_filtering", {}).get("total_reads", 0) / 1e6,
                     "q30": s.get("after_filtering", {}).get("q30_rate", np.nan), "dup": j.get("duplication", {}).get("rate", np.nan)})
    d = pd.DataFrame(recs)
    if d.empty: return False
    piv = d.pivot_table(index="sample", columns="lane", values="reads", aggfunc="sum").fillna(0)
    fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.16 * len(piv)), 6.5), sharex=True)
    bottom = np.zeros(len(piv))
    for i, lane in enumerate(piv.columns):
        axes[0].bar(piv.index, piv[lane], bottom=bottom, color=PAL[i % len(PAL)], width=0.8); bottom += piv[lane].values
    axes[0].set_ylabel("read pairs after trimming (M)"); axes[0].set_title("Yield per sample, stacked by lane")
    q = d.groupby("sample")["q30"].mean()
    axes[1].bar(q.index, q.values, color=PAL[0]); axes[1].set_ylabel("mean Q30 fraction"); axes[1].set_ylim(0.8, 1.0)
    axes[1].tick_params(axis="x", rotation=90, labelsize=6)
    save(fig, "lane_yield_q30", "Reads", "Lane yield and Q30", "fastp after-filtering totals; each colour is a lane-pair.")


@figure("Reads", "Base quality by cycle")
def fig_basequal(fastp_dir):
    files = sorted(glob.glob(os.path.join(fastp_dir, "*.fastp.json")))
    if not files: return False
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
    for f in files:
        j = json.load(open(f))
        for ax, key in ((axes[0], "read1_after_filtering"), (axes[1], "read2_after_filtering")):
            qc = j.get(key, {}).get("quality_curves", {}).get("mean")
            if qc: ax.plot(range(1, len(qc) + 1), qc, lw=0.6, alpha=0.4, color=PAL[1])
    for ax, t in zip(axes, ("read 1", "read 2")):
        ax.set_title(f"Mean base quality per cycle ({t})"); ax.set_xlabel("cycle"); ax.axhline(30, color=PAL[2], ls="--", lw=1)
    axes[0].set_ylabel("mean Phred quality")
    save(fig, "base_quality_by_cycle", "Reads", "Base quality by cycle", "One line per lane-pair (fastp, after trimming).")


@figure("Sample identity", "Contamination")
def fig_contam(qc):
    d = qc[qc["freemix"].notna()].sort_values("freemix")
    if d.empty: return False
    fig, ax = plt.subplots(figsize=(max(7, 0.16 * len(d)), 3.8))
    ax.bar(d["sample_id"], d["freemix"], color=[PAL[3] if v > ARGS.max_contam else PAL[0] for v in d["freemix"]])
    ax.axhline(ARGS.max_contam, color=PAL[2], ls="--", lw=1); ax.set_ylabel("FREEMIX"); ax.set_title("Estimated contamination (verifyBamID2)")
    ax.tick_params(axis="x", rotation=90, labelsize=6)
    save(fig, "contamination", "Sample identity", "Contamination", "FREEMIX above the dashed threshold is flagged.")


@figure("Sample identity", "Sex check")
def fig_sex(qc):
    d = qc[qc["X_ratio"].notna() & qc["Y_ratio"].notna()]
    if d.empty: return False
    fig, ax = plt.subplots(figsize=(5.5, 5))
    col = {"M": PAL[1], "F": PAL[6]}
    inferred = d["somalier_sex"].where(d["somalier_sex"].isin(["M", "F"]), d["sex_coverage"])
    for sx in ("M", "F", "other"):
        sel = d[(inferred == sx)] if sx != "other" else d[~inferred.isin(["M", "F"])]
        if sel.empty: continue
        ax.scatter(sel["X_ratio"], sel["Y_ratio"], s=28, color=col.get(sx, PAL[2]), label=f"inferred {sx} (n={len(sel)})", alpha=0.85)
    mis = d[d["qc_flags"].fillna("").str.contains("sex_mismatch|sex_ambiguous")]
    for _, r in mis.iterrows():
        ax.annotate(r["sample_id"], (r["X_ratio"], r["Y_ratio"]), fontsize=6, xytext=(3, 3), textcoords="offset points", color=PAL[3])
    ax.set_xlabel("chrX / autosomal coverage"); ax.set_ylabel("chrY / autosomal coverage"); ax.set_title("Sex inference from coverage"); ax.legend()
    save(fig, "sex_check", "Sample identity", "Sex check", "Labelled points disagree with the sample sheet or are ambiguous.")


@figure("Sample identity", "Relatedness")
def fig_relatedness(pairs_file):
    if not pairs_file or not os.path.isfile(pairs_file): return False
    d = pd.read_csv(pairs_file, sep="\t")
    d.columns = [c.lstrip("#") for c in d.columns]
    if "relatedness" not in d: return False
    samples = sorted(set(d["sample_a"]) | set(d["sample_b"])); idx = {s: i for i, s in enumerate(samples)}
    m = np.full((len(samples), len(samples)), np.nan)
    for _, r in d.iterrows(): m[idx[r["sample_a"]], idx[r["sample_b"]]] = m[idx[r["sample_b"]], idx[r["sample_a"]]] = r["relatedness"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1, 1.4]})
    axes[0].hist(d["relatedness"].clip(-0.5, 1), bins=60, color=PAL[0]); axes[0].set_xlabel("relatedness"); axes[0].set_ylabel("pairs")
    axes[0].axvline(0.2, color=PAL[2], ls="--", lw=1); axes[0].set_title("Pairwise relatedness (somalier)")
    im = axes[1].imshow(m, cmap="viridis", vmin=-0.2, vmax=1); axes[1].grid(False)
    axes[1].set_xticks(range(len(samples))); axes[1].set_xticklabels(samples, rotation=90, fontsize=4); axes[1].set_yticks(range(len(samples))); axes[1].set_yticklabels(samples, fontsize=4)
    fig.colorbar(im, ax=axes[1], label="relatedness"); axes[1].set_title("Relatedness matrix")
    hi = d[d["relatedness"] >= 0.2]
    save(fig, "relatedness", "Sample identity", "Relatedness", f"{len(hi)} pair(s) at or above 0.2 (duplicates ≈ 1, first-degree ≈ 0.5, second-degree ≈ 0.25).")


@figure("Population structure", "PCA")
def fig_pca(eigenvec, eigenval, sheet):
    if not eigenvec or not os.path.isfile(eigenvec): return False
    ev = pd.read_csv(eigenvec, sep=r"\s+")
    ev.columns = [c.lstrip("#") for c in ev.columns]
    idcol = "IID" if "IID" in ev else ev.columns[0]
    pcs = [c for c in ev.columns if c.startswith("PC")]
    if len(pcs) < 2: return False
    var = None
    if eigenval and os.path.isfile(eigenval):
        vals = np.loadtxt(eigenval); var = vals / vals.sum() * 100
    pheno = sheet.set_index("sample_id")["phenotype"].to_dict() if sheet is not None and "phenotype" in sheet else {}
    groups = ev[idcol].map(lambda s: pheno.get(s) or "all")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (a, b) in zip(axes[:2], ((pcs[0], pcs[1]), (pcs[2], pcs[3]) if len(pcs) >= 4 else (pcs[0], pcs[1]))):
        for i, g in enumerate(sorted(groups.unique())):
            sel = ev[groups == g]; ax.scatter(sel[a], sel[b], s=22, alpha=0.8, color=PAL[i % len(PAL)], label=f"{g} (n={len(sel)})")
        lab = lambda p: f"{p} ({var[int(p[2:]) - 1]:.1f}%)" if var is not None and int(p[2:]) <= len(var) else p
        ax.set_xlabel(lab(a)); ax.set_ylabel(lab(b)); ax.legend(fontsize=7)
    axes[0].set_title("PCA on LD-pruned common SNPs")
    if var is not None:
        axes[2].bar(range(1, len(var) + 1), var, color=PAL[1]); axes[2].set_xlabel("PC"); axes[2].set_ylabel("variance explained (%)"); axes[2].set_title("Scree")
    save(fig, "pca", "Population structure", "PCA", "plink2 PCA; colours are the phenotype column of the sample sheet when present.")


@figure("Variants", "Per-sample variant metrics")
def fig_psc(bcf_stats):
    if not bcf_stats or not os.path.isfile(bcf_stats): return False
    rows = stats_section(bcf_stats, "PSC")
    if not rows: return False
    d = pd.DataFrame([{"sample": r[1], "nRefHom": int(r[2]), "nNonRefHom": int(r[3]), "nHets": int(r[4]), "nTs": int(r[5]), "nTv": int(r[6]), "nIndels": int(r[7]), "nSingletons": int(r[9])} for r in rows]).sort_values("sample")
    d["nSNV"] = d["nNonRefHom"] + d["nHets"]; d["het_hom"] = d["nHets"] / d["nNonRefHom"].replace(0, np.nan); d["titv"] = d["nTs"] / d["nTv"].replace(0, np.nan)
    fig, axes = plt.subplots(2, 2, figsize=(max(9, 0.22 * len(d)), 7), sharex=True)
    for ax, col, lab in zip(axes.flat, ("nSNV", "nIndels", "het_hom", "titv"), ("non-reference SNVs", "indels", "het / hom-alt", "Ti/Tv")):
        ax.bar(d["sample"], d[col], color=PAL[0]); ax.set_ylabel(lab)
    axes[1, 1].axhline(2.0, color=PAL[2], ls="--", lw=1); axes[1, 0].tick_params(axis="x", rotation=90, labelsize=6); axes[1, 1].tick_params(axis="x", rotation=90, labelsize=6)
    axes[0, 0].set_title("Per-sample variant metrics (PASS, normalised)")
    save(fig, "variants_per_sample", "Variants", "Per-sample variant metrics", "Whole-genome Ti/Tv is expected around 2.0–2.1; outliers suggest excess false positives or contamination.")


@figure("Variants", "Site statistics")
def fig_site_stats(bcf_stats):
    if not bcf_stats or not os.path.isfile(bcf_stats): return False
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    idd = stats_section(bcf_stats, "IDD")
    if idd:
        x = [int(r[1]) for r in idd]; y = [int(r[2]) for r in idd]
        axes[0, 0].bar(x, y, color=PAL[1], width=1); axes[0, 0].set_xlim(-40, 40); axes[0, 0].set_xlabel("indel length (bp)"); axes[0, 0].set_ylabel("sites"); axes[0, 0].set_title("Indel length distribution"); axes[0, 0].set_yscale("log")
    st = stats_section(bcf_stats, "ST")
    if st:
        axes[0, 1].bar([r[1] for r in st], [int(r[2]) for r in st], color=PAL[0]); axes[0, 1].set_title("Substitution types"); axes[0, 1].tick_params(axis="x", rotation=90)
    qual = stats_section(bcf_stats, "QUAL")
    if qual:
        q = np.array([float(r[1]) for r in qual]); n = np.array([int(r[2]) + int(r[5]) for r in qual])
        axes[1, 0].plot(q, n, color=PAL[1]); axes[1, 0].set_xscale("log"); axes[1, 0].set_xlabel("QUAL"); axes[1, 0].set_ylabel("sites"); axes[1, 0].set_title("Site quality")
    dp = stats_section(bcf_stats, "DP")
    if dp:
        x = [int(r[1]) if r[1].isdigit() else 500 for r in dp]; y = [int(r[4]) for r in dp]
        axes[1, 1].plot(x, y, color=PAL[1]); axes[1, 1].set_xlim(0, 300); axes[1, 1].set_xlabel("site depth"); axes[1, 1].set_ylabel("sites"); axes[1, 1].set_title("Site depth distribution")
    save(fig, "site_statistics", "Variants", "Site statistics", "bcftools stats on the PASS, normalised cohort VCF.")


@figure("Variants", "VQSR tranches")
def fig_tranches(tr_dir):
    files = [f for f in glob.glob(os.path.join(tr_dir or ".", "*.tranches")) if os.path.getsize(f) > 0]
    if not files: return False
    fig, axes = plt.subplots(1, len(files), figsize=(5.5 * len(files), 4), squeeze=False)
    for ax, f in zip(axes[0], sorted(files)):
        d = pd.read_csv(f, comment="#")
        if "targetTruthSensitivity" not in d: continue
        d = d.sort_values("targetTruthSensitivity")
        ax.bar(d["targetTruthSensitivity"].astype(str), d["numNovel"], color=PAL[2], label="novel"); ax.bar(d["targetTruthSensitivity"].astype(str), d["numKnown"], bottom=d["numNovel"], color=PAL[1], label="known")
        ax2 = ax.twinx(); ax2.plot(d["targetTruthSensitivity"].astype(str), d["novelTiTv"], color=PAL[3], marker="o", label="novel Ti/Tv"); ax2.set_ylabel("novel Ti/Tv"); ax2.grid(False)
        ax.set_xlabel("truth sensitivity tranche (%)"); ax.set_ylabel("variants (cumulative)"); ax.set_title(os.path.basename(f).replace(".tranches", "") + " tranches"); ax.legend(loc="upper left", fontsize=7); ax.tick_params(axis="x", rotation=45)
    save(fig, "vqsr_tranches", "Variants", "VQSR tranches", "Variants retained at each truth-sensitivity level; novel Ti/Tv falls as lower-quality tranches are admitted.")


@figure("Annotation", "Allele frequency and annotation summary")
def fig_annotation(flat):
    if not flat or not os.path.isfile(flat): return False
    usecols = lambda c: re.sub(r"^#?\[\d+\]", "", c).strip() in ("CHROM", "ID", "AF", "Consequence", "IMPACT", "gnomAD_AF", "CLNSIG")
    af, known, cons, imp, cln, gn_both, n_novel, n_total = [], [], Counter(), Counter(), Counter(), [], 0, 0
    for chunk in pd.read_csv(flat, sep="\t", usecols=usecols, chunksize=500_000, dtype=str, na_values=["."], low_memory=False):
        chunk.columns = [c.split("]")[-1].lstrip("#").strip() for c in chunk.columns]
        a = pd.to_numeric(chunk.get("AF"), errors="coerce"); af.append(a.dropna().values)
        isk = chunk.get("ID", pd.Series(dtype=str)).fillna("").str.startswith("rs"); known.append(isk.values); n_total += len(chunk)
        if "Consequence" in chunk: cons.update(chunk["Consequence"].fillna("").str.split("&").str[0].value_counts().to_dict())
        if "IMPACT" in chunk: imp.update(chunk["IMPACT"].fillna("").value_counts().to_dict())
        if "CLNSIG" in chunk: cln.update(chunk["CLNSIG"].dropna().str.split("[|,/]").str[0].value_counts().to_dict())
        if "gnomAD_AF" in chunk:
            g = pd.to_numeric(chunk["gnomAD_AF"], errors="coerce"); m = g.notna() & a.notna(); n_novel += int((~g.notna()).sum())
            if m.any(): gn_both.append(np.column_stack([g[m].values, a[m].values]))
    af = np.concatenate(af) if af else np.array([]); known = np.concatenate(known) if known else np.array([])
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    if af.size:
        bins = np.linspace(0, 1, 41)
        axes[0, 0].hist(af[known[:len(af)]] if known.size == af.size else af, bins=bins, color=PAL[1], label="in dbSNP", alpha=0.9)
        if known.size == af.size: axes[0, 0].hist(af[~known], bins=bins, color=PAL[2], label="novel", alpha=0.8)
        axes[0, 0].set_yscale("log"); axes[0, 0].set_xlabel("cohort allele frequency"); axes[0, 0].set_ylabel("variants"); axes[0, 0].set_title(f"Site frequency spectrum ({n_total:,} alleles)"); axes[0, 0].legend()
    if cons:
        top = cons.most_common(18)[::-1]
        axes[0, 1].barh([k or "(none)" for k, _ in top], [v for _, v in top], color=PAL[0]); axes[0, 1].set_xscale("log"); axes[0, 1].set_title("Worst consequence per variant (VEP)")
    if gn_both:
        gb = np.concatenate(gn_both); h = axes[1, 0].hist2d(gb[:, 0], gb[:, 1], bins=50, cmap="viridis", norm=matplotlib.colors.LogNorm())
        axes[1, 0].plot([0, 1], [0, 1], color="white", lw=0.8, ls="--"); axes[1, 0].set_xlabel("gnomAD AF"); axes[1, 0].set_ylabel("cohort AF"); axes[1, 0].grid(False)
        axes[1, 0].set_title(f"Cohort vs gnomAD AF ({len(gb):,} shared; {n_novel:,} absent from gnomAD)"); fig.colorbar(h[3], ax=axes[1, 0])
    if cln or imp:
        items = [(k, v) for k, v in cln.most_common(10)] or []
        if items:
            axes[1, 1].barh([k for k, _ in items][::-1], [v for _, v in items][::-1], color=PAL[3]); axes[1, 1].set_xscale("log"); axes[1, 1].set_title("ClinVar clinical significance (annotated variants)")
        else:
            axes[1, 1].bar(list(imp), list(imp.values()), color=PAL[1]); axes[1, 1].set_title("VEP impact")
    save(fig, "annotation_summary", "Annotation", "Allele frequency and annotation summary", "From variants.flat.tsv.gz (worst-consequence transcript per variant).")


# ----------------------------------------------------------------------------- dashboard
def dashboard(qc, out_html, run_label):
    def img(name):
        with open(os.path.join(OUT, name), "rb") as fh: return "data:image/png;base64," + base64.b64encode(fh.read()).decode()
    sections = defaultdict(list)
    for sec, title, fn, cap in FIGS: sections[sec].append((title, fn, cap))
    flagged = qc[qc["qc_status"] != "PASS"] if "qc_status" in qc else pd.DataFrame()
    cols = [c for c in ("sample_id", "phenotype", "n_lanes", "mean_coverage", "pct_20x", "dup_rate", "pct_mapped", "freemix", "sex_sheet", "somalier_sex", "titv", "het_hom_ratio", "related_to", "qc_flags", "qc_status") if c in qc]
    h = [f"<!doctype html><html><head><meta charset='utf-8'><title>GenomeForge QC — {run_label}</title><style>"
         "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f7f8;color:#1c1e21}"
         "header{background:#15191c;color:#f2ece2;padding:22px 36px}header h1{margin:0;font-size:22px}header span{color:#74c4b4}"
         "main{max-width:1200px;margin:0 auto;padding:20px 36px}h2{border-bottom:2px solid #2FA38E;padding-bottom:4px;margin-top:36px}"
         "figure{background:#fff;border:1px solid #e3e5e8;border-radius:8px;padding:12px;margin:14px 0}figure img{max-width:100%}figcaption{color:#555;font-size:13px;margin-top:6px}"
         "table{border-collapse:collapse;font-size:12px;background:#fff}th,td{border:1px solid #e3e5e8;padding:4px 7px;text-align:left}th{background:#eef1f3}tr.flag td{background:#fdecea}"
         ".kpi{display:inline-block;background:#fff;border:1px solid #e3e5e8;border-radius:8px;padding:10px 16px;margin:6px 8px 6px 0}.kpi b{font-size:20px;display:block;color:#2FA38E}"
         "nav a{margin-right:14px;color:#2FA38E}</style></head><body>",
         f"<header><h1>Genome<span>Forge</span> cohort QC dashboard</h1><div>{run_label}</div></header><main>",
         "<nav>" + "".join(f"<a href='#{s.replace(' ', '_')}'>{s}</a>" for s in sections) + "<a href='#samples'>Samples</a></nav>"]
    n = len(qc); med = lambda c: f"{pd.to_numeric(qc[c], errors='coerce').median():.3g}" if c in qc else "n/a"
    h.append(f"<div><div class='kpi'><b>{n}</b>samples</div><div class='kpi'><b>{len(flagged)}</b>flagged</div><div class='kpi'><b>{med('mean_coverage')}x</b>median coverage</div>"
             f"<div class='kpi'><b>{med('dup_rate')}</b>median dup rate</div><div class='kpi'><b>{med('titv')}</b>median Ti/Tv</div></div>")
    for sec, items in sections.items():
        h.append(f"<h2 id='{sec.replace(' ', '_')}'>{sec}</h2>")
        for title, fn, cap in items:
            h.append(f"<figure><b>{title}</b><br><img src='{img(fn)}' alt='{title}'><figcaption>{cap} <i>({fn})</i></figcaption></figure>")
    h.append("<h2 id='samples'>Samples</h2><table><tr>" + "".join(f"<th>{c}</th>" for c in cols) + "</tr>")
    for _, r in qc.sort_values("sample_id").iterrows():
        h.append(f"<tr class='{'flag' if r.get('qc_status') != 'PASS' else ''}'>" + "".join(f"<td>{'' if pd.isna(r[c]) else r[c]}</td>" for c in cols) + "</tr>")
    h.append("</table>")
    if NOTES: h.append("<h2>Notes</h2><ul>" + "".join(f"<li>{n}</li>" for n in NOTES) + "</ul>")
    h.append("</main></body></html>")
    with open(out_html, "w") as fh: fh.write("\n".join(h))


def pdf_report(qc, out_pdf, run_label):
    """Multi-page PDF: title + KPIs, flagged samples, sample table, then one figure per page."""
    flagged = qc[qc["qc_status"] != "PASS"] if "qc_status" in qc else pd.DataFrame()
    med = lambda c: f"{pd.to_numeric(qc[c], errors='coerce').median():.3g}" if c in qc else "n/a"
    with PdfPages(out_pdf) as pdf:
        # ---- title page
        fig = plt.figure(figsize=(8.27, 11.69)); fig.patch.set_facecolor("white")
        fig.text(0.08, 0.90, "GenomeForge", fontsize=30, weight="bold", color="#15191c")
        fig.text(0.08, 0.86, "Cohort QC report", fontsize=18, color="#2FA38E")
        fig.text(0.08, 0.82, run_label, fontsize=9, color="#555")
        kpis = [("samples", str(len(qc))), ("flagged", str(len(flagged))), ("median coverage", med("mean_coverage") + "x"),
                ("median duplicate rate", med("dup_rate")), ("median mapped fraction", med("pct_mapped")),
                ("median FREEMIX", med("freemix")), ("median Ti/Tv", med("titv")), ("median het/hom", med("het_hom_ratio"))]
        for i, (k, v) in enumerate(kpis):
            y = 0.74 - i * 0.045
            fig.text(0.10, y, k, fontsize=11, color="#333"); fig.text(0.55, y, v, fontsize=11, weight="bold", color="#15191c")
        fig.text(0.08, 0.36, "Thresholds", fontsize=12, weight="bold")
        fig.text(0.10, 0.32, f"mean coverage ≥ {ARGS.min_depth}x · duplicate rate ≤ {ARGS.max_dup} · FREEMIX ≤ {ARGS.max_contam} · mapped ≥ {ARGS.min_mapped}", fontsize=9)
        fig.text(0.08, 0.26, "Contents", fontsize=12, weight="bold")
        secs = []
        for sec, _, _, _ in FIGS:
            if sec not in secs: secs.append(sec)
        for i, sec in enumerate(secs):
            fig.text(0.10, 0.22 - i * 0.03, f"{i + 2}. {sec}", fontsize=10)
        if NOTES:
            fig.text(0.08, 0.22 - len(secs) * 0.03 - 0.03, "Skipped figures: " + "; ".join(n.split(":")[0] for n in NOTES), fontsize=7, color="#9E2B25", wrap=True)
        pdf.savefig(fig); plt.close(fig)
        # ---- flagged samples + full table (paginated)
        cols = [c for c in ("sample_id", "n_lanes", "mean_coverage", "pct_20x", "dup_rate", "pct_mapped", "freemix", "sex_sheet", "somalier_sex", "titv", "qc_status") if c in qc]
        tab = qc.sort_values("sample_id")[cols].copy()
        for c in tab.columns:
            if tab[c].dtype.kind == "f": tab[c] = tab[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
        tab = tab.fillna("")
        per_page = 38
        for start in range(0, len(tab), per_page):
            fig, ax = plt.subplots(figsize=(8.27, 11.69)); ax.axis("off")
            ax.set_title("Sample QC table" + (f" ({start + 1}–{min(start + per_page, len(tab))} of {len(tab)})" if len(tab) > per_page else ""), loc="left", fontsize=12)
            chunk = tab.iloc[start:start + per_page]
            t = ax.table(cellText=chunk.values, colLabels=chunk.columns, loc="upper left", cellLoc="left")
            t.auto_set_font_size(False); t.set_fontsize(6); t.scale(1, 1.15)
            for (r, c), cell in t.get_celld().items():
                cell.set_edgecolor("#dddddd")
                if r == 0: cell.set_facecolor("#eef1f3"); cell.set_text_props(weight="bold")
                elif chunk.iloc[r - 1]["qc_status"] != "PASS": cell.set_facecolor("#fdecea")
            if start == 0 and len(flagged):
                ax.text(0, -0.02, "Flagged: " + "; ".join(f"{r['sample_id']} [{r['qc_flags']}]" for _, r in flagged.iterrows())[:1500],
                        transform=ax.transAxes, fontsize=6, va="top", wrap=True, color="#9E2B25")
            pdf.savefig(fig); plt.close(fig)
        # ---- one figure per page
        for sec, title, fn, cap in FIGS:
            img = mpimg.imread(os.path.join(OUT, fn))
            h, w = img.shape[:2]
            fig = plt.figure(figsize=(11.69, 8.27) if w > h * 1.15 else (8.27, 11.69))
            fig.text(0.05, 0.96, f"{sec} — {title}", fontsize=13, weight="bold", va="top")
            ax = fig.add_axes([0.04, 0.12, 0.92, 0.80]); ax.imshow(img); ax.axis("off")
            fig.text(0.05, 0.07, cap, fontsize=8, color="#444", wrap=True, va="top")
            fig.text(0.95, 0.02, fn, fontsize=6, color="#999", ha="right")
            pdf.savefig(fig); plt.close(fig)
        d = pdf.infodict(); d["Title"] = "GenomeForge cohort QC report"; d["Subject"] = run_label


def main():
    global ARGS, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-summary", required=True); ap.add_argument("--samplesheet")
    ap.add_argument("--mosdepth", default=""); ap.add_argument("--samtools", default=""); ap.add_argument("--fastp", default="")
    ap.add_argument("--somalier-pairs"); ap.add_argument("--pca-eigenvec"); ap.add_argument("--pca-eigenval")
    ap.add_argument("--bcftools-stats"); ap.add_argument("--tranches", default=""); ap.add_argument("--flat")
    ap.add_argument("--min-depth", type=float, default=20); ap.add_argument("--max-dup", type=float, default=0.25)
    ap.add_argument("--max-contam", type=float, default=0.03); ap.add_argument("--min-mapped", type=float, default=0.95)
    ap.add_argument("--out-dir", default="figures"); ap.add_argument("--dashboard", default="qc_dashboard.html")
    ap.add_argument("--pdf", default="qc_report.pdf"); ap.add_argument("--label", default="")
    ARGS = ap.parse_args(); OUT = ARGS.out_dir
    qc = pd.read_csv(ARGS.qc_summary, sep="\t", dtype=str)
    for c in ("mean_coverage", "pct_20x", "dup_rate", "pct_mapped", "pct_properly_paired", "freemix", "X_ratio", "Y_ratio", "titv", "het_hom_ratio"):
        if c in qc: qc[c] = pd.to_numeric(qc[c], errors="coerce")
    sheet = pd.read_csv(ARGS.samplesheet, dtype=str).drop_duplicates("sample_id") if ARGS.samplesheet and os.path.isfile(ARGS.samplesheet) else None
    fig_lanes(ARGS.fastp, sheet); fig_basequal(ARGS.fastp)
    fig_coverage(qc); fig_cum_coverage(ARGS.mosdepth); fig_chrom_heatmap(ARGS.mosdepth); fig_mapping(qc); fig_insert(ARGS.samtools); fig_gc(ARGS.samtools)
    fig_contam(qc); fig_sex(qc); fig_relatedness(ARGS.somalier_pairs)
    fig_pca(ARGS.pca_eigenvec, ARGS.pca_eigenval, sheet)
    fig_psc(ARGS.bcftools_stats); fig_site_stats(ARGS.bcftools_stats); fig_tranches(ARGS.tranches); fig_annotation(ARGS.flat)
    dashboard(qc, ARGS.dashboard, ARGS.label)
    try:
        pdf_report(qc, ARGS.pdf, ARGS.label)
    except Exception as e:
        NOTES.append(f"PDF report: failed ({type(e).__name__}: {e})")
    print(f"[plot_qc] {len(FIGS)} figures -> {OUT}; dashboard -> {ARGS.dashboard}; pdf -> {ARGS.pdf}" + (f"; skipped: {len(NOTES)}" if NOTES else ""))
    for n in NOTES: print("  - " + n)


if __name__ == "__main__":
    main()
