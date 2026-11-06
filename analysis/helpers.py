# Imports
from bed_reader import open_bed
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import stats



# Directory variables
ROOT = Path(__file__).resolve().parent.parent

LDREF_DIR = ROOT / "data" / "genotype" / "LDREF"
EXPRESSION_FILE = ROOT / "data" / "phenotype" / "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
GENE_ANNOT_FILE = ROOT / "data" / "annotation" / "gene_annot.txt.gz"
HEIGHT_GWAS_FILE = ROOT / "data" / "gwas" / "Height.gwas.txt.gz"

OUT_EQTL = ROOT / "output" / "part1_eqtl"
OUT_PRS = ROOT / "output" / "part2_prs"

PLINK2 = ROOT / "tools" / ("plink2.exe" if sys.platform.startswith("win") else "plink2")



# Helper functions
def read_expression(path=EXPRESSION_FILE):
    return pd.read_csv(path, sep="\t", compression="infer")


def read_gene_annotation(path=GENE_ANNOT_FILE):
    return pd.read_csv(path, sep="\t", compression="infer")


def read_bim(prefix):
    return pd.read_csv(f"{prefix}.bim", sep="\t", header=None, names=["chr", "snp", "cm", "pos", "a1", "a2"])


def read_fam(prefix):
    return pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, names=["fid", "iid", "pat", "mat", "sex", "pheno"])


def read_genotypes(prefix):
    bed = open_bed(f"{prefix}.bed")
    geno = bed.read(dtype="float32")
    iids = [str(s) for s in bed.iid]
    return geno, iids, read_bim(prefix)


def read_gwas_sumstats(path=HEIGHT_GWAS_FILE):
    return pd.read_csv(path, sep=r"\s+", compression="infer")


def read_sscore(path):
    return pd.read_csv(path, sep="\t").rename(columns={"#FID": "FID", "#IID": "IID"})


def ldref_prefix(chrom):
    return LDREF_DIR / f"1000G.EUR.{chrom}"


def plink_run(args, *, check=True, quiet=True):
    if not PLINK2.exists():
        raise FileNotFoundError(
            f"plink2 not found at {PLINK2}. Place the executable in {ROOT / 'tools'}."
        )
    out = subprocess.run([str(PLINK2), *map(str, args)], capture_output=True, text=True, check=check)
    if not quiet:
        print(out.stdout)
        if out.stderr:
            print(out.stderr)
    return out


def plink_clump(bfile, sumstats, out, *, snp_field="SNP", p_field="P", r2=0.1, kb=250, p1=1.0, extract=None):
    args = ["--bfile", bfile, "--clump", sumstats, "--clump-snp-field", snp_field, "--clump-field", p_field, "--clump-p1", p1, "--clump-r2", r2, "--clump-kb", kb, "--out", out]
    if extract is not None:
        args += ["--extract", extract]
    plink_run(args)
    return out.with_suffix(".clumps")


def plink_score(bfile, score_file, out, *, snp_col, a1_col, beta_col, q_range, q_values, extract=None):
    args = ["--bfile", bfile, "--score", score_file, snp_col, a1_col, beta_col, "header", "--q-score-range", q_range, q_values, "--out", out]
    if extract is not None:
        args += ["--extract", extract]
    plink_run(args)


def vectorized_ols(G, y):
    G = np.asarray(G, dtype=np.float64).copy()
    y = np.asarray(y, dtype=np.float64)

    col_means = np.nanmean(G, axis=0)
    nan_mask = np.isnan(G)
    if nan_mask.any():
        G[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    keep = np.where(G.var(axis=0) > 0)[0]
    if keep.size == 0:
        return pd.DataFrame(columns=["snp_idx", "beta", "se", "pval"])
    G = G[:, keep]

    n = G.shape[0]
    y_c = y - y.mean()
    G_c = G - G.mean(axis=0, keepdims=True)
    GtG = (G_c ** 2).sum(axis=0)
    Gty = G_c.T @ y_c
    beta = Gty / GtG
    rss = np.maximum(0.0, (y_c ** 2).sum() - beta * Gty)
    se = np.sqrt((rss / (n - 2)) / GtG)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = beta / se
    pval = 2 * stats.t.sf(np.abs(tstat), df=n - 2)
    return pd.DataFrame({"snp_idx": keep, "beta": beta, "se": se, "pval": pval})


def cis_snp_indices(bim, start, stop, window=500_000):
    pos = bim["pos"].to_numpy()
    return np.where((pos >= start - window) & (pos <= stop + window))[0]


def eqtl_for_gene(geno, bim, expr_vec, start, stop, *, gene= None, window=500_000):
    cols = cis_snp_indices(bim, start, stop, window)
    if cols.size == 0:
        return pd.DataFrame()
    res = vectorized_ols(geno[:, cols], expr_vec)
    if res.empty:
        return res
    real = cols[res["snp_idx"].to_numpy()]
    out = pd.DataFrame({
        "gene": gene,
        "snp": bim["snp"].to_numpy()[real],
        "chr": bim["chr"].to_numpy()[real],
        "pos": bim["pos"].to_numpy()[real],
        "a1": bim["a1"].to_numpy()[real],
        "beta": res["beta"].to_numpy(),
        "se": res["se"].to_numpy(),
        "pval": res["pval"].to_numpy(),
    })
    out["neglog10p"] = -np.log10(out["pval"])
    return out


def _impute(G, means):
    G = G.copy()
    nan = np.isnan(G)
    if nan.any():
        G[nan] = np.take(means, np.where(nan)[1])
    return G


def _ols_beta_p(G, y):
    n = G.shape[0]
    beta = np.full(G.shape[1], np.nan)
    pval = np.full(G.shape[1], np.nan)
    ok = G.var(axis=0) > 0
    if ok.any():
        Gc = G[:, ok] - G[:, ok].mean(axis=0)
        yc = y - y.mean()
        GtG = (Gc ** 2).sum(axis=0)
        Gty = Gc.T @ yc
        b = Gty / GtG
        rss = np.maximum(0.0, (yc ** 2).sum() - b * Gty)
        se = np.sqrt((rss / (n - 2)) / GtG)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = b / se
        beta[ok] = b
        pval[ok] = 2 * stats.t.sf(np.abs(t), df=n - 2)
    return beta, pval


def _greedy_clump(G, pos, pval, r2=0.1, kb=250):
    cand = np.where(np.isfinite(pval))[0]
    order = cand[np.argsort(pval[cand])]
    Gc = G - G.mean(axis=0)
    sd = Gc.std(axis=0)
    removed = np.zeros(G.shape[1], dtype=bool)
    kept = []
    span, n = kb * 1000, len(G)
    for i in order:
        if removed[i] or sd[i] == 0:
            continue
        kept.append(i)
        removed[i] = True
        for j in np.where(np.abs(pos - pos[i]) <= span)[0]:
            if removed[j] or sd[j] == 0:
                continue
            r = (Gc[:, i] @ Gc[:, j]) / (n * sd[i] * sd[j])
            if r * r > r2:
                removed[j] = True
    return kept


def cross_validate_gene(geno, bim, expr_vec, start, stop, *, thresholds=None, window=500_000, splits=(0.6, 0.3, 0.1), seed=42, null=True):
    if thresholds is None:
        thresholds = [0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

    cols = cis_snp_indices(bim, start, stop, window)
    if cols.size == 0:
        return None

    y = np.asarray(expr_vec, dtype=float)
    keep_ind = np.isfinite(y)
    G = geno[np.ix_(keep_ind, cols)].astype(float)
    y = y[keep_ind]
    pos = bim["pos"].to_numpy()[cols]
    n = len(y)
    if n < 20:
        return None

    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    n_tr, n_va = int(splits[0] * n), int(splits[1] * n)
    tr, va, te = idx[:n_tr], idx[n_tr:n_tr + n_va], idx[n_tr + n_va:]

    train_means = np.nanmean(G[tr], axis=0)
    train_means = np.where(np.isfinite(train_means), train_means, 0.0)

    beta, pval = _ols_beta_p(_impute(G[tr], train_means), y[tr])
    kept = _greedy_clump(_impute(G[tr], train_means), pos, pval)

    def r2_at(th, betas, subset):
        use = [j for j in kept if np.isfinite(pval[j]) and pval[j] <= th]
        if not use:
            return np.nan
        prs = _impute(G[np.ix_(subset, use)], train_means[use]) @ betas[use]
        yy = y[subset]
        m = np.isfinite(prs) & np.isfinite(yy)
        if m.sum() < 5 or np.std(prs[m]) == 0:
            return np.nan
        return float(np.corrcoef(prs[m], yy[m])[0, 1] ** 2)

    def best_on_val(betas):
        vals = np.array([r2_at(t, betas, va) for t in thresholds], dtype=float)
        if not np.any(np.isfinite(vals)):
            return 0, np.nan
        i = int(np.nanargmax(np.where(np.isfinite(vals), vals, -np.inf)))
        return i, vals[i]

    bi, val_r2 = best_on_val(beta)
    out = {"n_cis": int(cols.size), "n_clumped": len(kept),
           "best_threshold": thresholds[bi], "val_r2": val_r2,
           "test_r2": r2_at(thresholds[bi], beta, te)}

    if null:
        b_null = beta.copy()
        fin = np.where(np.isfinite(b_null))[0]
        b_null[fin] = rng.permutation(b_null[fin])
        ni, _ = best_on_val(b_null)
        out["null_test_r2"] = r2_at(thresholds[ni], b_null, te)

    return out


def write_range_list(path, thresholds=[0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]):
    Path(path).write_text("\n".join(f"{t} 0 {t}" for t in thresholds) + "\n")
    return Path(path)


def qc_sumstats(gwas, *, info_min=0.8, maf_min=0.01):
    d = gwas.copy()
    if "BETA" not in d.columns and "OR" in d.columns:
        d["BETA"] = np.log(d["OR"])
    a1, a2 = d["A1"].str.upper(), d["A2"].str.upper()
    d = d[~(a1.map({"A": "T", "T": "A", "C": "G", "G": "C"}) == a2)]
    if "INFO" in d.columns:
        d = d[d["INFO"] > info_min]
    if "MAF" in d.columns:
        d = d[d["MAF"] > maf_min]
    return d[~d["SNP"].duplicated()].reset_index(drop=True)


def clump_per_chromosome(sumstats_file, out_dir, *, extract=None, prefix="clump"):
    index_snps = []
    for chrom in range(1, 23):
        bfile = ldref_prefix(chrom)
        if not Path(f"{bfile}.bed").exists():
            continue
        clumps = plink_clump(bfile, sumstats_file, out_dir / f"{prefix}_chr{chrom}", extract=extract)
        if clumps.exists():
            index_snps.extend(pd.read_csv(clumps, sep=r"\s+")["ID"].tolist())
    return index_snps


def aggregate_scores(out_dir, pattern, threshold):
    total = None
    for chrom in range(1, 23):
        f = out_dir / f"{pattern}_chr{chrom}.{threshold}.sscore"
        if not f.exists():
            continue
        sc = read_sscore(f).assign(s=lambda x: x["SCORE1_AVG"] * x["ALLELE_CT"])[["IID", "s"]]
        if total is None:
            total = sc.rename(columns={"s": "prs"})
        else:
            total = total.merge(sc, on="IID", how="outer")
            total["prs"] = total["prs"].fillna(0) + total["s"].fillna(0)
            total = total.drop(columns="s")
    return total if total is not None else pd.DataFrame(columns=["IID", "prs"])