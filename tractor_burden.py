#!/usr/bin/env python3

import argparse
import logging
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import statsmodels.api as sm

logger = logging.getLogger("tractor_burden")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

META_COLS = ["CHROM", "POS", "ID", "REF", "ALT"]
ID_COLS = ["CHROM", "POS", "REF", "ALT"]


def norm_chrom(x):
    x = str(x).strip()
    return x if x.startswith("chr") else "chr" + x


def normalize_variant_df(df):
    df = df.copy()
    if "#CHROM" in df.columns and "CHROM" not in df.columns:
        df = df.rename(columns={"#CHROM": "CHROM"})
    df["CHROM"] = df["CHROM"].map(norm_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")
    df["REF"] = df["REF"].astype(str).str.strip()
    df["ALT"] = df["ALT"].astype(str).str.strip()
    return df


def make_key(row):
    return (
        norm_chrom(row["CHROM"]),
        int(row["POS"]),
        str(row["REF"]).strip(),
        str(row["ALT"]).strip(),
    )


def parse_ann_string(x):
    if pd.isna(x):
        return []
    s = str(x)
    for sep in ["&", ";", "|"]:
        s = s.replace(sep, ",")
    return [a.strip() for a in s.split(",") if a.strip()]


def parse_region_token(tok):
    tok = tok.strip()
    m = re.match(r"^(chr)?([^:]+):(\d+)-(\d+)$", tok)
    if m is None:
        return None

    chrom = norm_chrom(m.group(2))
    start = int(m.group(3))
    end = int(m.group(4))

    if start > end:
        start, end = end, start

    return chrom, start, end


def load_set_file(set_file):
    intervals_by_chrom = defaultdict(list)
    gene_chrom = {}

    with open(set_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)
            gene = parts[0]

            # RVTESTS set files can have comma-separated intervals in one field
            region_tokens = []
            for field in parts[1:]:
                region_tokens.extend([x for x in field.split(",") if x.strip()])

            for tok in region_tokens:
                parsed = parse_region_token(tok)
                if parsed is None:
                    continue

                chrom, start, end = parsed
                intervals_by_chrom[chrom].append((start, end, gene))
                gene_chrom[gene] = chrom.replace("chr", "")

    for chrom in intervals_by_chrom:
        intervals_by_chrom[chrom].sort(key=lambda x: x[0])

    logger.info(f"Loaded set file: {set_file}")
    logger.info(f"Genes in set file: {len(gene_chrom):,}")
    return intervals_by_chrom, gene_chrom


def genes_for_position(chrom, pos, intervals_by_chrom):
    genes = []
    for start, end, gene in intervals_by_chrom.get(chrom, []):
        if pos < start:
            break
        if start <= pos <= end:
            genes.append(gene)
    return genes


def load_annotation(annotation_file, ann_col, keep_annotations, set_file=None):
    logger.info(f"Loading annotation: {annotation_file}")
    ann = pd.read_csv(annotation_file, sep="\t", dtype=str)
    ann = normalize_variant_df(ann)

    required = {"CHROM", "POS", "REF", "ALT", ann_col}
    missing = required - set(ann.columns)
    if missing:
        raise SystemExit(f"Annotation missing required columns: {missing}")

    if keep_annotations:
        keep_set = set(keep_annotations)

        def keep_ann(x):
            return any(a in keep_set for a in parse_ann_string(x))

        ann = ann[ann[ann_col].apply(keep_ann)].copy()
        logger.info(f"Annotation rows after consequence filter: {len(ann):,}")
    else:
        logger.info("No consequence filter requested.")

    ann = ann.dropna(subset=["CHROM", "POS", "REF", "ALT"]).copy()

    v2g = defaultdict(set)
    gene_chrom = {}

    if set_file is None:
        logger.info("No --set-file provided; using Gene_Name from annotation.")
        if "Gene_Name" not in ann.columns:
            raise SystemExit("Gene_Name column required when --set-file is not provided.")

        ann = ann.dropna(subset=["Gene_Name"]).copy()
        ann["Gene_Name"] = ann["Gene_Name"].astype(str).str.strip()

        for _, r in ann.iterrows():
            key = make_key(r)
            gene = r["Gene_Name"]
            v2g[key].add(gene)
            gene_chrom[gene] = str(r["CHROM"]).replace("chr", "")

    else:
        logger.info("Using --set-file gene boundaries for variant-to-gene assignment.")
        intervals_by_chrom, gene_chrom_from_set = load_set_file(set_file)

        for _, r in ann.iterrows():
            key = make_key(r)
            chrom, pos = key[0], key[1]
            genes = genes_for_position(chrom, pos, intervals_by_chrom)

            for gene in genes:
                v2g[key].add(gene)
                gene_chrom[gene] = gene_chrom_from_set.get(gene, chrom.replace("chr", ""))

    logger.info(f"Annotated variants assigned to genes: {len(v2g):,}")
    logger.info(f"Genes after assignment: {len(gene_chrom):,}")
    return v2g, gene_chrom


def read_header_cols(path):
    return pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).str.strip().tolist()


def choose_sample_columns(dos_cols, hap_cols, phenotype_ids):
    dos_samples = [c for c in dos_cols if c not in META_COLS]
    hap_samples = [c for c in hap_cols if c not in META_COLS]

    hap_set = set(hap_samples)
    phenotype_ids = set(map(str, phenotype_ids))

    samples = [s for s in dos_samples if s in hap_set and s in phenotype_ids]

    if len(samples) == 0:
        raise SystemExit("No overlapping samples among phenotype, dosage, and hapcount files.")

    return samples


def compute_total_keep(
    v2g,
    ancestry_names,
    dosage_files,
    hapcount_files,
    phenotype_ids,
    min_mac,
    maf_lo,
    maf_hi,
    chunksize,
):
    logger.info("Computing total MAF keep-set restricted to phenotype samples.")

    total_ac = defaultdict(float)
    total_an = defaultdict(float)
    annotated_keys = set(v2g.keys())

    for anc, dos_path, hap_path in zip(ancestry_names, dosage_files, hapcount_files):
        dos_cols = read_header_cols(dos_path)
        hap_cols = read_header_cols(hap_path)
        sample_ids = choose_sample_columns(dos_cols, hap_cols, phenotype_ids)
        usecols = META_COLS + sample_ids

        dos_iter = pd.read_csv(dos_path, sep="\t", dtype=str, chunksize=chunksize, usecols=usecols)
        hap_iter = pd.read_csv(hap_path, sep="\t", dtype=str, chunksize=chunksize, usecols=usecols)

        for dos_chunk, hap_chunk in zip(dos_iter, hap_iter):
            dos_chunk.columns = dos_chunk.columns.astype(str).str.strip()
            hap_chunk.columns = hap_chunk.columns.astype(str).str.strip()

            dos_chunk = normalize_variant_df(dos_chunk)
            hap_chunk = normalize_variant_df(hap_chunk)

            keys = [
                (row.CHROM, int(row.POS), row.REF, row.ALT)
                for row in dos_chunk[ID_COLS].itertuples(index=False)
            ]

            hit_idx = [i for i, k in enumerate(keys) if k in annotated_keys]
            if not hit_idx:
                continue

            dos_sub = dos_chunk.iloc[hit_idx]
            hap_sub = hap_chunk.iloc[hit_idx]
            sub_keys = [keys[i] for i in hit_idx]

            dos_mat = dos_sub[sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
            hap_mat = hap_sub[sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)

            ac = np.nansum(dos_mat, axis=1)
            an = np.nansum(hap_mat, axis=1)

            for k, a, n in zip(sub_keys, ac, an):
                total_ac[k] += float(a)
                total_an[k] += float(n)

    keep = set()
    for k in annotated_keys:
        ac = total_ac.get(k, 0.0)
        an = total_an.get(k, 0.0)
        maf = ac / an if an > 0 else 0.0

        if ac >= float(min_mac) and maf >= float(maf_lo) and maf <= float(maf_hi):
            keep.add(k)

    logger.info(f"Total-MAF keep-set size: {len(keep):,}")
    return keep


def stream_burdens_one_ancestry(
    anc,
    dos_path,
    hap_path,
    v2g,
    phenotype_ids,
    maf_scope,
    total_keep,
    min_mac,
    maf_lo,
    maf_hi,
    chunksize,
):
    logger.info(f"[{anc}] Starting streaming burdens")
    logger.info(f"[{anc}] dosage:   {dos_path}")
    logger.info(f"[{anc}] hapcount: {hap_path}")

    dos_cols = read_header_cols(dos_path)
    hap_cols = read_header_cols(hap_path)

    sample_ids = choose_sample_columns(dos_cols, hap_cols, phenotype_ids)
    sample_ids = np.array(sample_ids, dtype=str)

    logger.info(f"[{anc}] phenotype-restricted sample count before burden computation: {len(sample_ids):,}")

    usecols = META_COLS + sample_ids.tolist()

    dos_iter = pd.read_csv(dos_path, sep="\t", dtype=str, chunksize=chunksize, usecols=usecols)
    hap_iter = pd.read_csv(hap_path, sep="\t", dtype=str, chunksize=chunksize, usecols=usecols)

    annotated_keys = set(v2g.keys())

    burdens = defaultdict(lambda: np.zeros(len(sample_ids), dtype=float))
    gene_mac = defaultdict(float)
    gene_n_variants = defaultdict(int)

    n_chunks = 0
    n_annotated_seen = 0
    n_kept = 0

    for dos_chunk, hap_chunk in zip(dos_iter, hap_iter):
        n_chunks += 1

        dos_chunk.columns = dos_chunk.columns.astype(str).str.strip()
        hap_chunk.columns = hap_chunk.columns.astype(str).str.strip()

        dos_chunk = normalize_variant_df(dos_chunk)
        hap_chunk = normalize_variant_df(hap_chunk)

        dos_keys = [
            (row.CHROM, int(row.POS), row.REF, row.ALT)
            for row in dos_chunk[ID_COLS].itertuples(index=False)
        ]
        hap_keys = [
            (row.CHROM, int(row.POS), row.REF, row.ALT)
            for row in hap_chunk[ID_COLS].itertuples(index=False)
        ]

        if dos_keys != hap_keys:
            raise SystemExit(f"[{anc}] dosage/hapcount chunks are not row-aligned.")

        hit_idx = [i for i, k in enumerate(dos_keys) if k in annotated_keys]
        if not hit_idx:
            continue

        n_annotated_seen += len(hit_idx)

        dos_sub = dos_chunk.iloc[hit_idx]
        hap_sub = hap_chunk.iloc[hit_idx]
        sub_keys = [dos_keys[i] for i in hit_idx]

        dos_mat = dos_sub[sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
        hap_mat = hap_sub[sample_ids].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)

        ac = np.nansum(dos_mat, axis=1)

        if maf_scope == "none":
            keep = ac >= float(min_mac)

        elif maf_scope == "ancestry":
            an = np.nansum(hap_mat, axis=1)
            maf = np.divide(ac, an, out=np.zeros_like(ac, dtype=float), where=an > 0)
            keep = (ac >= float(min_mac)) & (maf >= float(maf_lo)) & (maf <= float(maf_hi))

        elif maf_scope == "total":
            if total_keep is None:
                raise SystemExit("maf_scope=total but total_keep is None.")
            keep = np.array([k in total_keep for k in sub_keys], dtype=bool)

        else:
            raise SystemExit(f"Unknown maf_scope: {maf_scope}")

        if not np.any(keep):
            continue

        kept_keys = [k for k, ok in zip(sub_keys, keep) if ok]
        kept_dos = dos_mat[keep, :]
        kept_ac = ac[keep]

        n_kept += len(kept_keys)

        for row_idx, key in enumerate(kept_keys):
            variant_dos = kept_dos[row_idx, :]
            variant_ac = float(kept_ac[row_idx])

            for gene in v2g[key]:
                burdens[gene] += variant_dos
                gene_mac[gene] += variant_ac
                gene_n_variants[gene] += 1

        if n_chunks % 25 == 0:
            logger.info(
                f"[{anc}] processed chunks={n_chunks:,}, "
                f"annotated_seen={n_annotated_seen:,}, kept={n_kept:,}, "
                f"genes_with_burden={len(burdens):,}"
            )

    gene_n_carriers = {g: int(np.sum(burdens[g] > 0)) for g in burdens.keys()}

    logger.info(
        f"[{anc}] Finished streaming: chunks={n_chunks:,}, "
        f"annotated_seen={n_annotated_seen:,}, kept={n_kept:,}, "
        f"genes_with_burden={len(burdens):,}"
    )

    return sample_ids, dict(burdens), dict(gene_mac), gene_n_carriers, dict(gene_n_variants)


def is_binary_y(y):
    vals = sorted(pd.Series(y).dropna().unique())
    return set(vals).issubset({0, 1}) and len(vals) <= 2


def run_regression(y, X):
    X = sm.add_constant(X, has_constant="add")

    if is_binary_y(y):
        model = sm.Logit(y, X, missing="drop")
        fit = model.fit(disp=False, maxiter=100)
    else:
        model = sm.OLS(y, X, missing="drop")
        fit = model.fit()

    return fit


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--annotation-file", required=True)
    p.add_argument("--ann-col", required=True)
    p.add_argument("--set-file", default=None, help="Optional RVTESTS-style refFlat.set gene boundary file.")

    p.add_argument("--ancestry-names", nargs="+", required=True)
    p.add_argument("--dosage-files", nargs="+", required=True)
    p.add_argument("--hapcount-files", nargs="+", required=True)

    p.add_argument("--phenotype-file", required=True)
    p.add_argument("--out-tsv", required=True)

    p.add_argument("--keep-annotations", nargs="*", default=[])

    p.add_argument("--min-mac", type=float, default=1)
    p.add_argument("--maf-scope", choices=["none", "ancestry", "total"], default="none")
    p.add_argument("--maf-lo", type=float, default=0.0)
    p.add_argument("--maf-hi", type=float, default=0.01)

    p.add_argument("--covariates", nargs="*", default=[])
    p.add_argument("--chunksize", type=int, default=256)

    return p.parse_args()


def main():
    args = parse_args()

    if not (
        len(args.ancestry_names)
        == len(args.dosage_files)
        == len(args.hapcount_files)
    ):
        raise SystemExit("ancestry-names, dosage-files, and hapcount-files must have same length.")

    logger.info(f"Loading phenotype: {args.phenotype_file}")
    pheno = pd.read_csv(args.phenotype_file, sep="\t", dtype={"IID": str})

    if "y" not in pheno.columns:
        raise ValueError(f"'y' not found. Available columns: {list(pheno.columns)}")

    if "IID" not in pheno.columns or "y" not in pheno.columns:
        raise SystemExit("Phenotype file must contain IID and y.")

    pheno["IID"] = pheno["IID"].astype(str).str.strip()
    pheno["y"] = pd.to_numeric(pheno["y"], errors="coerce")

    for c in args.covariates:
        if c not in pheno.columns:
            raise SystemExit(f"Covariate {c} not found in phenotype file.")
        pheno[c] = pd.to_numeric(pheno[c], errors="coerce")

    model_cols = ["IID", "y"] + args.covariates
    pheno_model = pheno[model_cols].dropna().copy()
    phenotype_ids = set(pheno_model["IID"].astype(str))

    logger.info(f"Phenotype rows usable for model: {len(pheno_model):,}")

    v2g, gene_chrom = load_annotation(
        annotation_file=args.annotation_file,
        ann_col=args.ann_col,
        keep_annotations=args.keep_annotations,
        set_file=args.set_file,
    )

    total_keep = None
    if args.maf_scope == "total":
        total_keep = compute_total_keep(
            v2g=v2g,
            ancestry_names=args.ancestry_names,
            dosage_files=args.dosage_files,
            hapcount_files=args.hapcount_files,
            phenotype_ids=phenotype_ids,
            min_mac=args.min_mac,
            maf_lo=args.maf_lo,
            maf_hi=args.maf_hi,
            chunksize=args.chunksize,
        )

    all_burdens = {}
    all_gene_mac = {}
    all_n_carriers = {}
    all_n_variants = {}

    sample_ids_master = None

    for anc, dos_path, hap_path in zip(args.ancestry_names, args.dosage_files, args.hapcount_files):
        sample_ids, burdens, gene_mac, gene_n_carriers, gene_n_variants = stream_burdens_one_ancestry(
            anc=anc,
            dos_path=dos_path,
            hap_path=hap_path,
            v2g=v2g,
            phenotype_ids=phenotype_ids,
            maf_scope=args.maf_scope,
            total_keep=total_keep,
            min_mac=args.min_mac,
            maf_lo=args.maf_lo,
            maf_hi=args.maf_hi,
            chunksize=args.chunksize,
        )

        if sample_ids_master is None:
            sample_ids_master = sample_ids
        else:
            if not np.array_equal(sample_ids_master, sample_ids):
                raise SystemExit("Sample order differs across ancestries after phenotype restriction.")

        all_burdens[anc] = burdens
        all_gene_mac[anc] = gene_mac
        all_n_carriers[anc] = gene_n_carriers
        all_n_variants[anc] = gene_n_variants

    sample_ids = sample_ids_master
    logger.info(f"Final phenotype-restricted sample count used in burdens/models: {len(sample_ids):,}")

    ph = pheno_model.set_index("IID").reindex(sample_ids)

    valid_mask = ph["y"].notna().to_numpy()
    for c in args.covariates:
        valid_mask &= ph[c].notna().to_numpy()

    logger.info(f"Samples with non-missing y/covariates: {int(valid_mask.sum()):,}")

    y = ph.loc[valid_mask, "y"].astype(float)

    covar_df = pd.DataFrame(index=ph.index)
    for c in args.covariates:
        covar_df[c] = ph[c].astype(float)
    covar_df = covar_df.loc[valid_mask]

    genes = sorted(set().union(*[set(b.keys()) for b in all_burdens.values()]))
    logger.info(f"Genes with any burden: {len(genes):,}")

    results = []

    for gene in genes:
        X = covar_df.copy()

        terms = []
        for anc in args.ancestry_names:
            term = f"burden_{anc}"
            bvec = all_burdens.get(anc, {}).get(gene, np.zeros(len(sample_ids), dtype=float))
            X[term] = bvec[valid_mask]
            terms.append((anc, term))

        if all(X[term].sum() == 0 for _, term in terms):
            continue

        try:
            fit = run_regression(y, X)
        except Exception as e:
            logger.warning(f"Regression failed for {gene}: {e}")
            continue

        for anc, term in terms:
            if term not in fit.params.index:
                continue

            pval = float(fit.pvalues.get(term, np.nan))
            estimate = float(fit.params.get(term, np.nan))

            results.append({
                "chrom": gene_chrom.get(gene, np.nan),
                "gene": gene,
                "term": term,
                "estimate": estimate,
                "pval": pval,
                "mac": float(all_gene_mac.get(anc, {}).get(gene, 0.0)),
                "n_carriers": int(all_n_carriers.get(anc, {}).get(gene, 0)),
                "n_variants": int(all_n_variants.get(anc, {}).get(gene, 0)),
            })

    res = pd.DataFrame(results)

    if res.empty:
        logger.warning("No results generated.")
        res.to_csv(args.out_tsv, sep="\t", index=False)
        return

    m_by_term = {}
    for anc in args.ancestry_names:
        term = f"burden_{anc}"
        m_by_term[term] = sum(
            int(all_n_carriers.get(anc, {}).get(g, 0)) > 0
            for g in genes
        )

    res["m_genes_tested"] = res["term"].map(m_by_term).astype(int)
    res["neglog10p"] = -np.log10(res["pval"].replace(0, np.nextafter(0, 1)))

    cols = [
        "chrom", "gene", "term", "estimate", "pval", "neglog10p",
        "mac", "n_carriers", "n_variants", "m_genes_tested"
    ]

    res = res[cols].sort_values("pval", ascending=True)

    logger.info(f"Writing results: {args.out_tsv}")
    res.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
