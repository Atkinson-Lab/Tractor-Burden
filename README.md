# Tractor-Burden

Tractor-Burden is an ancestry-aware rare variant burden testing framework for admixed populations. Similar to the intuition behind Tractor, Tractor-Burden incorporates local ancestry information into rare variant association testing by aggregating ancestry-specific rare variant dosages within user-defined genomic regions and performing ancestry-specific burden tests.

This framework supports:

- Multi-ancestry burden testing
- Ancestry-specific burden models
- Binary or quantitative phenotypes
- Annotation-based variant filtering
- Region-based aggregation using either:
  - annotation `Gene_Name`
  - RVTESTS-style `.set` files
- User-defined genomic intervals 

---

## Overview

For each ancestry, Tractor-Burden:

1. Loads ancestry-specific dosage and hapcount files
2. Filters variants based on annotation class and MAF thresholds
3. Aggregates rare variant dosages within genomic regions
4. Fits association models:
   - Logistic regression for binary traits
   - Linear regression for quantitative traits
5. Returns ancestry-specific burden statistics

---

# Required Inputs

## 1. Annotation File

Tab-delimited annotation file containing variant information.

### Required Columns

| Column | Description |
|----------|-------------|
| CHROM | Chromosome |
| POS | Position |
| REF | Reference allele |
| ALT | Alternate allele |
| Annotation column | Functional annotation specified with `--ann-col` |

### Optional Columns

| Column | Description |
|----------|-------------|
| Gene_Name | Gene assignment used for gene-based aggregation |

### Example

```text
CHROM   POS      REF ALT Gene_Name Annotation
1       12345    A   G   LDLR      missense_variant
1       12500    T   C   LDLR      frameshift_variant
1       13000    G   A   APOB      splice_acceptor_variant
```

Example annotations:

```text
missense_variant
frameshift_variant
splice_acceptor_variant
splice_donor_variant
stop_gained
```

---

## 2. Dosage Files

Ancestry-specific dosage files generated from Tractor.

Example:

```text
chr1.anc0.dosage.gz
chr1.anc1.dosage.gz
```

Each dosage file contains ancestry-specific minor allele dosages.

---

## 3. Hapcount Files

Ancestry-specific haplotype count files generated from Tractor.

Example:

```text
chr1.anc0.hapcount.gz
chr1.anc1.hapcount.gz
```

These are used to calculate ancestry-specific allele frequencies and apply rare variant filters.

---

## 4. Phenotype File

Tab-delimited phenotype file.

### Required Columns

```text
IID
PHENO
```

### Example

```text
IID     PHENO
ID001   1
ID002   0
ID003   1
```

For quantitative traits:

```text
IID     PHENO
ID001   112.3
ID002   97.4
ID003   105.2
```

---

## 5. Covariate File

Tab-delimited covariate file.

Example:

```text
IID age sex AFR_prop
ID001 52 1 0.80
ID002 44 0 0.65
ID003 61 1 0.92
```

Typical covariates include:

- age
- sex
- global ancestry proportions
- principal components

---

# Aggregation Methods

## Option 1: Gene-Based Aggregation

Use the `Gene_Name` column in the annotation file.

Example:

```bash
--gene-col Gene_Name
```

Variants assigned to the same gene will be collapsed into a single burden score.

---

## Option 2: RVTESTS-Style Set File

Use a custom region definition file.

Example:

```text
LDLR chr19:11000000-11200000
APOB chr2:21000000-21200000
```

Run with:

```bash
--set-file gene_boundaries.set
```

---

# Flexible Region-Based Aggregation

Although Tractor-Burden is commonly used for gene-based rare variant aggregation, the framework can aggregate variants across **any user-defined genomic regions**.

Examples include:

- genes
- enhancers
- promoters
- regulatory elements
- chromatin interaction regions
- sliding windows
- pathway intervals
- custom genomic regions

Example:

```text
Region1 chr1:100000-150000
EnhancerA chr2:250000-275000
CustomWindow chr5:1000000-1100000
```

This enables Tractor-Burden analyses beyond coding variation into noncoding and regulatory regions.

---

# Running Tractor-Burden


Tractor-Burden automatically detects whether the phenotype column `y` is binary or quantitative. Users do **not** need to specify the trait type.

- If `y` contains only `0/1` values, Tractor-Burden runs logistic regression.
- Otherwise, Tractor-Burden runs linear regression.

## Example Run

```bash
python tractor_burden_final.py \
  --annotation-file /path/to/annotations.tsv \
  --ann-col consolidated_annotation \
  --set-file /path/to/regions.refFlat.set \
  --ancestry-names EUR AFR \
  --dosage-files \
    /path/to/chr19.anc0.dosage.txt.gz \
    /path/to/chr19.anc1.dosage.txt.gz \
  --hapcount-files \
    /path/to/chr19.anc0.hapcount.txt.gz \
    /path/to/chr19.anc1.hapcount.txt.gz \
  --phenotype-file /path/to/phenotype.tsv \
  --out-tsv /path/to/tractor_burden_results.tsv \
  --keep-annotations \
    missense_variant \
    frameshift_variant \
    splice_acceptor_variant \
    splice_donor_variant \
    stop_gained \
    start_lost \
  --min-mac 1 \
  --maf-scope none \
  --covariates global_ancestry_AFR age sex \
  --chunksize 256
```

---

# Output

The output contains ancestry-specific burden association statistics.

Example:

```text
gene    ancestry    estimate    pval
LDLR    AFR         0.41        3.5e-07
LDLR    EUR         0.09        0.11
APOB    AFR         0.28        1.4e-04
APOB    EUR         0.05        0.42
```

Columns:

| Column | Description |
|----------|-------------|
| gene | Aggregated region |
| ancestry | Tested ancestry |
| estimate | Regression coefficient |
| pval | Association p-value |

Additional columns may include:

```text
neglog10p
mac
n_carriers
n_variants
m_genes_tested
```

---

# Recommended Variant Filters

### Protein-Truncating Variants (pLoF)

```text
frameshift_variant
stop_gained
splice_acceptor_variant
splice_donor_variant
start_lost
transcript_ablation
```

### Damaging Missense Variants

```text
missense_variant
```

(optionally filtered by REVEL, CADD, or other pathogenicity scores)

### Combined Coding

```text
missense_variant
frameshift_variant
stop_gained
splice_acceptor_variant
splice_donor_variant
```

---

# Example Workflow

![Tractor-Burden workflow](figures/tractor_burden_workflow.png)

---

# Additional Resources

For detailed explanations of phasing, local ancestry painting, and extracting tracts, refer to:

**[Tractor Tutorial](https://github.com/Atkinson-Lab/Tractor-tutorial/tree/main)**

---

# Citation

If you use this pipeline in your research, please cite:

**XYZ (bioRxiv)**


Please direct questions to: pragati.kore@bcm.edu or elizabeth.atkinson@bcm.edu.
