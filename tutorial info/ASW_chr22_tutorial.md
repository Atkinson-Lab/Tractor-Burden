# ASW chr22 Tractor-Burden Tutorial

This tutorial demonstrates how to run Tractor-Burden on publicly available chromosome 22 ASW samples from the Tractor tutorial. The workflow uses local ancestry tracts inferred by RFMix and projects ancestry labels onto rare variants before performing ancestry-specific burden testing.

## Required Files

Download the following files:

```text
ASW.phased.vcf.gz
ASW.deconvoluted.msp.tsv
ASW.deconvoluted.rfmix.Q
chr22.refFlat.set
ASW.pheno.covars.tsv
tractor_burden.py
```

## Step 1: Create a Rare Variant VCF

Filter the phased VCF to variants with MAF ≤ 1%.

```bash
bcftools view \
    -i 'MAF<=0.01' \
    ASW.phased.vcf.gz \
    -Oz \
    -o ASW.rare.vcf.gz

tabix -p vcf ASW.rare.vcf.gz
```

This file will contain only rare variants and serves as the input for ancestry projection.

---

## Step 2: Extract Ancestry-Specific Rare Variant Dosages

Use Tractor's `extract_tracts.py` script to assign local ancestry labels to rare variants using the RFMix MSP file.

```bash
python extract_tracts.py \
    --vcf ASW.rare.vcf.gz \
    --msp ASW.deconvoluted.msp.tsv \
    --num-ancs 2 \
    --output-dir .
```

The MSP file was generated using common variants during local ancestry inference. Rare variants are assigned ancestry labels based on the ancestry tract in which they reside.

---

## Step 3: Generate Tractor-Burden Input Files

Convert the Tractor outputs into Tractor-Burden dosage and hapcount files.

Final files used by Tractor-Burden:

```text
ASW.rare.anc0.dosage.txt
ASW.rare.anc1.dosage.txt

ASW.rare.anc0.hapcount.txt
ASW.rare.anc1.hapcount.txt
```

Each file contains:

```text
CHROM  POS  ID  REF  ALT  sample1  sample2  sample3 ...
```

where dosage files contain ancestry-specific alternate allele dosages and hapcount files contain ancestry-specific haplotype counts.

---

## Step 4: Phenotype and Covariate File

The tutorial phenotype file is:

```text
ASW.pheno.covars.tsv
```

with the following columns:

| Column | Description                                 |
| ------ | ------------------------------------------- |
| IID    | Sample identifier                           |
| y      | Simulated quantitative phenotype            |
| AFR    | Global African ancestry proportion          |
| EUR    | Global European ancestry proportion         |
| sex    | Binary sex covariate (0 = female, 1 = male) |

The phenotype is simulated for demonstration purposes and is not intended to represent a real biological trait.

---

## Step 5: Run Tractor-Burden

```bash
python tractor_burden.py \
    --set-file chr22.refFlat.set \
    --ancestry-names AFR EUR \
    --dosage-files \
        ASW.rare.anc0.dosage.txt \
        ASW.rare.anc1.dosage.txt \
    --hapcount-files \
        ASW.rare.anc0.hapcount.txt \
        ASW.rare.anc1.hapcount.txt \
    --phenotype-file ASW.pheno_covars.tsv \
    --covariates AFR sex \
    --maf-scope none \
    --min-mac 1 \
    --out-tsv ASW.results.tsv
```

---

## Step 6: Examine Results

View the top associations:

```bash
head ASW.results.tsv
```

The output contains:

| Column         | Description                                       |
| -------------- | --------------------------------------------------|
| chrom          | Chromosome                                        |
| gene           | Gene or region tested                             |
| term           | Ancestry-specific burden term                     |
| estimate       | Effect estimate                                   |
| pval           | Association p-value                               |
| neglog10p      | −log10(p-value)                                   |
| mac            | Minor allele count                                |
| n_carriers     | Number of carriers                                |
| n_variants     | Number of variants contributing to the burden     |
| m_genes_tested | Number of genes tested for that ancestry per chr  |

---

## Step 7: Interpreting Burden Results

A highly significant burden association is not necessarily driven by many variants.

For example, some genes may exhibit extremely strong associations while containing only a single contributing variant (`n_variants = 1`). In such cases, the signal may reflect a single-variant association rather than a classical multi-variant burden effect.

When prioritizing associations, users should evaluate:

* P-value
* Effect size
* Minor allele count (MAC)
* Number of carriers
* Number of contributing variants
* Biological plausibility

This illustrates an important downstream analysis step: burden test results should always be interpreted alongside burden composition metrics.

---

## Expected Output

This tutorial demonstrates the complete Tractor-Burden workflow:

```text
ASW.phased.vcf.gz
        ↓
bcftools MAF filter for rare variants
        ↓
ASW.rare.vcf.gz
        ↓
extract_tracts.py
        ↓
ASW.rare.anc0.dosage.txt
ASW.rare.anc1.dosage.txt
ASW.rare.anc0.hapcount.txt
ASW.rare.anc1.hapcount.txt
        ↓
Tractor-Burden
```

The resulting analysis produces several ancestry-specific burden associations and illustrates how local ancestry information can be incorporated into rare variant association testing.
