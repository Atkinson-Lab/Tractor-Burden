# Tractor-Burden Tutorial

The example dataset used in this tutorial consists of chromosome 22 data from 61 African American (ASW) individuals from the 1000 Genomes Project. These individuals are two-way admixed with African (AFR) and European (EUR) ancestral components. To demonstrate the Tractor-Burden workflow, we generated a simulated quantitative phenotype and accompanying covariates for use in ancestry-specific rare variant burden testing.

The phased genotype data and local ancestry inference (LAI) outputs used in this tutorial were obtained from the public [Tractor tutorial](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Rfmix.md). Local ancestry tracts were inferred using RFMix and are provided as MSP files, which are subsequently used to assign ancestry labels to rare variants prior to burden testing.

  
## Required Files

Download and unzip this example dataset to follow along the tutorial.

```text
ASW.phased.vcf.gz
ASW.deconvoluted.msp.tsv
ASW.deconvoluted.rfmix.Q
chr22.refFlat.set
ASW.pheno.covars.tsv
tractor_burden.py
```

## Step 1: Create a Rare Variant VCF

Index the phased VCF, then filter to variants with MAF ≤ 1%, and then index the resulting rare variant VCF.

```bash
tabix -p vcf ASW.phased.vcf.gz

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
> **Note:** This script requires Python 3. If your environment defaults to Python 2, you may encounter syntax errors. Verify your Python version before proceeding:

```bash
python3 --version
```

Then run:

```bash
python3 extract_tracts.py \
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
## Software Requirements

Before running Tractor-Burden, we recommend creating a dedicated Conda environment with the required Python packages. This helps avoid `ModuleNotFoundError` issues on shared computing clusters where packages such as `numpy`, `pandas`, `scipy`, or `statsmodels` may not be installed by default.

```bash
conda create -n tractor-burden python=3.12 numpy pandas scipy statsmodels -y
conda activate tractor-burden
```

If you are working on an HPC or shared server, you may need to load Conda or Anaconda first. For example:

```bash
module load anaconda3
conda activate tractor-burden
```

You can confirm that the required packages are available by running:

```bash
python -c "import numpy, pandas, scipy, statsmodels; print('All required packages loaded successfully')"
```

After activating the `tractor-burden` Conda environment:
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

