# tcga

### CNV

#### Download data

```bash
mkdir -p /mnt/d/Documents/Data/TCGA/CNV

gdc-client download -d /mnt/d/Documents/Data/TCGA/CNV -m manifests/cnv_manifest.txt
```

#### Notes

- Should use the column `copy_number`
- Can do the circle chart for ploting missing values

### DNA Methylation

#### Download data

```bash
mkdir -p /mnt/d/Documents/Data/TCGA/DNA_Methylation

gdc-client download -d /mnt/d/Documents/Data/TCGA/DNA_Methylation -m manifests/dna_methylation_manifest.txt
```

### Gene expression

#### Download data

```bash
mkdir -p /mnt/d/Documents/Data/TCGA/gene_expression

gdc-client download -d /mnt/d/Documents/Data/TCGA/gene_expression -m manifests/gene_expression_quantification_manifest.txt
```

#### Notes

- Apply `z-score` normalization for `tpm_unstranded` for every samples.

### miRNA

#### Download data

```bash
mkdir -p /mnt/d/Documents/Data/TCGA/miRNA

gdc-client download -d /mnt/d/Documents/Data/TCGA/miRNA -m manifests/miRNA_manifest.txt
```

#### Notes

- Using `reads_per_million_miRNA_mapped`
- Should only get `cross-mapped` with value is `N` because `miRNAs` marked as `Y` can introduce noise (they might not be uniquely mapped to a single `miRNA`).



### Notes

- Using chi-square to reduce features
- Calculate risk scores to determine the features selection result