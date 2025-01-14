#!/bin/bash

# Preprocess clinical data
python3 scripts/preprocess_clinical.py

# Preprocess CNV data
python3 scripts/preprocess_cnv.py --data /mnt/d/Documents/Data/TCGA/CNV

# Preprocess DNA Methylation data
python3 scripts/preprocess_dna_methylation.py --data /mnt/d/Documents/Data/TCGA/DNA_Methylation

# Preprocess Gene Expression data

# Preprocess miRNA data