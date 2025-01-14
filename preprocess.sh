#!/bin/bash

docker run -itv .:/workspace namxle/py-annotate:22.04 bash

pip install scikit-learn category_encoders lifelines

# Server
cd /data/GL/master/TCGA/tcga/

# Preprocess clinical data
python3 scripts/preprocess_clinical.py

# Preprocess CNV data
# Local
python3 scripts/preprocess_cnv.py --data /mnt/d/Documents/Data/TCGA/CNV
# Server
python3 scripts/preprocess_cnv.py --data /data/GL/master/TCGA/data/TCGA/CNV

# Preprocess DNA Methylation data
# Local
python3 scripts/preprocess_dna_methylation.py --data /mnt/d/Documents/Data/TCGA/DNA_Methylation
# Server
python3 scripts/preprocess_dna_methylation.py --data /data/GL/master/TCGA/data/DNA_Methylation

# Preprocess Gene Expression data
# Local
python3 scripts/preprocess_gene_expression.py --data /mnt/d/Documents/Data/TCGA/gene_expression
# Server
python3 scripts/preprocess_gene_expression.py --data /data/GL/master/TCGA/data/gene_expression

# Preprocess miRNA data
# Local
python3 scripts/preprocess_miRNA.py --data /mnt/d/Documents/Data/TCGA/miRNA
# Server
python3 scripts/preprocess_miRNA.py --data /data/GL/master/TCGA/data/miRNA