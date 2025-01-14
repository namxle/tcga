import pandas as pd
import numpy as np
import category_encoders as ce
import logging
import argparse
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer
from sklearn.impute import KNNImputer
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import StandardScaler
from functools import reduce
from sklearn.preprocessing import OrdinalEncoder

from lifelines import CoxPHFitter
import pandas as pd


# Create a logger
logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

processed_dir = f"data"

logger.info("Loading data")

# Clinical data
df_clinical_data = pd.read_csv(f"{processed_dir}/clinical.tsv", sep="\t")
list_cases = df_clinical_data["case_id"].to_list()
df_clinical_data = df_clinical_data[[
    "case_id",
    "age_at_index",
    "ethnicity",
    "gender",
    "race",
    "morphology",
    "ajcc_pathologic_m",
    "ajcc_pathologic_n",
    "ajcc_pathologic_stage",
    "primary_diagnosis",
    "treatment_type",
    "event",
    "days_to_event"
]]

# logger.info(df_clinical_data["treatment_type"].unique())

# Standardize age
scaler = StandardScaler()
df_clinical_data['age_at_index'] = scaler.fit_transform(df_clinical_data[['age_at_index']])

# Label encode gender
encoder = LabelEncoder()
df_clinical_data['gender'] = encoder.fit_transform(df_clinical_data['gender'])
df_clinical_data['ethnicity'] = encoder.fit_transform(df_clinical_data['ethnicity'])
df_clinical_data['race'] = encoder.fit_transform(df_clinical_data['race'])
df_clinical_data['morphology'] = encoder.fit_transform(df_clinical_data['morphology'])
df_clinical_data['primary_diagnosis'] = encoder.fit_transform(df_clinical_data['primary_diagnosis'])
df_clinical_data['treatment_type'] = encoder.fit_transform(df_clinical_data['treatment_type'])


# Ordinal encoder
df_clinical_data['ajcc_pathologic_m'] = df_clinical_data['ajcc_pathologic_m'].fillna("Unknown")
order = [['Unknown', 'M0', 'MX', 'M1', 'M1a', 'M1b']]
encoder = OrdinalEncoder(categories=order)
df_clinical_data['ajcc_pathologic_m'] = encoder.fit_transform(df_clinical_data[['ajcc_pathologic_m']])

df_clinical_data['ajcc_pathologic_n'] = df_clinical_data['ajcc_pathologic_n'].fillna("Unknown")
order = [['Unknown', 'N0', 'NX', 'N1', 'N2', 'N3', ]]
encoder = OrdinalEncoder(categories=order)
df_clinical_data['ajcc_pathologic_n'] = encoder.fit_transform(df_clinical_data[['ajcc_pathologic_n']])

df_clinical_data['ajcc_pathologic_stage'] = df_clinical_data['ajcc_pathologic_stage'].fillna("Unknown")
order = [["Unknown", 'Stage I', 'Stage IA', 'Stage IB', 'Stage II', 'Stage IIA', 'Stage IIB', 'Stage IIIA', 'Stage IIIB', 'Stage IV']]
encoder = OrdinalEncoder(categories=order)
df_clinical_data['ajcc_pathologic_stage'] = encoder.fit_transform(df_clinical_data[['ajcc_pathologic_stage']])


# CNV data
df_cnv_data = pd.read_csv(f"{processed_dir}/cnv.tsv", sep="\t")
# df_cnv_merged = pd.merge(df_clinical_data, df_cnv_data, on='case_id', how='left')

# DNA methylation data
# df_dna_methylation_data = pd.read_csv(f"{processed_dir}/dna_methylation.tsv", sep="\t")

# Gene expression data
df_gene_expression_data = pd.read_csv(f"{processed_dir}/gene_expression.tsv", sep="\t")
# df_gene_expression_merged = pd.merge(df_clinical_data, df_gene_expression_data, on='case_id', how='left')

# miRNA data
df_miRNA_data =  pd.read_csv(f"{processed_dir}/miRNA.tsv", sep="\t")
# df_miRNA_merged = pd.merge(df_clinical_data, df_miRNA_data, on='case_id', how='left')

logger.info("Done load data")

logger.info("Merging data")

dfs = [df_cnv_data, df_gene_expression_data, df_miRNA_data]
omics_merged_df = reduce(lambda left, right: pd.merge(left, right, on='case_id', how='outer'), dfs)
merged_df = pd.merge(df_clinical_data, omics_merged_df, on='case_id', how='left')
merged_df = merged_df.fillna(0)

# Fill NaN values in selected columns with 0
cols = [col for col in merged_df.columns if not col.startswith('CNV')]
merged_df[cols] = merged_df[cols].fillna(0)

cols = [col for col in merged_df.columns if col.startswith('CNV')]
merged_df[cols] = merged_df[cols].fillna(2)

# Standard scale
scaler = StandardScaler()

# Standard scale gene expression
cols = merged_df.filter(like='Gene_exp_').columns
merged_df[cols] = scaler.fit_transform(merged_df[cols])

# Standard scale miRNA
cols = merged_df.filter(like='miRNA_').columns
merged_df[cols] = scaler.fit_transform(merged_df[cols])

merged_df = merged_df.reset_index(drop=True)
merged_df = merged_df.drop(columns=['case_id'])

logger.info("Done merging")

logger.info(merged_df.shape)
logger.info(merged_df.head())

# Initialize VarianceThreshold with a threshold (e.g., 0.1)
selector = VarianceThreshold(threshold=0.1)

# Fit and transform the DataFrame (removes low-variance columns)
merged_df_high_variance = pd.DataFrame(selector.fit_transform(merged_df), columns=merged_df.columns[selector.get_support()])

logger.info(f"merged_df.shape: {merged_df.shape}")
logger.info(f"merged_df_high_variance.shape: {merged_df_high_variance.shape}")

# Fit Cox model
cph = CoxPHFitter()
cph.fit(merged_df_high_variance, duration_col='days_to_event', event_col='event',)
cph.print_summary()

# Calculate risk scores
risk_scores = cph.predict_partial_hazard(merged_df)
print(risk_scores)