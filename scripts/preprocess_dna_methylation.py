import pandas as pd
import numpy as np
import category_encoders as ce
import argparse
import logging
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# Create a logger
logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Input
parser = argparse.ArgumentParser(description="Preprocess DNA Methylation")
parser.add_argument('--data', type=str, required=True, help='Data dir')
args = parser.parse_args()

# Directories
processed_dir = "data"
sample_sheets_dir = "sample_sheets"
data_dir = args.data

# Load sample sheet
df_sample_sheet = pd.read_csv(f"{sample_sheets_dir}/gdc_dna_methylation_sample_sheet.tsv", sep="\t")

# Check for duplication cases
value_counts = df_sample_sheet['Case ID'].value_counts()
logger.info("Duplicated rows:")
logger.info(df_sample_sheet[df_sample_sheet['Case ID'].isin(value_counts[value_counts > 1].index)])

# Keeping only Primary Tumor
# df_sample_sheet = df_sample_sheet[df_sample_sheet['Sample Type'] == 'Primary Tumor']

# Separate tumor and normal samples
df_sample_sheet_tumor = df_sample_sheet[df_sample_sheet["Sample Type"] == "Primary Tumor"].reset_index(drop=True)
df_sample_sheet_normal = df_sample_sheet[df_sample_sheet["Sample Type"] == "Solid Tissue Normal"].reset_index(drop=True)

# Keep the first Case ID only
# df_sample_sheet = df_sample_sheet.drop_duplicates(subset="Case ID", keep='first')
# df_sample_sheet = df_sample_sheet.reset_index(drop=True)

logger.info(f"df_sample_sheet_tumor.shape: {df_sample_sheet_tumor.shape}")
logger.info(f"df_sample_sheet_normal.shape: {df_sample_sheet_normal.shape}")
# logger.info(df_sample_sheet.head(5))

def process_samples(df_sheet):
    # Loop through all rows
    case_data = []
    column_names = ["id", "value"]

    logger.info(f"df_sheet.shape:{df_sheet.shape}")

    for index, row in df_sheet.iterrows():
        file_path = f"{data_dir}/{row['File ID']}/{row['File Name']}"
        case_id = row["Case ID"]
        # Check if file exist
        if os.path.exists(file_path):
            # logger.info(f"Processing case: {row['Case ID']} at {dna_methylation_data_file_path}")
            df_dna_methylation_data = pd.read_csv(file_path, header=None, names=column_names, sep="\t")
            
            # Include only cg probes
            df_dna_methylation_data = df_dna_methylation_data[df_dna_methylation_data['id'].astype(str).str.startswith('cg')]

            # Exclude NaN value
            df_dna_methylation_data = df_dna_methylation_data[df_dna_methylation_data['value'].notna()]
        
            logger.info(f"At: {index}. Shape: {df_dna_methylation_data.shape}")

            df_case_data = df_dna_methylation_data.set_index('id').T
            df_case_data['case_id'] = case_id
            df_case_data = df_case_data.set_index('case_id')

            case_data.append(df_case_data)
        else:
            logger.error(f"File {file_path} not exists!")
    
    df = pd.concat(case_data, axis=0)
    df.columns.name = None

    logger.info("Done concat")

    # Calculate the percentage of missing values per column
    missing_percentage = df.isnull().mean() * 100

    # Summarize the number of columns in different ranges of missing values
    summary = {
        '0-10%': (missing_percentage <= 10).sum(),
        '10-30%': ((missing_percentage > 10) & (missing_percentage <= 30)).sum(),
        '30-50%': ((missing_percentage > 30) & (missing_percentage <= 50)).sum(),
        '50-100%': (missing_percentage > 50).sum()
    }
    logger.info("Summary of missing value percentages:")
    logger.info(summary)

    # Define a threshold for missing values
    threshold = 0.2  # 20% threshold

    # Drop columns with more than 20% missing values
    # df_processed = df.loc[:, missing_percentage <= (threshold * 100)]
    # logger.info(f"Drop columns with more than 20% missing values: {df_processed.shape}")

    # Drop columns with low variance
    # df_processed = df_processed.loc[df_processed.var(axis=1) > 0.01]
    # logger.info(f"Drop columns with low variance: {df_processed.shape}")


    return df

# Process tumor samples
df_tumor = process_samples(df_sample_sheet_tumor)

# Process normal samples
df_normal = process_samples(df_sample_sheet_normal)

# Save data
df_tumor.to_csv(f"{processed_dir}/dna_methylation_tumor.tsv", sep="\t")
df_normal.to_csv(f"{processed_dir}/dna_methylation_normal.tsv", sep="\t")

logger.info("Data saved")


def perform_t_test(df_tumor, df_normal, load=False):
    # Exclude 'case_id' from the CpG columns
    tumor_cpgs = df_tumor.loc[:, df_tumor.columns != 'case_id']
    normal_cpgs = df_normal.loc[:, df_normal.columns != 'case_id']
    
    # Find common CpG sites
    cpgs = tumor_cpgs.columns.intersection(normal_cpgs.columns)
    tumor_cpgs_only = tumor_cpgs.columns.difference(normal_cpgs.columns).tolist()

    logger.info(len(cpgs))

    # Calculate means for tumor and normal in a vectorized manner
    tumor_means = tumor_cpgs[cpgs].mean(axis=0)
    normal_means = normal_cpgs[cpgs].mean(axis=0)

    if not load:
        # Vectorized log-fold change calculation
        log_fold_changes = np.log2(tumor_means + 1e-5) - np.log2(normal_means + 1e-5)

        logger.info(f"log_fold_changes: {log_fold_changes}")

        # Perform t-test in a vectorized manner
        _, p_values = ttest_ind(tumor_cpgs[cpgs], normal_cpgs[cpgs], nan_policy='omit', axis=0)
        logger.info(f"p_values: {p_values}")

        # Create results DataFrame
        results = pd.DataFrame({
            "cpg_site": cpgs,
            "log_fold_change": log_fold_changes,
            "p_value": p_values
        })

        # Adjust p-values using Benjamini-Hochberg method
        results = results.sort_values(by="p_value").reset_index(drop=True)
        results["adj_p_value"] = results["p_value"] * len(results) / (np.arange(1, len(results) + 1))
        results["adj_p_value"] = results["adj_p_value"].clip(upper=1.0)  # Ensure p-values are capped at 1.0

        # Display top significant results
        significant_probes = results[(results["adj_p_value"] < 0.05) & (abs(results["log_fold_change"]) > 0.5)]

        significant_probes.to_csv(f"{processed_dir}/significant_probes.tsv", sep="\t")
        results.to_csv(f"{processed_dir}/dna_methylation_probes.tsv", sep="\t")

        scpgs = significant_probes["cpg_site"].to_list()

        logger.info(significant_probes.head())
        logger.info(f"significant_probes.shape: {significant_probes.shape}")
    else:
        significant_probes = pd.read_csv(f"{processed_dir}/significant_probes.tsv", sep="\t")
        scpgs = significant_probes["cpg_site"].to_list()

    return scpgs, tumor_cpgs_only

logger.info('Perform t-test')
scpgs, tumor_cpgs_only = perform_t_test(df_tumor, df_normal, False)

all_cpgs = list(set(scpgs + tumor_cpgs_only))

logger.info(f"All cpgs: {len(all_cpgs)}")

final_columns = all_cpgs

df_tumor_final = df_tumor[final_columns]
df_tumor_final.to_csv(f"{processed_dir}/dna_methylation.tsv", sep="\t")

logger.info(f"df_tumor_final.shape: {df_tumor_final.shape}")

# Impute missing values with the mean of each probe
# df_processed = df_processed.fillna(df_processed.mean())

# Add prefix for every field
# df_processed = df_processed.add_prefix("DNA_meth_")

# df_processed.to_csv(f"{processed_dir}/dna_methylation.tsv", sep="\t")