import pandas as pd
import numpy as np
import category_encoders as ce
import logging
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
import argparse
from sklearn.impute import KNNImputer
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import StandardScaler


# Create a logger
logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Input
parser = argparse.ArgumentParser(description="Preprocess miRNA")
parser.add_argument('--data', type=str, required=True, help='Data dir')
args = parser.parse_args()


# Directories
processed_dir = "data"
sample_sheets_dir = "sample_sheets"
data_dir = args.data

# Load sample sheet
df_sample_sheet = pd.read_csv(f"{sample_sheets_dir}/gdc_miRNA_sample_sheet.tsv", sep="\t")

# Check for duplication cases
value_counts = df_sample_sheet['Case ID'].value_counts()
logger.info(df_sample_sheet[df_sample_sheet['Case ID'].isin(value_counts[value_counts > 1].index)])

# Remove duplication cases
logger.info(df_sample_sheet.shape)

# Keeping only Primary Tumor
df_sample_sheet = df_sample_sheet[df_sample_sheet['Sample Type'] == 'Primary Tumor']

# Keep the first Case ID only
df_sample_sheet = df_sample_sheet.drop_duplicates(subset="Case ID", keep='first')
df_sample_sheet = df_sample_sheet.reset_index(drop=True)

logger.info(df_sample_sheet.shape)
logger.info(df_sample_sheet.head(5))


# Loop through all rows
case_data = []

logger.info(df_sample_sheet.shape)

for index, row in df_sample_sheet.iterrows():
    file_path = f"{data_dir}/{row['File ID']}/{row['File Name']}"
    case_id = row["Case ID"]
    logger.info(f"At: {index}.")
    # Check if file exist
    if os.path.exists(file_path):
        # logger.info(f"Processing case: {row['Case ID']} at {dna_methylation_data_file_path}")
        df_miRNA_data = pd.read_csv(file_path, sep="\t")
        
        # Include only protein coding genes
        df_miRNA_data = df_miRNA_data[df_miRNA_data['cross-mapped'] == 'N'][["miRNA_ID", "reads_per_million_miRNA_mapped"]]

        # Exclude NaN value
        df_miRNA_data = df_miRNA_data[df_miRNA_data['reads_per_million_miRNA_mapped'].notna()]
    
        # logger.info(f"At: {index}. Shape: {df_miRNA_data.shape}")

        df_case_data = df_miRNA_data.set_index('miRNA_ID').T
        df_case_data['case_id'] = case_id
        df_case_data = df_case_data.set_index('case_id')

        case_data.append(df_case_data)
    else:
        logger.error(f"File {file_path} not exists!")
logger.info("Done")


# Merge all cases
df = pd.concat(case_data, axis=0)
df.columns.name = None
df.to_csv(f"{processed_dir}/miRNA_raw.tsv", sep="\t")

logger.info(df.shape)

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
df_processed = df.loc[:, missing_percentage <= (threshold * 100)]
logger.info(df_processed.shape)

# Add prefix for every field
df_processed = df_processed.add_prefix("miRNA_")

# Normalization
# imputer = KNNImputer(n_neighbors=5)
# df_miRNA_imputed = pd.DataFrame(imputer.fit_transform(df_processed), columns=df_processed.columns, index=df_processed.index)

# scaler = StandardScaler()
# df_miRNA_scaled = pd.DataFrame(scaler.fit_transform(df_miRNA_imputed), columns=df_miRNA_imputed.columns, index=df_miRNA_imputed.index)

df_processed.to_csv(f"{processed_dir}/miRNA.tsv", sep="\t")

logger.info(df_processed.head())