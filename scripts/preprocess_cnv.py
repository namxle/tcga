import pandas as pd
import numpy as np
import category_encoders as ce
import logging
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
import argparse

# Create a logger
logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Input
parser = argparse.ArgumentParser(description="Preprocess DNA Methylation")
parser.add_argument('--data', type=str, required=True, help='Data dir')
args = parser.parse_args()

# Directories
processed_dir = f"data"
sample_sheets_dir = f"sample_sheets"
data_dir = args.data

# Load sample sheet
df_sample_sheet = pd.read_csv(f"{sample_sheets_dir}/gdc_cnv_sample_sheet.tsv", sep="\t")

# Load clinical data
df_clinical_data = pd.read_csv(f"{processed_dir}/clinical.tsv", sep="\t")
list_cases = df_clinical_data["case_id"].to_list()

# Remove duplication cases
df_sample_sheet = df_sample_sheet.drop_duplicates(subset="Case ID", keep='first')
df_sample_sheet = df_sample_sheet.reset_index(drop=True)
df_sample_sheet.shape
df_sample_sheet.head(5)

# Check for duplication cases
value_counts = df_sample_sheet['Case ID'].value_counts()
logger.info(df_sample_sheet[df_sample_sheet['Case ID'].isin(value_counts[value_counts > 1].index)][["File ID", "File Name", "Case ID"]])

# Get cnv data of case exist in clinical data
df_sample_sheet = df_sample_sheet[df_sample_sheet["Case ID"].isin(list_cases)]
logger.info(df_sample_sheet.shape)

# Loop through all rows
case_data = []
for index, row in df_sample_sheet.iterrows():
    cnv_data_file_path = f"{data_dir}/{row['File ID']}/{row['File Name']}"
    case_id = row["Case ID"]
    # Check if file exist
    if os.path.exists(cnv_data_file_path):
        # logger.info(f"Processing case: {row['Case ID']} at {cnv_data_file_path}")
        df_cnv_data = pd.read_csv(cnv_data_file_path, sep="\t")
        # print(df_cnv_data)
        df_cnv_data = df_cnv_data[df_cnv_data['copy_number'].notna()][['chromosome','gene_name', 'copy_number']]
        df_cnv_data = df_cnv_data.groupby('gene_name', as_index=False)['copy_number'].mean()
        
        # value_counts = df_cnv_data['gene_name'].value_counts()
        # print(df_cnv_data[df_cnv_data['gene_name'].isin(value_counts[value_counts > 1].index)][["gene_name", "copy_number"]])
        df_case_data = df_cnv_data.set_index('gene_name').T
        df_case_data['case_id'] = case_id
        df_case_data = df_case_data.set_index('case_id')
    
        case_data.append(df_case_data)
    else:
        logger.error(f"File {cnv_data_file_path} not exists!")

# Merge all cases
df = pd.concat(case_data, axis=0)
df.columns.name = None
df.to_csv(f"{processed_dir}/cnv_raw.tsv", sep="\t")

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

# Fill NA with CN = 2 which is normal for human
df_processed = df_processed.fillna(2)

# Add prefix for every field
df_processed = df_processed.add_prefix("CNV_gene_")

logger.info(df_processed.head())

logger.info(df_processed.shape)

# Export data
df_processed.to_csv(f"{processed_dir}/cnv.tsv", sep="\t")





