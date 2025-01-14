import pandas as pd
import numpy as np
import category_encoders as ce
import logging
import os
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt

# Create a logger
logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Directories
processed_dir = "data"
sample_sheets_dir = "sample_sheets"

df = pd.read_csv(f"{sample_sheets_dir}/gdc_clinical.tsv", sep="\t")

# Select only some fields
selected_cols = [
    "case_submitter_id",
    "age_at_index",
    "days_to_death",
    "days_to_last_follow_up",
    "morphology",
    "ethnicity",
    "gender",
    "race",
    "vital_status",
    # "year_of_birth",
    # "year_of_death",
    "ajcc_pathologic_m",
    "ajcc_pathologic_n",
    "ajcc_pathologic_stage",
    "primary_diagnosis",
    "treatment_or_therapy",
    "treatment_type",
]

df = df[selected_cols]

# Replace '--' with NaN
df = df.replace("'--", np.nan)

logger.info(df.shape)

# Add new column for Pharmaceutical Therapy, NOS
# Add new column for Radiation Therapy, NOS
df['pharmaceutical_treatment'] = ((df['treatment_type'] == 'Pharmaceutical Therapy, NOS') & (df["treatment_or_therapy"] == "yes")).astype(int)
df['radiation_treatment'] = ((df['treatment_type'] == 'Radiation Therapy, NOS') & (df["treatment_or_therapy"] == "yes")).astype(int)

# Group by 'case_submitter_id' and aggregate
df = df.groupby('case_submitter_id').agg({
    "age_at_index": "first",
    "days_to_death": "first",
    "days_to_last_follow_up": "first",
    "ethnicity": "first",
    "gender": "first",
    "race": "first",
    "age_at_index": "first",
    "morphology": "first",
    "vital_status": "first",
    "ajcc_pathologic_m": "first",
    "ajcc_pathologic_n": "first",
    "ajcc_pathologic_stage": "first",
    "primary_diagnosis": "first",
    "pharmaceutical_treatment": "max",             # Update to 1 if any row has 1
    "radiation_treatment": "max"                   # Update to 1 if any row has 1
}).reset_index()

# Update 'treatment_type' based on 'treatment_or_therapy'
df['treatment_type'] = np.where(df['pharmaceutical_treatment'] == 1, 'Pharmaceutical Therapy, NOS', 'None')
df['treatment_type'] = np.where(df['radiation_treatment'] == 1, 'Radiation Therapy, NOS', df['treatment_type'])

# Remove rows where 'vital_status' is NaN
df = df.dropna(subset=['vital_status'])

# Remove rows where 'days_to_death' or 'days_to_last_follow_up' is NaN
df = df[(df["days_to_death"].notna()) | (df["days_to_last_follow_up"].notna())]

# Convert vital_status to a binary event indicator (1 = Dead, 0 = Alive)
df['event'] = (df['vital_status'] == 'Dead').astype(int)

# Get 'days_to_event'
df['days_to_event'] = np.where((df['days_to_death'].notna()) & (df['event'] == 1), df['days_to_death'], df['days_to_last_follow_up']).astype(float)

df = df[df['days_to_event'] > 0]

# Rename some columns
df = df.rename(columns={'case_submitter_id': 'case_id'})

# Drop some columns
df = df.drop(columns=['pharmaceutical_treatment', 'radiation_treatment'])

logger.info(f"Total samples: {df.shape[0]}")
logger.info(f"Total clincal features: { df.shape[1]}")

# Export data
df.to_csv(f"{processed_dir}/clinical.tsv", sep="\t")
