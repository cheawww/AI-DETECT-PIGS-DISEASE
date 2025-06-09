import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# File paths
fai_file = "GCF_000003025.6_Sscrofa11.1_genomic.fna.fai"
csv_file = "feature_weights_cfd_from_us (2).csv"

# Validate if files exist
if not os.path.exists(fai_file):
    print(f"Error: The file '{fai_file}' does not exist.")
    exit()

if not os.path.exists(csv_file):
    print(f"Error: The file '{csv_file}' does not exist.")
    exit()

# Validate if CSV file is empty
if os.stat(csv_file).st_size == 0:
    print(f"Error: The file '{csv_file}' is empty.")
    exit()

# Load the FAI file
try:
    fai_data = pd.read_csv(
        fai_file, sep="\t", header=None,
        names=["Chromosome", "Length", "Offset", "LineBases", "LineWidth"]
    )
except Exception as e:
    print(f"Error reading FAI file: {e}")
    exit()

# Load the filtered gRNA results CSV
try:
    csv_data = pd.read_csv(csv_file)
except Exception as e:
    print(f"Error reading CSV file: {e}")
    exit()

# Validate CSV file content
if csv_data.empty:
    print(f"Error: The file '{csv_file}' has no data.")
    exit()

# Create a figure with multiple subplots (tracks)
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15, 10), sharex=True)

# Track 1: Chromosome lengths as bars
axes[0].bar(fai_data["Chromosome"], fai_data["Length"], color="lightblue", edgecolor="black")
axes[0].set_ylabel("Chromosome Length")
axes[0].set_title("Chromosome Lengths")
axes[0].tick_params(axis='x', labelrotation=45)

# Ensure that the CSV file has sufficient columns
if csv_data.shape[1] < 8:
    print("Error: CSV file does not have enough columns to plot tracks.")
    exit()

# Track 2: A sample feature (e.g., column 5 in the CSV file)
axes[1].plot(csv_data.index, csv_data.iloc[:, 5], color="blue", label="Score 1")
axes[1].set_ylabel("Score 1")
axes[1].legend()

# Track 3: Another sample feature (e.g., column 6 in the CSV file)
axes[2].plot(csv_data.index, csv_data.iloc[:, 6], color="orange", label="Score 2")
axes[2].set_ylabel("Score 2")
axes[2].legend()

# Track 4: Additional features (e.g., column 7 in the CSV file)
axes[3].plot(csv_data.index, csv_data.iloc[:, 7], color="green", label="Normalized gRNA")
axes[3].set_ylabel("Normalized gRNA")
axes[3].legend()

# Adjust layout
plt.xlabel("Genomic Position (Base Pair Index)")
plt.xticks(rotation=45)
plt.tight_layout()

# Show the final plot
plt.show()
