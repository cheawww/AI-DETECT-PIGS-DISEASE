from BCBio import GFF
from Bio import SeqIO

fna_file = "GCF_000003025.6_Sscrofa11.1_genomic.fna"
gff_file = "GCF_000003025.6_Sscrofa11.1_genomic.gff"
target_region = (100000, 200000)  # Example range for verification

# Parse GFF
with open(gff_file) as gff_handle:
    for rec in GFF.parse(gff_handle):
        print(f"Chromosome: {rec.id}")
        for feature in rec.features:
            if feature.location.start >= target_region[0] and feature.location.end <= target_region[1]:
                print(f"Feature: {feature.type}, Start: {feature.location.start}, End: {feature.location.end}")

# Verify sequence in FNA
seq_record = SeqIO.read(fna_file, "fasta")
target_seq = seq_record.seq[target_region[0]:target_region[1]]
print(f"Target Sequence: {target_seq[:100]}...")  # Preview first 100 bases

import pandas as pd

gff_file = "GCF_000003025.6_Sscrofa11.1_genomic.gff"
gff_data = pd.read_csv(
    gff_file,
    sep="\t",
    comment="#",
    header=None,
    names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
)

features_of_interest = gff_data[gff_data["feature"].isin(["gene", "CDS"])]
print(features_of_interest.head())  # Preview extracted regions

gRNA_positions = [110468283, 3099452, 164620117, 24584244, 3347057]

for gRNA_pos in gRNA_positions:
    overlaps = features_of_interest[
        (features_of_interest["start"] <= gRNA_pos) & 
        (features_of_interest["end"] >= gRNA_pos)
    ]
    
    if not overlaps.empty:
        print(f"gRNA at position {gRNA_pos} overlaps with the following features:")
        print(overlaps[["seqname", "feature", "start", "end", "attributes"]])
    else:
        print(f"gRNA at position {gRNA_pos} does not overlap with any feature.")


    from Bio import SeqIO
    
fna_file = "GCF_000003025.6_Sscrofa11.1_genomic.fna"
record = SeqIO.read(fna_file, "fasta")

def extract_region(seq, position, window=20):
    start = max(0, position - window)  # Prevent negative indexing
    end = position + window
    return seq[start:end]

for gRNA_pos in gRNA_positions:
    extracted_seq = extract_region(record.seq, gRNA_pos, window=20)
    print(f"Sequence around gRNA at {gRNA_pos}: {extracted_seq}")

