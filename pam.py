from Bio import SeqIO
import pandas as pd

genome_fasta = "data/GCF_000003025.6_Sscrofa11.1_genomic.fna"
gff_file = "data/GCF_000003025.6_Sscrofa11.1_genomic.gff"
valid_PAMs = {"GGG", "GGA", "GGT", "GGC"}
gRNAs = [
    {"sequence": "TGAGTAGAGTCAGTCATGGA", "position": 35405522},
    {"sequence": "GAGTAGAGTCAGTCATGGAG", "position": 35405523},
    # Add more gRNAs here
]

def extract_pam_context(gRNA, genome_file, flank=20):
    for record in SeqIO.parse(genome_file, "fasta"):
        start = max(0, gRNA["position"] - flank)
        end = gRNA["position"] + len(gRNA["sequence"]) + flank
        region = str(record.seq[start:end].upper())
        pam = region[flank + len(gRNA["sequence"]):flank + len(gRNA["sequence"]) + 3]
        return region, pam
        
def load_gff(file):
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_data = pd.read_csv(file, sep="\t", comment="#", names=columns)
    return gff_data[gff_data["type"] == "gene"]

def check_overlap(gRNA, gff_data):
    overlaps = gff_data[(gff_data["start"] <= gRNA["position"]) & (gff_data["end"] >= gRNA["position"])]
    return overlaps
    
gff_data = load_gff(gff_file)
results = []
for gRNA in gRNAs:
    region, pam = extract_pam_context(gRNA, genome_fasta)
    is_valid_pam = pam in valid_PAMs
    overlaps = check_overlap(gRNA, gff_data)
    result = {
        "gRNA": gRNA["sequence"],
        "Position": gRNA["position"],
        "PAM": pam,
        "Is_Valid_PAM": is_valid_pam,
        "Context": region,
        "Overlapping_Genes": overlaps["attributes"].tolist() if not overlaps.empty else None,
    }
    results.append(result)
    
for result in results:
    print(f"gRNA: {result['gRNA']}, Position: {result['Position']}, PAM: {result['PAM']} ({'Valid' if result['Is_Valid_PAM'] else 'Invalid'})")
    print(f"Context: {result['Context']}")
    if result["Overlapping_Genes"]:
        print(f"Overlapping Genes: {', '.join(result['Overlapping_Genes'])}")
    else:
        print("No overlapping genes.")
    print("-" * 50)
    
output_file = "gRNA_validation_results.csv"
pd.DataFrame(results).to_csv(output_file, index=False)
print(f"Results saved to {output_file}")

