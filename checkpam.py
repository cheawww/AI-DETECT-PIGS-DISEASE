from Bio import SeqIO

# Input files
fna_file = "GCF_000003025.6_Sscrofa11.1_genomic.fna"
gRNAs = [
    {"sequence": "TTGCTTTGGTCTTTAAACTC", "position": 123991, "chromosome": "NC_010450.4"},
    {"sequence": "TCTTTAAACTCAGGAGCCAC", "position": 124000, "chromosome": "NC_010450.4"},
    {"sequence": "AACATCCTGCCCAGAACAGA", "position": 124029, "chromosome": "NC_010450.4"},
]
sequences = {record.id: record.seq for record in SeqIO.parse(fna_file, "fasta")}
for gRNA in gRNAs:
    chrom = gRNA["chromosome"]
    position = gRNA["position"]
    sequence = gRNA["sequence"]
    if chrom not in sequences:
        print(f"Chromosome {chrom} not found in the genome.")
        continue
    start = max(0, position - 10)
    end = position + 23
    context = sequences[chrom][start:end]
    pam_start = 20  # PAM is immediately downstream of the gRNA (20 bases)
    pam = context[pam_start:pam_start + 3]  # Extract 3 bases for PAM
    pam_valid = pam[0] in "AGCT" and pam[1:] == "GG"
    
    print(f"gRNA: {sequence}, Position: {position}")
    print(f"Context: {context}")
    print(f"PAM: {pam} {'(Valid)' if pam_valid else '(Invalid)'}")
    print("-" * 50)
