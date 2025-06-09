from pyfaidx import Fasta

fasta = Fasta("GCF_000003025.6_Sscrofa11.1_genomic.fna")
print(list(fasta.keys()))  # This prints all sequence names in the FASTA file


chrom = "chr1"
if chrom not in fasta.keys():
    print(f"Error: {chrom} not found in FASTA file.")
else:
    region_seq = fasta[chrom][100000:200000]
