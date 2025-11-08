from Bio import SeqIO
import re
fna_file = "data/GCF_000003025.6_Sscrofa11.1_genomic.fna"
gff_file = "data/GCF_000003025.6_Sscrofa11.1_genomic.gff"

genome = {record.id: record.seq for record in SeqIO.parse(fna_file, "fasta")}

def find_gRNAs(sequence, pam="NGG", gRNA_length=20):
    gRNAs = []
    pam_regex = pam.replace("N", ".")
    for match in re.finditer(pam_regex, str(sequence)):
        pam_start = match.start()
        if pam_start >= gRNA_length:
            gRNA = sequence[pam_start - gRNA_length:pam_start]
            gRNAs.append((str(gRNA), pam_start - gRNA_length))
    return gRNAs

with open(gff_file) as gff:
    for line in gff:
        if not line.startswith("#"):
            parts = line.strip().split("\t")
            if parts[2] == "gene":
                seq_id = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                gene_id = re.search(r"ID=([^;]+)", attributes).group(1)

                if seq_id in genome:
                    region_seq = genome[seq_id][start - 1:end]
                    if strand == "-":
                        region_seq = region_seq.reverse_complement()

                    gRNAs = find_gRNAs(region_seq)
                    print(f"Gene: {gene_id}, gRNAs found: {len(gRNAs)}")
                    for gRNA, pos in gRNAs:
                        print(f"gRNA: {gRNA}, Position: {pos}")

def gc_content(sequence):
    return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
