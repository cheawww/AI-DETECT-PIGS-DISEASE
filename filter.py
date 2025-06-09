import matplotlib.pyplot as plt

# Sample Data
positions = filtered_df['Position']
gRNAs = filtered_df['gRNA']

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(positions, range(len(positions)), color="blue", label="gRNA Positions")
plt.yticks(range(len(gRNAs)), gRNAs, fontsize=8)
plt.xlabel("Genomic Position")
plt.ylabel("gRNA")
plt.title("gRNA Positions and Overlapping Genes")
plt.legend()
plt.grid()
plt.tight_layout()

# Save Plot
plt.savefig("gRNA_positions_plot.png")
plt.show()
