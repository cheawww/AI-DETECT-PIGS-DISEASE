import matplotlib.pyplot as plt
import pandas as pd

# Example data for gRNA (could be modified as per your CSV data)
grna_data = pd.DataFrame({
    'gRNA Sequence': ['AGCT', 'GATC', 'TACG', 'CGTA'],
    'Start Position': [1, 50, 100, 150],
    'End Position': [30, 80, 130, 180],
    'Score': [80, 95, 75, 85]
})

# Create a bar plot of gRNA scores
plt.figure(figsize=(8, 6))
plt.bar(grna_data['gRNA Sequence'], grna_data['Score'], color='skyblue')
plt.xlabel('gRNA Sequence')
plt.ylabel('Score')
plt.title('gRNA Scores')
plt.tight_layout()

# Save the plot as an image
plt.savefig('grna_scores.png')
