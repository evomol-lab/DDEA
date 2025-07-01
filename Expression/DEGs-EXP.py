import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the gene expression datafile:
df_expression = pd.read_csv("GSE205748 - SL-PsA vs Healthy.tsv", sep='\t')

# Load the gene list
with open("GeneList.txt", "r") as f:
    gene_list = [line.strip() for line in f]

# Define the correct column names based on the previous inspection
gene_col = 'Symbol'
logFC_col = 'log2(fold change)(SkinLesion-PsA vs HealthySkin)'
neg_log10_p_value_col = '-log10(Pvalue)(SkinLesion-PsA vs HealthySkin)'

# Filter for genes in the gene list
df_filtered = df_expression[df_expression[gene_col].isin(gene_list)]

# Convert -log10(Pvalue) to actual P-value
df_filtered['P_value_actual'] = 10**(-df_filtered[neg_log10_p_value_col])

# Filter for differentially expressed genes (e.g., P-value < 0.05)
p_value_threshold = 0.05
df_diff_expressed = df_filtered[df_filtered['P_value_actual'] < p_value_threshold]

# Check if df_diff_expressed is empty
if df_diff_expressed.empty:
    print(f"\nNo differentially expressed genes found in the provided list with P-value < {p_value_threshold}.")
else:
    # Sort by absolute logFC for better visualization of most changed genes
    df_diff_expressed['logFC_abs'] = df_diff_expressed[logFC_col].abs()
    df_diff_expressed = df_diff_expressed.sort_values(by='logFC_abs', ascending=False)

    # Plot the logFC values
    plt.figure(figsize=(12, 8))
    plt.bar(df_diff_expressed[gene_col], df_diff_expressed[logFC_col], color='skyblue')
    plt.xlabel('Gene Name')
    plt.ylabel('Log2 Fold Change (SkinLesion-PsA vs HealthySkin)')
    plt.title(f'Differentially Expressed Genes (Log2 Fold Change) from Gene List (P-value < {p_value_threshold})')
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('GSE205748 - SL-PsA vs Healthy.png')

    # Describe the p-value for each gene
    print("\nDifferentially Expressed Genes and their p-values:")
    for index, row in df_diff_expressed.iterrows():
        print(f"Gene: {row[gene_col]}, Log2 Fold Change: {row[logFC_col]:.4f}, P-value: {row['P_value_actual']:.4e}")