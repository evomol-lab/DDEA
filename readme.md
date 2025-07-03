![](DDEA-small.png)
# DDEA - Diagonal Differential Expression Alley

DDEA is a simple Streamlit application to analyze Case x Control gene expression output data from GEO.
Allows users to upload a ```.tsv``` file, a custom gene list (optional),
and specify a comparison string, with controls in a sidebar.
DDEA then plots differentially expressed genes using Plotly and displays their p-values.

## How to Use the App

### App Sections and Functionality

The dashboard is divided into two main areas: a Sidebar for input and parameters, and a Main Content Area for data display and visualizations.

### Sidebar 

#### Upload Data & Set Parameters

This section allows you to upload your data files and configure the analysis parameters.

- Upload your DEG ```.tsv``` file:
    - Click the "Browse files" button to upload your gene expression data.
    - This file is expected to be in TSV (Tab-Separated Values) format from GEO2R.
        - It must contain a Symbol column for gene names, and columns for log2(fold change) and -log10(Pvalue) that include a specific comparison string in their headers (e.g., log2(fold change)(Case vs Control)).

#### Upload your custom gene list ```.txt``` file (optional):

You can optionally upload a plain text file (```.txt```) containing a list of gene symbols, with one gene symbol per line. If provided, the analysis will focus only on these genes. If no gene list is uploaded, the app will analyze the entire dataset and display the top 20 most differentially expressed genes. You can also specifiy the number of DEGs.

- Enter the comparison string (e.g., 'Case vs Control'):
    - This is a critical input. You must enter the exact string that appears within the parentheses in your log2(fold change) and -log10(Pvalue) column headers in your TSV file.
    - Example: If your column is log2(fold change)(Case vs Control), you should enter Case vs Control. The input is case-sensitive and space-sensitive.

### Filtering Thresholds:

**P-value threshold:** Set the maximum P-value for a gene to be considered differentially expressed. Genes with a P-value less than this threshold will be included. The default is 0.05.

**Minimum absolute Log2(FoldChange):** Set the minimum absolute Log2(FoldChange) value for a gene to be considered differentially expressed. Genes with an absolute Log2(FoldChange) greater than or equal to this value will be included. The default is 0.0, meaning no fold change filtering by default.

**Plot Customization:**

Maximum number of genes to plot (0 for all): Control how many genes are displayed in the bar plots.

- Enter 0 (default) to show all genes that meet the filtering criteria.
- Enter a positive integer to limit the plot to that number of top differentially expressed genes (based on absolute Log2(FoldChange)).

### Main Content Area: Data and Visualizations

This area displays the analysis results once the necessary files are uploaded and parameters are set.

1. Preview of your Raw Expression Data:

Shows the first few rows of your uploaded TSV file, allowing you to quickly inspect its format and content.

- Analysis Summary:
    - Provides key metrics about your data and the filtering process:
        - Genes from Custom List Found in Dataset / Total Unique Genes in Dataset: Indicates how many of your specified genes were found in the uploaded data, or the total number of unique genes if no custom list was provided.
        - Differentially Expressed Genes (Filtered): Shows the total count of genes that meet both the P-value and absolute Log2(FoldChange) thresholds.

2. Differentially Expressed Genes (Log2 Fold Change) for [Your Comparison]:

This is the main interactive bar plot generated using Plotly.

- X-axis: Gene Symbol.
- Y-axis: Log2(FoldChange).
- Title: Dynamically updates to include your specified comparison string and the applied filtering thresholds.
- Color: Bars are colored based on their Log2(FoldChange) value, typically showing a gradient from green (upregulated) to red (downregulated).
- Interactivity: Hover over bars to see exact values, zoom, and pan.
- Info Box: A helpful message clarifies that positive Log2(FoldChange) indicates upregulation and negative indicates downregulation relative to the control.

3. Details of Differentially Expressed Genes:

A table displaying the detailed information for the genes shown in the main plot.

- Includes: Gene Symbol, Log2(FoldChange), and -log10(Pvalue) (the original P-value format from your TSV).

4. Upregulated Genes (Log2 Fold Change) for [Your Comparison]:

A dedicated bar plot showing only the genes that are upregulated (Log2(FoldChange) > 0) and meet the filtering criteria. Sorted by Log2(FoldChange) in descending order.

5. Downregulated Genes (Log2 Fold Change) for [Your Comparison]:

A dedicated bar plot showing only the genes that are downregulated (Log2(FoldChange) < 0) and meet the filtering criteria. Sorted by absolute Log2(FoldChange) in descending order (meaning the most downregulated genes appear first).

6. Error Handling
The app includes specific error messages to guide you if:
- Required columns are not found in your TSV file (often due to an incorrect "comparison string").
- No genes from your custom list are found in the expression data.
- No differentially expressed genes are found after applying the filtering thresholds.
- An unexpected error occurs during data processing.

Remember to ensure your TSV file is correctly formatted and that the "comparison string" exactly matches your column headers for accurate analysis.