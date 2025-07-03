# DDEA - Diagonal Differential Expression Alley

import streamlit as st
import pandas as pd
import plotly.express as px # Import Plotly Express
import numpy as np
import os # Import os module to check for file existence

def run_gene_expression_analysis():
    """
    DDEA is a simple Streamlit application to analyze Case x Control gene expression output data from GEO.
    Allows users to upload a TSV file, a custom gene list (optional),
    and specify a comparison string, with controls in a sidebar.
    Plots differentially expressed genes using Plotly and displays their p-values.
    Includes a separate tab for documentation read from a readme.md file.
    """
    st.set_page_config(layout="wide") # Use wide layout for the dashboard

    st.title("Diagonal Differential Expression Alley")

    # Create tabs
    tab1, tab2 = st.tabs(["Analysis Dashboard", "Documentation"])

    with tab1:
        st.write("Upload your gene expression data and customize analysis parameters using the sidebar.")

        # --- Sidebar for File Uploads and User Inputs ---
        with st.sidebar:
            # Add the logo at the top of the sidebar
            st.image("DDEA-small.png", use_container_width=True) # Use use_column_width to make it responsive
            st.header("Upload Data & Set Parameters")

            use_example_files = st.checkbox("Use Example Files", value=False)

            if use_example_files:
                st.info("Using example files: 'Case vs Control.tsv' and 'GeneList.txt'")
                uploaded_tsv_file = "Examples/Case vs Control.tsv" # Path to example TSV
                uploaded_gene_list_file = "Examples/GeneList.txt" # Path to example gene list
                comparison_string = "Case vs Control" # Correct comparison for example TSV
                # Read example gene list content for processing
                with open(uploaded_gene_list_file, "rb") as f: # Open in binary mode for .decode()
                    example_gene_list_content = f.readlines()

            else:
                uploaded_tsv_file = st.file_uploader("Upload your DEG .tsv file", type=["tsv"])
                uploaded_gene_list_file = st.file_uploader("Upload your custom gene list .txt file (optional)", type=["txt"])

                comparison_string = st.text_input(
                    "Enter the comparison string (e.g., 'Case vs Control'). "
                    "**This must exactly match the string within the parentheses in your .tsv file's column headers** "
                    "for log2(fold change) and -log10(Pvalue). For example, if your column is 'log2(fold change)(Case vs Control)', you should enter 'Case vs Control'."
                )

            st.markdown("---") # Separator for clarity in sidebar

            st.subheader("Filtering Thresholds")
            p_value_threshold = st.number_input(
                "P-value threshold:",
                min_value=0.0,
                max_value=1.0,
                value=0.05,
                step=0.001,
                format="%.4f"
            )

            logFC_threshold = st.number_input(
                "Minimum absolute Log2(FoldChange):",
                min_value=0.0,
                value=0.0, # Default to 0, meaning no logFC filtering by default
                step=0.1,
                format="%.2f"
            )

            st.markdown("---") # Separator for clarity in sidebar

            st.subheader("Plot Customization")
            max_genes_to_plot = st.number_input(
                "Maximum number of genes to plot (0 for all):",
                min_value=0,
                value=0, # Default to 0, meaning show all genes
                step=1,
                help="Enter 0 to display all differentially expressed genes. Otherwise, specify a positive integer."
            )


        # --- Main Content Area for Data Preview, Plot, and Table ---
        # Check if files are available based on whether example files are used or uploaded
        if (use_example_files and os.path.exists(uploaded_tsv_file) and (os.path.exists(uploaded_gene_list_file) or True)) or \
           (not use_example_files and uploaded_tsv_file is not None and comparison_string):
            try:
                # Load the gene expression data
                # If using example files, load directly from path. Otherwise, use uploaded file object.
                if use_example_files:
                    df_expression = pd.read_csv(uploaded_tsv_file, sep='\t')
                else:
                    df_expression = pd.read_csv(uploaded_tsv_file, sep='\t')


                st.subheader("1. Preview of your Raw Expression Data:")
                st.dataframe(df_expression.head())

                # Define column names
                gene_col = 'Symbol'
                logFC_col = f'log2(fold change)({comparison_string})'
                neg_log10_p_value_col = f'-log10(Pvalue)({comparison_string})'

                # Validate if the expected columns exist in the DataFrame
                required_columns = [gene_col, logFC_col, neg_log10_p_value_col]
                missing_columns = [col for col in required_columns if col not in df_expression.columns]

                if missing_columns:
                    st.error(
                        f"Error: The following required columns were not found in your TSV file: "
                        f"{', '.join(missing_columns)}. "
                        f"Please check your 'Comparison String' input in the sidebar and ensure "
                        f"it exactly matches the column headers in your uploaded file."
                    )
                    st.write("Available columns in your file:", df_expression.columns.tolist())
                    return # Stop execution if columns are missing

                # Clean and standardize the gene symbols in the DataFrame
                df_expression[gene_col] = df_expression[gene_col].astype(str).str.strip().str.upper()

                # Initialize df_initial_filtered and count of genes from list
                df_initial_filtered = None
                genes_in_dataset_from_list = 0

                # Determine the DataFrame to filter based on gene list presence
                if use_example_files:
                    # Use content from example_gene_list_content
                    gene_list = [line.decode('utf-8').strip().upper() for line in example_gene_list_content]
                    df_initial_filtered = df_expression[df_expression[gene_col].isin(gene_list)].copy()
                    genes_in_dataset_from_list = df_initial_filtered[gene_col].nunique() # Count unique genes from list found in data
                elif uploaded_gene_list_file is not None:
                    gene_list = [uploaded_gene_list_file.read().decode('utf-8').strip().upper() for line in uploaded_gene_list_file.readlines()]
                    df_initial_filtered = df_expression[df_expression[gene_col].isin(gene_list)].copy()
                    genes_in_dataset_from_list = df_initial_filtered[gene_col].nunique() # Count unique genes from list found in data

                    if df_initial_filtered.empty:
                        st.warning("No genes from your custom list were found in the expression data after matching. Please double-check for typos, casing, or extra spaces in your gene list or the 'Symbol' column of your TSV file.")
                        return
                else:
                    df_initial_filtered = df_expression.copy()
                    genes_in_dataset_from_list = df_initial_filtered[gene_col].nunique() # Count all unique genes in the dataset
                    st.info("No custom gene list uploaded. Analyzing the entire dataset for differentially expressed genes.")


                # Convert -log10(Pvalue) to actual P-value
                df_initial_filtered['P_value_actual'] = 10**(-df_initial_filtered[neg_log10_p_value_col])

                # Filter for differentially expressed genes based on P-value and Log2(FoldChange) thresholds
                df_diff_expressed = df_initial_filtered[
                    (df_initial_filtered['P_value_actual'] < p_value_threshold) &
                    (df_initial_filtered[logFC_col].abs() >= logFC_threshold)
                ].copy()

                # --- Display Metrics ---
                st.subheader("Analysis Summary")
                col1, col2 = st.columns(2)
                with col1:
                    if uploaded_gene_list_file is not None or use_example_files: # Check if a list was used (either uploaded or example)
                        st.metric(label="Genes from Custom/Example List Found in Dataset", value=genes_in_dataset_from_list)
                    else:
                        st.metric(label="Total Unique Genes in Dataset", value=genes_in_dataset_from_list)
                with col2:
                    st.metric(label="Differentially Expressed Genes (Filtered)", value=len(df_diff_expressed))


                if df_diff_expressed.empty:
                    st.warning(
                        f"No differentially expressed genes found in your selection "
                        f"with a P-value less than {p_value_threshold} "
                        f"AND an absolute Log2(FoldChange) of at least {logFC_threshold}. "
                        f"Try adjusting the thresholds or checking if these genes are truly differentially expressed."
                    )
                else:
                    # Sort by absolute logFC for consistent ordering across all plots/tables
                    df_diff_expressed['logFC_abs'] = df_diff_expressed[logFC_col].abs()
                    df_diff_expressed = df_diff_expressed.sort_values(by='logFC_abs', ascending=False)

                    # Apply max_genes_to_plot if specified for the main plot and table
                    df_plot_main = df_diff_expressed.copy()
                    if max_genes_to_plot > 0:
                        df_plot_main = df_plot_main.head(max_genes_to_plot)

                    st.subheader(f"2. Differentially Expressed Genes (Log2 Fold Change) for {comparison_string}")
                    st.info("The Log2 Fold Change values represent the expression relative to the control. Positive values indicate upregulation in the treatment group, while negative values indicate downregulation.")

                    # Plot the main logFC values using Plotly Express
                    fig_main = px.bar(
                        df_plot_main,
                        x=gene_col,
                        y=logFC_col,
                        title=f'All Differentially Expressed Genes (P-value < {p_value_threshold}, |LogFC| ≥ {logFC_threshold})',
                        labels={gene_col: 'Gene Symbol', logFC_col: 'Log2(FoldChange)'},
                        color=logFC_col,
                        color_continuous_scale=px.colors.sequential.Viridis
                    )
                    fig_main.update_xaxes(tickangle=90)
                    fig_main.update_layout(height=600)
                    st.plotly_chart(fig_main, use_container_width=True)

                    st.subheader("3. Details of Differentially Expressed Genes:")
                    st.dataframe(df_plot_main[[gene_col, logFC_col, neg_log10_p_value_col]].rename(columns={
                        gene_col: 'Gene Symbol',
                        logFC_col: 'Log2(FoldChange)',
                        neg_log10_p_value_col: '-log10(Pvalue)'
                    }))

                    # --- Upregulated Genes Plot ---
                    df_upregulated = df_diff_expressed[df_diff_expressed[logFC_col] > 0].copy()
                    if max_genes_to_plot > 0:
                        df_upregulated = df_upregulated.head(max_genes_to_plot)

                    if not df_upregulated.empty:
                        st.subheader(f"4. Upregulated Genes (Log2 Fold Change) for {comparison_string}")
                        fig_up = px.bar(
                            df_upregulated,
                            x=gene_col,
                            y=logFC_col,
                            title=f'Upregulated Genes (P-value < {p_value_threshold}, LogFC ≥ {logFC_threshold})',
                            labels={gene_col: 'Gene Symbol', logFC_col: 'Log2(FoldChange)'},
                            color_discrete_sequence=px.colors.sequential.Greens_r # Green for upregulation
                        )
                        fig_up.update_xaxes(tickangle=90)
                        fig_up.update_layout(height=600)
                        st.plotly_chart(fig_up, use_container_width=True)
                    else:
                        st.info("No upregulated genes found based on the current filters.")


                    # --- Downregulated Genes Plot ---
                    df_downregulated = df_diff_expressed[df_diff_expressed[logFC_col] < 0].copy()
                    if max_genes_to_plot > 0:
                        df_downregulated = df_downregulated.head(max_genes_to_plot)

                    if not df_downregulated.empty:
                        st.subheader(f"5. Downregulated Genes (Log2 Fold Change) for {comparison_string}")
                        fig_down = px.bar(
                            df_downregulated,
                            x=gene_col,
                            y=logFC_col,
                            title=f'Downregulated Genes (P-value < {p_value_threshold}, LogFC ≤ -{logFC_threshold})',
                            labels={gene_col: 'Gene Symbol', logFC_col: 'Log2(FoldChange)'},
                            color_discrete_sequence=px.colors.sequential.Reds # Red for downregulation
                        )
                        fig_down.update_xaxes(tickangle=90)
                        fig_down.update_layout(height=600)
                        st.plotly_chart(fig_down, use_container_width=True)
                    else:
                        st.info("No downregulated genes found based on the current filters.")


            except Exception as e:
                st.error(f"An unexpected error occurred during data processing: {e}")
                st.write("Please ensure your TSV file is correctly formatted and matches the expected column names and data types. If using example files, ensure they are present in the same directory as the app.")

        else:
            st.info("Please upload your .tsv file and enter the comparison string in the sidebar to begin the analysis. You can optionally upload a custom gene list, or select 'Use Example Files'.")

    with tab2:
        st.header("App Documentation")
        # Check if readme.md exists and display its content
        readme_file_path = "readme.md"
        if os.path.exists(readme_file_path):
            with open(readme_file_path, "r", encoding="utf-8") as f:
                documentation_content = f.read()
            st.markdown(documentation_content)
        else:
            st.warning(f"Documentation file '{readme_file_path}' not found in the same directory as the app. Please ensure it's present.")
            st.markdown("Please create a `readme.md` file in the same directory as this Streamlit app and paste the documentation content into it.")


if __name__ == '__main__':
    run_gene_expression_analysis()
