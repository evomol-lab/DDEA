# DDEA - Diagonal Differential Expression Alley

import streamlit as st
import pandas as pd
import plotly.express as px # Import Plotly Express
import numpy as np
import os # Import os module to check for file existence
import re # Import regular expression module

def run_gene_expression_analysis():
    """
    DDEA is a Streamlit application to analyze gene expression data.
    Allows users to upload a TSV file, provide a custom gene list via text input (optional),
    or upload a gene list file (optional). The comparison string is automatically detected.
    Plots differentially expressed genes using Plotly and displays their p-values.
    Includes a separate tab for documentation read from a readme.md file.
    Also includes an option to use example files.
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
            st.image("DDEA-small.png", use_container_width=True)
            st.header("Upload Data & Set Parameters")

            use_example_files = st.checkbox("Use Example Files", value=False)

            # Initialize variables for file paths and content
            uploaded_tsv_file_obj = None
            uploaded_gene_list_file_obj = None
            custom_gene_list_input = ""
            current_comparison_string = ""
            example_gene_list_content = None

            if use_example_files:
                st.info("Using example files: 'Case vs Control.tsv' and 'GeneList.txt'")
                uploaded_tsv_file_path = "Examples/Case vs Control.tsv"
                uploaded_gene_list_file_path = "Examples/GeneList.txt"
                # For example files, we know the comparison string
                current_comparison_string = "Case vs Control"

                if not os.path.exists(uploaded_tsv_file_path):
                    st.error(f"Example TSV file '{uploaded_tsv_file_path}' not found. Please ensure it's in the same directory as the app.")
                    return
                if os.path.exists(uploaded_gene_list_file_path):
                    with open(uploaded_gene_list_file_path, "rb") as f:
                        example_gene_list_content = f.readlines()
                else:
                    st.warning(f"Example gene list file '{uploaded_gene_list_file_path}' not found. The analysis will proceed without a custom gene list.")

            else:
                uploaded_tsv_file_obj = st.file_uploader("Upload your DEG .tsv file", type=["tsv"])
                
                uploaded_gene_list_file_obj = st.file_uploader("Upload your custom gene list .txt file (optional)", type=["txt"])

                custom_gene_list_input = st.text_area(
                    "Or paste your custom gene list here (one gene symbol per line, optional)",
                    height=200,
                    help="If no gene list is provided (neither file nor pasted), the top 20 most differentially expressed genes from the dataset will be shown."
                )
                # Removed the manual comparison_string input here


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
                value=0.0,
                step=0.1,
                format="%.2f"
            )

            st.markdown("---") # Separator for clarity in sidebar

            st.subheader("Plot Customization")
            max_genes_to_plot = st.number_input(
                "Maximum number of genes to plot (0 for all):",
                min_value=0,
                value=0,
                step=1,
                help="Enter 0 to display all differentially expressed genes. Otherwise, specify a positive integer."
            )


        # --- Main Content Area for Data Processing and Display ---
        df_expression = None
        if use_example_files:
            if os.path.exists(uploaded_tsv_file_path):
                df_expression = pd.read_csv(uploaded_tsv_file_path, sep='\t')
        elif uploaded_tsv_file_obj is not None:
            df_expression = pd.read_csv(uploaded_tsv_file_obj, sep='\t')

        if df_expression is not None: # Check if TSV is loaded
            try:
                st.subheader("1. Preview of your Raw Expression Data:")
                st.dataframe(df_expression.head())

                # Automatically detect comparison string from column names
                gene_col = 'Symbol'
                logFC_col_pattern = r'log2\(fold change\)\((.*)\)'
                neg_log10_p_value_col_pattern = r'-log10\(Pvalue\)\((.*)\)'

                # Find the actual logFC column and extract the comparison string
                found_logFC_col = None
                found_neg_log10_p_value_col = None
                
                for col in df_expression.columns:
                    match_logFC = re.match(logFC_col_pattern, col)
                    if match_logFC:
                        found_logFC_col = col
                        # Extract the comparison string from the regex group
                        current_comparison_string = match_logFC.group(1)
                        
                    match_pvalue = re.match(neg_log10_p_value_col_pattern, col)
                    if match_pvalue:
                        found_neg_log10_p_value_col = col

                if not found_logFC_col or not found_neg_log10_p_value_col:
                    st.error(
                        "Error: Could not automatically detect the comparison string from your TSV file. "
                        "Please ensure your 'log2(fold change)' and '-log10(Pvalue)' columns follow the format "
                        "'log2(fold change)(Your Comparison)' and '-log10(Pvalue)(Your Comparison)'."
                    )
                    st.write("Available columns in your file:", df_expression.columns.tolist())
                    return

                # Use the detected column names
                logFC_col = found_logFC_col
                neg_log10_p_value_col = found_neg_log10_p_value_col

                st.info(f"Automatically detected comparison: **{current_comparison_string}**")


                # Clean and standardize the gene symbols in the DataFrame
                df_expression[gene_col] = df_expression[gene_col].astype(str).str.strip().str.upper()

                # Initialize df_initial_filtered and count of genes from list
                df_initial_filtered = None
                genes_in_dataset_from_list = 0
                using_custom_list = False
                gene_list = []

                # Determine the DataFrame to filter based on gene list presence
                if use_example_files and example_gene_list_content:
                    gene_list = [line.decode('utf-8').strip().upper() for line in example_gene_list_content]
                    using_custom_list = True
                elif uploaded_gene_list_file_obj is not None: # Prioritize uploaded file over text area
                    uploaded_gene_list_file_obj.seek(0)
                    gene_list = [line.decode('utf-8').strip().upper() for line in uploaded_gene_list_file_obj.readlines()]
                    using_custom_list = True
                elif custom_gene_list_input.strip(): # Check if text area has content
                    gene_list = [s.strip().upper() for s in custom_gene_list_input.split('\n') if s.strip()]
                    using_custom_list = True

                if using_custom_list:
                    df_initial_filtered = df_expression[df_expression[gene_col].isin(gene_list)].copy()
                    genes_in_dataset_from_list = df_initial_filtered[gene_col].nunique()
                    if df_initial_filtered.empty:
                        st.warning("No genes from your custom list were found in the expression data after matching. Please double-check for typos, casing, or extra spaces in your gene list or the 'Symbol' column of your TSV file.")
                        return
                else:
                    df_initial_filtered = df_expression.copy()
                    genes_in_dataset_from_list = df_initial_filtered[gene_col].nunique()
                    st.info("No custom gene list provided. Analyzing the entire dataset for differentially expressed genes.")


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
                    if using_custom_list:
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

                    # If no custom list was provided AND max_genes_to_plot is 0, take top 20 DEGs
                    # Otherwise, apply max_genes_to_plot if specified for the main plot and table
                    df_plot_main = df_diff_expressed.copy()
                    if not using_custom_list and max_genes_to_plot == 0:
                        df_plot_main = df_plot_main.head(20)
                    elif max_genes_to_plot > 0:
                        df_plot_main = df_plot_main.head(max_genes_to_plot)


                    st.subheader(f"2. Differentially Expressed Genes (Log2 Fold Change) for {current_comparison_string}")
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
                        st.subheader(f"4. Upregulated Genes (Log2 Fold Change) for {current_comparison_string}")
                        fig_up = px.bar(
                            df_upregulated,
                            x=gene_col,
                            y=logFC_col,
                            title=f'Upregulated Genes (P-value < {p_value_threshold}, LogFC ≥ {logFC_threshold})',
                            labels={gene_col: 'Gene Symbol', logFC_col: 'Log2(FoldChange)'},
                            color_discrete_sequence=px.colors.sequential.Blues_r # Green for upregulation
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
                        st.subheader(f"5. Downregulated Genes (Log2 Fold Change) for {current_comparison_string}")
                        fig_down = px.bar(
                            df_downregulated,
                            x=gene_col,
                            y=logFC_col,
                            title=f'Downregulated Genes (P-value < {p_value_threshold}, LogFC ≤ -{logFC_threshold})',
                            labels={gene_col: 'Gene Symbol', logFC_col: 'Log2(FoldChange)'},
                            color_discrete_sequence=px.colors.sequential.Reds_r # Red for downregulation
                        )
                        fig_down.update_xaxes(tickangle=90)
                        fig_down.update_layout(height=600)
                        st.plotly_chart(fig_down, use_container_width=True)
                    else:
                        st.info("No downregulated genes found based on the current filters.")


            except Exception as e:
                st.error(f"An unexpected error occurred during data processing: {e}")
                st.write(f"Please check your input files and parameters. Specific error: {e}")

        else:
            st.info("Please upload your .tsv file to begin the analysis. The comparison string will be automatically detected. You can optionally paste a custom gene list, upload a .txt gene list file, or select 'Use Example Files'.")

    with tab2:
        st.header("App Documentation")
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
