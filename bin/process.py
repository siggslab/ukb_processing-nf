import pandas as pd
import re
import argparse
import os

def extract_and_clean_phenotypic_data(variant_csv, phenotype_csv_files, output_file):
    """
    Extracts phenotypic data for the Participant ID in variant CSV, cleans arrays and instances, 
    and combines all into one CSV.

    Args:
        variant_csv (str): Path to the variant CSV.
        phenotype_csv_files (list): List of "must-use" phenotypic CSV file paths.
        output_file (str): Output csv.
    """
    # Step 1: Load the variant CSV and get the Participant ID
    variant_df = pd.read_csv(variant_csv, sep='\t')
    participant_id = variant_df['Participant ID'].iloc[0]  # Assume one Participant ID

    # Initialize a list to store cleaned dataframes
    cleaned_dataframes = []

    # Step 2: Process each "must-use" phenotype CSV
    for phenotype_csv in phenotype_csv_files:
        df = pd.read_csv(phenotype_csv)
        if 'Participant ID' not in df.columns:
            print("HEWJADJAWDBJ")
            print(df.columns)
            raise KeyError("The column 'Participant ID' is missing in the variant CSV.")
        # Extract rows matching the Participant ID
        participant_data = df[df['Participant ID'] == participant_id]
        
        if participant_data.empty:
            print(f"No matching data found in {phenotype_csv} for Participant ID: {participant_id}")
            continue

        # Clean up instances and arrays
        cleaned_df = clean_instances_and_arrays(participant_data)
        cleaned_dataframes.append(cleaned_df)

    # Step 3: Merge all cleaned phenotypic dataframes on Participant ID
    if cleaned_dataframes:
        combined_phenotypic_df = cleaned_dataframes[0]
        for df in cleaned_dataframes[1:]:
            combined_phenotypic_df = combined_phenotypic_df.merge(df, on='Participant ID', how='outer')

    else:
        print("No phenotypic data found to combine.")
        # Create a blank phenotype DataFrame with columns from variant CSV and placeholders
        combined_phenotypic_df = pd.DataFrame(columns=cleaned_dataframes.columns)
        blank_row = {col: "" for col in combined_phenotypic_df.columns}
        blank_row["Participant ID"] = participant_id  # Ensure Participant ID is filled
        combined_phenotypic_df = combined_phenotypic_df.append(blank_row, ignore_index=True)
        
    # Step 4: Combine phenotypic data with the variant data
    final_combined_df = variant_df.merge(combined_phenotypic_df, on='Participant ID', how='left')
    # Save the final combined CSV
    final_combined_df.to_csv(output_file, index=False)
    print(f"Final combined CSV saved to {output_file}")
    
    if cleaned_dataframes is None or participant_id not in cleaned_dataframes:
            log_missing_phenotype(participant_id)

def log_missing_phenotype(participant_id, log_file_path="missing_phenotypes.txt"):
    """
    Logs participant IDs with missing phenotype data.

    Args:
        participant_id (str): The ID of the participant missing phenotype data.
        log_file_path (str): Path to the log file.
    """
    # Check if the log file exists
    file_exists = os.path.exists(log_file_path)

    # Open the file in append mode and write the log
    with open(log_file_path, "a") as log_file:
        # Write a header if the file is new
        if not file_exists:
            log_file.write("Missing Phenotype Log\n")
            log_file.write("Participant ID\n")
            log_file.write("====================\n")
        # Append the participant ID
        log_file.write(f"{participant_id}\n")

def clean_instances_and_arrays(df):
    # Step 1: Initialize final DataFrame and columns to maintain order
    final_columns = []
    new_df = pd.DataFrame()

    for col in df.columns:
        # Always keep 'Participant ID' as the first column
        if col == 'Participant ID':
            new_df['Participant ID'] = df['Participant ID'].fillna("")
            final_columns.append('Participant ID')

        # For instance columns: Keep only Instance 0
        elif 'Instance 0' in col and 'Array' not in col:
            clean_col_name = re.sub(r'\s*\|\s*Instance 0', '', col)  # Clean the column name
            if clean_col_name not in new_df.columns:
                new_df[clean_col_name] = df[col].fillna("")
                final_columns.append(clean_col_name)
        
        # For array columns: Concatenate arrays into a single column
        elif 'Array' in col:
            # If this column contains "Array", we treat it as an array column
            base_col_name = re.sub(r'\s*\|\s*(Instance \d+|\s*Array \d+)+', '', col)  # Base column name
            
            # If this is the first time encountering this base column name
            if base_col_name not in new_df.columns:
                # Find all relevant array columns for this base column (whether Instance 0 or not)
                array_columns = [c for c in df.columns if base_col_name in c]
                
                # Combine array columns into a single column
                combined_array = df[array_columns].apply(
                    lambda row: [val for val in row if pd.notna(val) and val != ""], axis=1
                )
                new_df[base_col_name] = combined_array.apply(
                    lambda x: x if len(x) > 0 else ""  # Replace empty lists with an empty string
                )
                final_columns.append(base_col_name)

        # Non-instance columns: Retain as-is and replace NaN with ""
        elif 'Instance' not in col and 'Array' not in col:
            new_df[col] = df[col].fillna("")
            final_columns.append(col)
    
    # Step 2: Write to a final CSV while preserving column order
    new_df = new_df[final_columns]  
    return new_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_csv", required=True)
    parser.add_argument("--phenotype_csv_files", nargs='+', required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()

    extract_and_clean_phenotypic_data(args.variant_csv, args.phenotype_csv_files, args.output_file)
    