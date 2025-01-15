import pandas as pd
import re
import argparse
import os


def extract_and_clean_phenotypic_data(variant_csv, phenotype_csv_files, output_file):
    """
    Extracts phenotypic data for the Participant ID in the variant CSV, cleans arrays and instances,
    and combines all into one CSV.

    Args:
        variant_csv (str): Path to the variant CSV.
        phenotype_csv_files (list): List of "must-use" phenotypic CSV file paths.
        output_file (str): Output CSV path.
    """
    # Step 1: Load the variant CSV and extract the Participant ID
    variant_df = pd.read_csv(variant_csv, sep='\t')
    participant_id = variant_df['Participant ID'].iloc[0]  # Assume one Participant ID

    # Initialize a dictionary to store combined data for each Participant ID
    combined_data = {}
    #missing_from_files = {}

    # Step 2: Process each "must-use" phenotype CSV
    for phenotype_csv in phenotype_csv_files:
        df = pd.read_csv(phenotype_csv)
        # Extract rows matching the Participant ID
        participant_data = df[df['Participant ID'] == participant_id]
        if participant_data.empty:
            # If no data is found, create an empty row with all columns from the CSV
            placeholder_row = {col: "" for col in phenotype_csv.columns} 
            placeholder_row["Participant ID"] = participant_id
            participant_data = pd.DataFrame([placeholder_row])
            '''
            # Log that the Participant ID is missing from this phenotype CSV
            if phenotype_csv not in missing_from_files:
                missing_from_files[phenotype_csv] = []
            missing_from_files[phenotype_csv].append(participant_id)
            '''
        # Clean up instances and arrays
        cleaned_df = clean_instances_and_arrays(participant_data)

        # Update combined data with cleaned DataFrame rows
        for _, row in cleaned_df.iterrows():
            pid = row['Participant ID']
            if pid not in combined_data:
                combined_data[pid] = row.to_dict()  # Initialize with row data
            else:
                # Update existing entry, ensuring no data is overwritten
                for col, val in row.items():
                    if col not in combined_data[pid] or pd.isna(combined_data[pid][col]) or combined_data[pid][col] == "":
                        combined_data[pid][col] = val

    # Step 3: Convert combined data into a DataFrame
    if combined_data:
        combined_phenotypic_df = pd.DataFrame.from_dict(combined_data, orient='index')
    else:
        print("No phenotypic data found to combine.")
        # Create a blank phenotype DataFrame with columns from variant CSV and placeholders
        combined_phenotypic_df = pd.DataFrame(columns=variant_df.columns)
        blank_row = {col: "" for col in combined_phenotypic_df.columns}
        blank_row["Participant ID"] = participant_id  # Ensure Participant ID is filled
        combined_phenotypic_df = combined_phenotypic_df.append(blank_row, ignore_index=True)

    # Merge phenotypic data with the variant data
    final_combined_df = variant_df.merge(combined_phenotypic_df, on='Participant ID', how='left')
    # Save the final combined CSV
    final_combined_df.to_csv(output_file,sep='\t', index=False)
    print(f"Final combined CSV saved to {output_file}")
    '''
    # Log any missing phenotype data
    if participant_id not in combined_data:
        log_missing_phenotype(participant_id,args.log_file_path)

    # Log missing Participant ID from specific phenotype CSVs
    for phenotype_csv, missing_ids in missing_from_files.items():
        for pid in missing_ids:
            log_missing_phenotype(pid,args.log_file_path,phenotype_csv=phenotype_csv)
    '''
def log_missing_phenotype(participant_id, log_file_path, phenotype_csv=None):
    """
    Logs participant IDs with missing phenotype data and the corresponding phenotype CSV file.

    Args:
        participant_id (str): The ID of the participant missing phenotype data.
        log_file_path (str): Path to the log file.
        phenotype_csv (str): Name of the phenotype CSV where the ID is missing (optional).
    """
    # Check if the log file exists
    file_exists = os.path.exists(log_file_path)

    # Open the file in append mode and write the log
    with open(log_file_path, "a") as log_file:
        # Write a header if the file is new
        if not file_exists:
            log_file.write("Missing Phenotype Log\n")
            log_file.write("Participant ID, Phenotype CSV\n")
            log_file.write("==============================\n")
        
        if phenotype_csv:
            # Log the missing participant ID and which phenotype CSV it was missing from
            log_file.write(f"{participant_id}, {phenotype_csv}\n")
        else:
            # Log the missing participant ID with no phenotype CSV (e.g., if no data was found for the participant)
            log_file.write(f"{participant_id}, N/A\n")

def clean_instances_and_arrays(df):
    """
    Cleans phenotypic data by processing instance and array columns.

    Args:
        df (pd.DataFrame): DataFrame containing participant phenotypic data.

    Returns:
        pd.DataFrame: Cleaned DataFrame.
    """
    # Initialize final DataFrame
    cleaned_df = pd.DataFrame()

    for col in df.columns:
        if col == 'Participant ID':
            # Always retain Participant ID as is
            cleaned_df['Participant ID'] = df['Participant ID'].fillna("")
        elif 'Instance 0' in col and 'Array' not in col:
            # Keep only Instance 0 and clean column name
            clean_col_name = re.sub(r'\s*\|\s*Instance 0', '', col)
            cleaned_df[clean_col_name] = df[col].fillna("")
        elif 'Array' in col:
            # Combine array columns into a single column
            base_col_name = re.sub(r'\s*\|\s*(Instance \d+|\s*Array \d+)+', '', col)
            if base_col_name not in cleaned_df.columns:
                array_cols = [c for c in df.columns if base_col_name in c]
                combined_array = df[array_cols].astype(str).apply(
                    lambda row: [x for x in row if x != 'nan' and x != ''], axis=1
                )
                cleaned_df[base_col_name] = combined_array.apply(lambda x: x if x else "")
        elif 'Instance' not in col and 'Array' not in col:
            # Retain non-instance, non-array columns as-is
            cleaned_df[col] = df[col].fillna("")

    return cleaned_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_csv", required=True)
    parser.add_argument("--phenotype_csv_files", nargs='+', required=True)
    parser.add_argument("--output_file", required=True)
    #parser.add_argument("--log_file_path", required=True)
    args = parser.parse_args()

    extract_and_clean_phenotypic_data(args.variant_csv, args.phenotype_csv_files, args.output_file)
