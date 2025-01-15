import pandas as pd
import re
import argparse


def extract_and_clean_phenotypic_data(variant_csv, phenotype_csv_files, output_file):
    """
    Extracts phenotypic data for the Participant ID in the variant CSV, cleans arrays and instances,
    and combines all into one tab-delimited CSV.

    Args:
        variant_csv (str): Path to the variant tab-delimited CSV.
        phenotype_csv_files (list): List of comma-delimited phenotypic CSV file paths.
        output_file (str): Output tab-delimited CSV path.
    """
    # Step 1: Load the variant CSV (tab-delimited) and extract the Participant ID
    variant_df = pd.read_csv(variant_csv, sep="\t")
    participant_id = variant_df['Participant ID'].iloc[0]  # Assume one Participant ID

    # Initialize a dictionary to store combined data for each Participant ID
    combined_data = {}
    participant_found = False  # Flag to track if participant ID was found

    # Step 2: Process each phenotype CSV
    for phenotype_csv in phenotype_csv_files:
        for chunk in pd.read_csv(phenotype_csv, sep=',', chunksize=10000):  # Explicitly comma-delimited
            # Filter for the relevant Participant ID
            participant_data = chunk[chunk['Participant ID'] == participant_id]
            if not participant_data.empty:
                participant_found = True
                # Clean up instances and arrays
                cleaned_df = clean_instances_and_arrays(participant_data)
                for _, row in cleaned_df.iterrows():
                    pid = row['Participant ID']
                    if pid not in combined_data:
                        combined_data[pid] = row.to_dict()
                    else:
                        for col, val in row.items():
                            if col not in combined_data[pid] or pd.isna(combined_data[pid][col]) or combined_data[pid][col] == "":
                                combined_data[pid][col] = val

    # If participant not found in any phenotype file, add placeholder row
    if not participant_found:
        placeholder_row = {col: "" for col in variant_df.columns}
        placeholder_row["Participant ID"] = participant_id
        combined_data[participant_id] = placeholder_row

    # Step 3: Convert combined data into a DataFrame
    combined_phenotypic_df = pd.DataFrame.from_dict(combined_data, orient='index')

    # Merge phenotypic data with the variant data
    final_combined_df = variant_df.merge(combined_phenotypic_df, on='Participant ID', how='left')

    # Save the final combined CSV as tab-delimited
    final_combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Final combined CSV saved to {output_file}")


def clean_instances_and_arrays(df):
    """
    Cleans phenotypic data by processing instance and array columns.

    Args:
        df (pd.DataFrame): DataFrame containing participant phenotypic data.

    Returns:
        pd.DataFrame: Cleaned DataFrame.
    """
    cleaned_df = pd.DataFrame()
    for col in df.columns:
        if col == 'Participant ID':
            cleaned_df['Participant ID'] = df['Participant ID'].fillna("")
        elif 'Instance 0' in col and 'Array' not in col:
            clean_col_name = re.sub(r'\s*\|\s*Instance 0', '', col)
            cleaned_df[clean_col_name] = df[col].fillna("")
        elif 'Array' in col:
            base_col_name = re.sub(r'\s*\|\s*(Instance \d+|\s*Array \d+)+', '', col)
            if base_col_name not in cleaned_df.columns:
                array_cols = [c for c in df.columns if base_col_name in c]
                combined_array = df[array_cols].astype(str).apply(
                    lambda row: [x for x in row if x != 'nan' and x != ''], axis=1
                )
                cleaned_df[base_col_name] = combined_array.apply(lambda x: x if x else "")
        elif 'Instance' not in col and 'Array' not in col:
            cleaned_df[col] = df[col].fillna("")

    return cleaned_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--variant_csv", required=True, help="Path to the tab-delimited variant CSV file.")
    parser.add_argument("--phenotype_csv_files", nargs='+', required=True, help="List of comma-delimited phenotype CSV files.")
    parser.add_argument("--output_file", required=True, help="Path to save the final tab-delimited output file.")
    args = parser.parse_args()

    extract_and_clean_phenotypic_data(args.variant_csv, args.phenotype_csv_files, args.output_file)