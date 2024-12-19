import pandas as pd
import re
import os

def extract_and_clean_phenotypic_data(variant_csv, out_csv_files, output_dir):
    """
    Extracts phenotypic data for the Participant ID in variant CSV, cleans arrays and instances, 
    and combines all into one CSV.

    Args:
        variant_csv (str): Path to the variant CSV.
        out_csv_files (list): List of "must-use" phenotypic CSV file paths.
        output_dir (str): Directory where intermediate and final files will be saved.
    """
    # Step 1: Load the variant CSV and get the Participant ID
    variant_df = pd.read_csv(variant_csv, sep='\t')
    participant_id = variant_df['Participant ID'].iloc[0]  # Assume one Participant ID

    # Initialize a list to store cleaned dataframes
    cleaned_dataframes = []

    # Step 2: Process each "must-use" CSV
    for out_csv in out_csv_files:
        print(f"Processing {out_csv}...")
        df = pd.read_csv(out_csv)

        # Extract rows matching the Participant ID
        participant_data = df[df['Participant ID'] == participant_id]
        
        if participant_data.empty:
            print(f"No matching data found in {out_csv} for Participant ID: {participant_id}")
            continue

        # Clean up instances and arrays
        cleaned_df = clean_instances_and_arrays(participant_data)

        # Save intermediate cleaned CSV
        cleaned_csv_path = os.path.join(output_dir, f"cleaned_{os.path.basename(out_csv)}")
        cleaned_df.to_csv(cleaned_csv_path, index=False)
        print(f"Saved cleaned data to {cleaned_csv_path}")

        cleaned_dataframes.append(cleaned_df)

    # Step 3: Merge all cleaned phenotypic dataframes on Participant ID
    if cleaned_dataframes:
        combined_phenotypic_df = cleaned_dataframes[0]
        for df in cleaned_dataframes[1:]:
            combined_phenotypic_df = combined_phenotypic_df.merge(df, on='Participant ID', how='outer')

        # Step 4: Combine phenotypic data with the variant data
        final_combined_df = variant_df.merge(combined_phenotypic_df, on='Participant ID', how='left')

        # Save the final combined CSV
        final_output_path = os.path.join(output_dir, "final_combined_data.csv")
        final_combined_df.to_csv(final_output_path, index=False)
        print(f"Final combined CSV saved to {final_output_path}")
    else:
        print("No phenotypic data found to combine.")

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


# Example usage
variant_csv = "UKB_4033866_301565080_annotations.csv"   # Replace with your variant CSV
out_csv_files = ["test5.csv", "test1.csv", "test2.csv", "test3.csv", "test4.csv"]  # Must-use CSV files
output_dir = "output_files"  # Directory to save cleaned and final files

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Run the process
extract_and_clean_phenotypic_data(variant_csv, out_csv_files, output_dir)
