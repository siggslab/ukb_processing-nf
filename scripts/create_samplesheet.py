import os
import csv
import argparse

def generate_vcf_csv(input_directory, output_csv):
    # List to store the rows for the CSV file
    rows = []

    # Walk through the input directory and subdirectories
    for subdir, _, files in os.walk(input_directory):
        for file in files:
            # Only process .vcf.gz files
            if file.endswith("-filtered.vcf.gz"):
                # Extract the sample ID from the filename (before the first '_')
                sample_id = file.split('_')[0]
                
                # Get the absolute path of the VCF file
                vcf_filepath = os.path.abspath(os.path.join(subdir, file))
                
                # Append the sample ID and VCF path to the rows list
                rows.append([sample_id, vcf_filepath])

    # Write the rows to the CSV file
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["sample_id", "sample_vcf"])  # Write the header
        writer.writerows(rows)  # Write the data

    print(f"CSV file created at: {output_csv}")

# Argument parser for command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a CSV with sample IDs and VCF file paths.")
    parser.add_argument("input_directory", help="The directory to search for VCF files.")
    parser.add_argument("output_csv", help="The output CSV file path.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the function with the parsed arguments
    generate_vcf_csv(args.input_directory, args.output_csv)