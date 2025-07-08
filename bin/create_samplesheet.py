import os
import csv
import argparse

def generate_vcf_csv(input_directory, output_csv):
    rows = []

    # Walk through the input directory and subdirectories
    for subdir, _, files in os.walk(input_directory):
        for file in files:
            if file.endswith("-filtered.vcf.gz"):
                sample_id = file.split('_')[0]
                vcf_filepath = os.path.abspath(os.path.join(subdir, file))
                rows.append([sample_id, vcf_filepath])

    # Determine how to split the file
    if len(rows) > 5000:
        part1_csv = output_csv.replace(".csv", "_part1.csv")
        part2_csv = output_csv.replace(".csv", "_part2.csv")

        with open(part1_csv, mode='w', newline='') as f1, open(part2_csv, mode='w', newline='') as f2:
            writer1 = csv.writer(f1)
            writer2 = csv.writer(f2)

            writer1.writerow(["sample_id", "sample_vcf"])
            writer2.writerow(["sample_id", "sample_vcf"])

            writer1.writerows(rows[:5000])
            writer2.writerows(rows[5000:])

        print(f"CSV files created: {part1_csv}, {part2_csv}")
    else:
        with open(output_csv, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["sample_id", "sample_vcf"])
            writer.writerows(rows)
        print(f"CSV file created at: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a CSV with sample IDs and VCF file paths.")
    parser.add_argument("input_directory", help="The directory to search for VCF files.")
    parser.add_argument("output_csv", help="The output CSV file path.")
    args = parser.parse_args()
    generate_vcf_csv(args.input_directory, args.output_csv)
