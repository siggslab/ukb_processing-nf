import argparse
import re

def fix_bed_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split line using any whitespace or inconsistent delimiters
            fields = re.split(r'\s+|\t+', line.strip())
            
            # Ensure there are at least 3 columns
            if len(fields) < 3:
                continue
            
            # Add "chr" to the chromosome name if not already present
            chrom = fields[0]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            
            # Extract and format the first three columns
            chrom_start = fields[1]
            chrom_end = fields[2]
            fixed_line = f"{chrom}\t{chrom_start}\t{chrom_end}\n"
            
            outfile.write(fixed_line)

def main():
    parser = argparse.ArgumentParser(description="Fix a BED file for bcftools compatibility.")
    parser.add_argument("input_file", help="Path to the input BED file")
    parser.add_argument("output_file", help="Path to the output fixed BED file")
    args = parser.parse_args()
    
    fix_bed_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()