import argparse
import re

# List of genes to include
GENE_LIST = {
    "STING1", "LSM11", "RNU7-1", "CDC42", "STAT2", 
    "ATAD3A", "C2orf69", "RIPK1", "NCKAP1L", 
    "HCK1", "PSMB9", "IKBKG", "TBK1", "UBA1"
}

def fix_bed_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split line using any whitespace or inconsistent delimiters
            fields = re.split(r'\s+|\t+', line.strip())
            
            # Ensure there are at least 4 columns
            if len(fields) < 4:
                continue
            
            # Extract the gene name (the 4th column)
            gene = fields[3]
            
            # Skip lines that don't correspond to the genes in the list
            if gene not in GENE_LIST:
                continue
            
            # Add "chr" to the chromosome name if not already present
            chrom = fields[0]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            
            # Extract and format the first three columns
            chrom_start = fields[1]
            chrom_end = fields[2]
            fixed_line = f"{chrom}\t{chrom_start}\t{chrom_end}\t{gene}\n"
            
            outfile.write(fixed_line)

def main():
    parser = argparse.ArgumentParser(description="Fix a BED file for bcftools compatibility.")
    parser.add_argument("input_file", help="Path to the input BED file")
    parser.add_argument("output_file", help="Path to the output fixed BED file")
    args = parser.parse_args()
    
    fix_bed_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()