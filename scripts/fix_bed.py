import argparse
import re

# List of genes to include
GENE_LIST = {
    "STING1", "LSM11", "RNU7-1", "CDC42", "STAT2", 
    "ATAD3A", "C2orf69", "RIPK1", "NCKAP1L", 
    "HCK1", "PSMB9", "IKBKG", "TBK1", "UBA1", "ADA2",
    "TREX1", "SAMHD1", "IFIH1", "DNASE2", "LSM11", "ACP5",
    "POLA1", "USP18", "OAS1", "MEFV", "MVK", "NLRP3",
    "NLRP12", "NLRC4", "PLCG2", "NLRP1", "TNFRSF1A",
    "PSTPIP1", "NOD2", "ADAM17", "LPIN2", "IL1RN",
    "IL36RN","SLC29A3", "CARD14", "SH3BP2", "PSMB8",
    "PSMG2", "COPA", "OTULIN", "TNFAIP3", "AP1S3",
    "ALPI", "TRIM22", "HAVCR2", "SYK", "HCK"
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