import argparse
import re

def transform_csq(input_file, output_file):
    # Extract the CSQ header description
    csq_header_line = None
    with open(input_file, "r") as file:
        for line in file:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_header_line = line
                break

    if not csq_header_line:
        raise ValueError("No CSQ header found in the VCF file!")

    # Parse the CSQ field headers from the description
    csq_format = re.search(r"Format: (.+)>", csq_header_line).group(1)
    csq_headers = csq_format.split("|")

    # Process the VCF file
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            # Write header lines as is
            if line.startswith("#"):
                outfile.write(line)
                continue

            # Process data lines
            columns = line.strip().split("\t")
            info_field = columns[7]

            # Replace CSQ data with key-value pairs
            csq_data_match = re.search(r"CSQ=([^;]+)", info_field)
            if csq_data_match:
                csq_data = csq_data_match.group(1)
                csq_entries = csq_data.split(",")

                # Convert each CSQ entry to key-value format
                transformed_csq_entries = []
                for entry in csq_entries:
                    csq_values = entry.split("|")
                    transformed_entry = "|".join(
                        f"{key}={value}" for key, value in zip(csq_headers, csq_values)
                    )
                    transformed_csq_entries.append(transformed_entry)

                # Update INFO field
                transformed_csq = "CSQ=" + ",".join(transformed_csq_entries)
                info_field = re.sub(r"CSQ=[^;]+", transformed_csq, info_field)
                columns[7] = info_field

            # Write the updated line to the output
            outfile.write("\t".join(columns) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform VCF CSQ field into key-value format.")
    parser.add_argument("--input", required=True, help="Input annotated VEP VCF file.")
    parser.add_argument("--output", required=True, help="Output formatted VCF file.")
    args = parser.parse_args()

    transform_csq(args.input, args.output)