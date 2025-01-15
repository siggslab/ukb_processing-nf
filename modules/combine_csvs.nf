process combine_csvs {
    publishDir params.results_dir, mode: 'copy'

    input:
    path participant_csvs  // List of CSV files

    output:
    path "*UKB_mosaic_variants.tsv"

    script:
    """
    # Combine all CSV files into one, keeping only the header from the first file
    # Extract the first file from the list
    first_file=\$(echo "${participant_csvs}" | cut -d' ' -f1)  # Get the first file in the list

    # Initialize the output file by writing the header from the first file
    head -n 1 "\$first_file" > UKB_mosaic_variants.tsv

    # Append the content of all files to the output, excluding headers from the second file onward
    for file in ${participant_csvs}; do
        tail -n +2 "\$file" >> UKB_mosaic_variants.tsv
    done
    """
}