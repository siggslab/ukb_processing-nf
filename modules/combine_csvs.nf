process combine_csvs {
    publishDir params.results_dir, mode: 'copy'

    input:
    path tsv_files

    output:
    path "UKB_mosaic_variants.tsv"

    script:
    """
    # Extract the header from the first file
    header=\$(head -n 1 \$(ls ${tsv_files} | head -n 1))
    echo -e "\$header" > UKB_mosaic_variants.tsv

    # Append all TSV content, skipping headers except the first
    tail -n +2 -q ${tsv_files} >> UKB_mosaic_variants.tsv
    """
}