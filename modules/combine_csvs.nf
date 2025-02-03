process combine_csvs {
    publishDir params.results_dir, mode: 'copy'

    input:
    path tsv_files
    val batch_number

    output:
    path "${batch_number}_UKB_mosaic_variants.tsv"

    script:
    """
    # Extract the header from the first file
    header=\$(head -n 1 \$(ls ${tsv_files} | head -n 1))
    echo -e "\$header" > ${batch_number}_UKB_mosaic_variants.tsv

    # Append all TSV content, skipping headers except the first
    tail -n +2 -q ${tsv_files} >> ${batch_number}_UKB_mosaic_variants.tsv
    """
}