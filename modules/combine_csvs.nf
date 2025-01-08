process combine_csvs {
    publishDir params.results_dir, mode: 'copy'

    input:
    file participant_csvs  // List of CSV files

    output:
    path "UKB_mosaic_variants.csv"

    script:
    """
    # Use awk to combine all CSV files into one, keeping the header from the first file only
    awk 'NR==1{print; next} FNR>1' ${participant_csvs} > UKB_mosaic_variants.csv
    """
}