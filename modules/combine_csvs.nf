process combine_csvs {
    publishDir params.results_dir, mode: 'copy'

    input:
    file participant_csvs  // List of CSV files
    val batch_name
    output:
    path "*_UKB_mosaic_variants.csv"

    script:
    """
    # Combine all CSV files into one, keeping only the header from the first file
    awk 'FNR==1 && NR!=1 {next} {print}' ${participant_csvs} > ${batch_name}_UKB_mosaic_variants.csv
    """
}