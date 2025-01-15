process merge_phenotypes {
    container "${params.containers.pandas}"
    publishDir params.results_dir, mode: 'copy'
    input:
    tuple val(id), path(variant_csv)
    path csv_list_file
    path python_script

    output:
    path "output_*.csv"

    script:
    """
    phenotypic_csvs=\$(cat ${csv_list_file} | tr -d '\\r' | tr '\\n' ' ')
    python3 ${python_script} --variant_csv ${variant_csv} \\
                            --phenotype_csv_files \${phenotypic_csvs} \\
                            --output_file output_${id}.csv \\
    """
}