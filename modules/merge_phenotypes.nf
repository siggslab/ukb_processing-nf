process merge_phenotypes {
    container "${params.containers.pandas}"
    publishDir params.results_dir, mode: 'copy'
    input:
    path variant_csv
    path csv_list_file
    path python_script

    output:
    path "*.csv"

    script:
    """
    participant_id=\$(awk -F'\\t' 'NR==2 {print \$1}' ${variant_csv})
    phenotypic_csvs=\$(cat ${csv_list_file} | tr -d '\\r' | tr '\\n' ' ')
    python3 ${python_script} --variant_csv ${variant_csv} \\
                            --phenotype_csv_files \${phenotypic_csvs} \\
                            --output_file \${participant_id}.csv
    """
}