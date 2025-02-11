process merge_phenotypes {
    container "${params.containers.pandas}"
    //publishDir params.results_dir, mode: 'copy'
    input:
    path(variant_csvs)
    path csv_list_file
    path python_script

    output:
    path "output_*.tsv"

    script:
    """
    # Loop over each VCF in the batch
    for csv in ${variant_csvs}
    do
        # Define the output file name based on the sample name
        id=\$(basename \${csv} | cut -d. -f1)
        phenotypic_csvs=\$(cat ${csv_list_file} | tr -d '\\r' | tr '\\n' ' ')
        python3 ${python_script} --variant_csv \${csv} \\
                                --phenotype_csv_files \${phenotypic_csvs} \\
                                --output_file output_\${id}.tsv 
    done
    """
}