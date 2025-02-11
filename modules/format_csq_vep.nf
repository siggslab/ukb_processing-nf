process format_csq_vep {
    container "${params.containers.pandas}"
    //publishDir params.results_dir, mode: 'copy'
    
    input:
    path(input_veps)
    path python_script

    output:
    path("*.formatted_vep.vcf")

    script:
    """
    # Loop over each VCF in the batch
    for vcf in ${input_veps}
    do
        # Define the output file name based on the sample name
        id=\$(basename \${vcf} | cut -d. -f1)
        formatted_output_file="\${id}.formatted_vep.vcf"

        # Run the Python script to reformat the CSQ field
        python3 ${python_script} --input "\${vcf}" \\
                                --output "\${formatted_output_file}"
    done
    """
}