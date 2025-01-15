process format_csq_vep {
    container "${params.containers.pandas}"
    //publishDir params.results_dir, mode: 'copy'
    
    input:
    tuple val(id), path(input_vep)
    path python_script

    output:
    tuple val(id), path("*.formatted_vep.vcf")

    script:
    """
    formatted_output_file="${id}.formatted_vep.vcf"
    # Run the Python script to reformat the CSQ field
    python3 ${python_script} --input "${input_vep}" \\
                            --output "\${formatted_output_file}"
    """
}