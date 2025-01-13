process format_vep {
    container "${params.containers.pandas}"

    input:
    path input_vep
    path python_script

    output:
    path "*.formatted_vep.vcf"

    script:
    """
    # Extract the sample name from the VEP input
    sample=\$(basename "${input_vep}" | awk -F'.' '{print \$1}')
    formatted_output_file="\${sample}.formatted_vep.vcf"
    # Run the Python script to reformat the CSQ field
    python3 ${python_script} --input "${input_vep}" \\
                            --output "\${formatted_output_file}"
    """
}