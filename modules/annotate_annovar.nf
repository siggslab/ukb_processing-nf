process annotate_annovar {
    container "${params.containers.annovar}"
    //publishDir params.intermediate_dir, mode: 'copy'

    input:
    tuple val(id), path(input_vcf)    

    output:
    tuple val(id), path("*.annotated_annovar.vcf")

    script:
    """
    # Define the output file name based on the sample name
    output_file="${id}.annotated_annovar.vcf"

    # Run Annovar
    table_annovar.pl "${input_vcf}" \\
    "${params.annovar_db}" \\
    -buildver hg38 \\
    -out ${id} \\
    -protocol dbnsfp42a \\
    -operation f \\
    -remove \\
    -vcfinput \\
    -nastring . \\

    mv ${id}.hg38_multianno.vcf "\$output_file"
    """
}