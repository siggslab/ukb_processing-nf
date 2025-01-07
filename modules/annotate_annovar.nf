process annotate_annovar {
    container "${params.containers.annovar}"
    //publishDir params.results_dir, mode: 'copy'

    input:
    path input_vcf    

    output:
    path "*.annotated_annovar.vcf"

    script:
    """
    # Extract the sample name from the VCF file
    sample_name=\$(basename "${input_vcf}" | sed -E 's/\\.filtered\\.vcf\\.gz\$|\\.vcf\\.gz\$//')

    # Define the output file name based on the sample name
    output_file="\${sample_name}.annotated_annovar.vcf"

    # Run Annovar
    table_annovar.pl "${input_vcf}" \\
    "${params.annovar_db}" \\
    -buildver hg38 \\
    -out annotation \\
    -protocol dbnsfp42a \\
    -operation f \\
    -remove \\
    -vcfinput \\
    -nastring . \\

    mv annotation.hg38_multianno.vcf "\$output_file"
    """
}