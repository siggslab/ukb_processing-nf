process annotate_vep {
    container "${params.containers.ensembl_vep}"
    publishDir params.results_dir, mode: 'copy'

    input:
    path input_vcf    

    output:
    path "*.annotated.vcf", optional: true
    path "*.annotated.vcf_summary.html", optional: true

    script:
    """
    # Extract the sample name from the VCF file
    sample_name=\$(basename "${input_vcf}" | sed -E 's/\\.filtered\\.vcf\\.gz\$|\\.vcf\\.gz\$//')

    # Define the output file name based on the sample name
    output_file="\${sample_name}.annotated.vcf"

    # Run VEP annotation
    vep \\
        --cache \\
        --offline \\
        --dir_cache "${params.vep_cache}" \\
        --assembly GRCh38 \\
        --input_file "${input_vcf}" \\
        --output_file "\${output_file}" \\
        --vcf \\
        --force_overwrite \\
    """
}
