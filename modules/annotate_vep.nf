process annotate_vep {
    container "${params.containers.ensembl_vep}"

    input:
    path(input_vcfs)

    output:
    path("*.annotated_vep.vcf")
    
    script:
    """
    # Loop over each VCF in the batch
    for vcf in ${input_vcfs}
    do
        # Define the output file name based on the sample name
        id=\$(basename \${vcf} | cut -d. -f1)
        output_file="\${id}.annotated_vep.vcf"

        # Run VEP annotation
        vep \\
            --cache \\
            --offline \\
            --dir_cache "${params.vep_cache}" \\
            --assembly GRCh38 \\
            --input_file "\${vcf}" \\
            --output_file "\${output_file}" \\
            --vcf \\
            --force_overwrite \\
            --pick \\
            --everything 
    done
    """
}
