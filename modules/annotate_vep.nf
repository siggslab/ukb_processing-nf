process annotate_vep {
    container "${params.containers.ensembl_vep}"
    //publishDir params.intermediate_dir, mode: 'copy'

    input:
    tuple val(id), path(input_vcf)

    output:
    tuple val(id), path("*.annotated_vep.vcf")
    
    script:
    """
    # Define the output file name based on the sample name
    output_file="${id}.annotated_vep.vcf"

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
        --pick \\
        --everything \\
    """
}
