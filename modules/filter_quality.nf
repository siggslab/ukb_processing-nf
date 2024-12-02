process filter_quality {
    container "${params.containers.bcftools}"
    publishDir params.results_dir, mode: 'copy'
    input:
    path input_vcf    

    output:
    path "*.filtered.vcf.gz", optional: true

    script:
    
    """

    if [ -z "${input_vcf}" ]; then
        exit 0
    fi

    # Extract the sample name from the VCF file using bcftools
    sample_name=\$(bcftools query -l "${input_vcf}" | head -n 1)

    # Define the output file name based on the sample name
    output_file="\${sample_name}.filtered.vcf.gz"

    # Run bcftools filter with thresholds
    bcftools filter -i 'FORMAT/DP >=10' "${input_vcf}"\\
        -o "\${output_file}" -Oz 

    # Index the output VCF
    bcftools index "\${output_file}"
    """
}
