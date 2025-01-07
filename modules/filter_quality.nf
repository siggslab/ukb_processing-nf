process filter_quality {
    container "${params.containers.bcftools}"
    //publishDir params.results_dir, mode: 'copy'
    input:
    path input_vcf    

    output:
    path "*.filtered.vcf.gz", optional: true

    script:
    
    """
    # Extract the sample name from the VCF file using bcftools
    sample_name=\$(bcftools query -l "${input_vcf}" | head -n 1)

    # Define the output file name based on the sample name
    output_file="\${sample_name}.filtered.vcf.gz"

    # Run bcftools filter with thresholds
    bcftools filter \\
        -i 'FORMAT/DP >= 20 && FORMAT/AD[0:1] >= 3 && FORMAT/F1R2[*] >= 1 && FORMAT/F2R1[*] >= 1 && N_ALT=1' \\
        "${input_vcf}" \\
        -o "output.vcf.gz" -Oz

    # Check if there are any variants (excluding the header)
    variant_count=\$(bcftools view -H "output.vcf.gz" | wc -l)

    # If no variants, do not create the output file
    if [ "\${variant_count}" -eq 0 ]; then
        exit 0
    else
        # Save the filtered VCF output if variants are present
        mv "output.vcf.gz" "\${output_file}"
        bcftools index "\${output_file}"
    fi
    """
}
