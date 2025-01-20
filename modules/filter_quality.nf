process filter_quality {
    container "${params.containers.bcftools}"
    //publishDir params.results_dir, mode: 'copy'
    input:
    tuple val(id), path(input_vcf)

    output:
    tuple val(id), path("*.filtered.vcf.gz"), optional: true

    script:
    
    """
    # Define the output file name based on the sample name
    output_file="${id}.filtered.vcf.gz"

    # Run bcftools filter with thresholds
    bcftools filter \\
        -i 'FORMAT/DP >= 20 && FORMAT/AD[0:1] >= 3 && FORMAT/F1R2[*] >= 1 && FORMAT/F2R1[*] >= 1 & N_ALT=1 && FORMAT/DP != "." && FORMAT/AD !="." && FORMAT/F1R2 !="." && FORMAT/F2R1 !="."' \\
        "${input_vcf}" \\
        -o "${id}.output_quality.vcf.gz" -Oz

    # Check if there are any variants (excluding the header)
    variant_count=\$(bcftools view -H "${id}.output_quality.vcf.gz" | wc -l)

    # If no variants, do not create the output file
    if [ "\${variant_count}" -eq 0 ]; then
        exit 0
    else
        # Save the filtered VCF output if variants are present
        mv "${id}.output_quality.vcf.gz" "\${output_file}"
        bcftools index -t "\${output_file}"
    fi
    """
}
