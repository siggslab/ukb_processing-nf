process filter_regions {
    container "${params.containers.bcftools}"
    input:
    tuple val(id), path(vcf_file)      
    path bed_file       

    output:
    tuple val(id), path("*.filtered_regions.vcf.gz"), optional: true

    script:
    """

    # Define the output file name based on the sample name
    output_file="${id}.filtered_regions.vcf.gz"

    bcftools index -t "${vcf_file}"
    bcftools view -R "${bed_file}" -i 'FILTER="PASS"' "${vcf_file}" -o "${id}.output.vcf.gz" -Oz

    # Check if there are any variants (excluding the header)
    variant_count=\$(bcftools view -H "${id}.output.vcf.gz" | wc -l)

    # If no variants, do not create the output file
    if [ "\${variant_count}" -eq 0 ]; then
        exit 0
    else
        # Save the filtered VCF output if variants are present
        mv "${id}.output.vcf.gz" "\${output_file}"
        bcftools index -t "\${output_file}"  # Index the VCF if it contains variants
    fi
    """
}

