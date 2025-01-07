process filter_regions {
    container "${params.containers.bcftools}"
    input:
    path vcf_file       
    path bed_file       

    output:
    path "filtered.vcf.gz", optional: true

    script:
    """
    output_file="filtered.vcf.gz"

    bcftools index "${vcf_file}"
    bcftools view -R "${bed_file}" -i 'FILTER="PASS"' "${vcf_file}" -o "output.vcf.gz" -Oz

    # Check if there are any variants (excluding the header)
    variant_count=\$(bcftools view -H "output.vcf.gz" | wc -l)

    # If no variants, do not create the output file
    if [ "\${variant_count}" -eq 0 ]; then
        exit 0
    else
        # Save the filtered VCF output if variants are present
        mv "output.vcf.gz" "\${output_file}"
        bcftools index "\${output_file}"  # Index the VCF if it contains variants
    fi
    """
}

