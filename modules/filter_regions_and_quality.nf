process filter_regions_and_quality {
    container "${params.containers.bcftools}"
    
    input:
    //[[1,2,3],[1.vcf, 2.vcf, 3.vcf]]
    tuple val(ids), path(vcf_files)
    path bed_file

    output:
    path("*.filtered.vcf.gz"), optional: true

    script:
    if (ids instanceof List) {
        ids_str = ids.join(" ")
    } else {
        ids_str = ids
    }

    """
    IDS=(\$(echo $ids_str))
    VCF_FILES=(\$(echo $vcf_files))
    
    # Loop over each VCF in the batch
    for i in "\${!VCF_FILES[@]}"
    do
        id=\${IDS[\$i]}
        vcf_file=\${VCF_FILES[\$i]}

        # Define intermediate and final output file names
        regions_filtered_vcf="\${id}.filtered_regions.vcf.gz"
        final_filtered_vcf="\${id}.filtered.vcf.gz"

        # Step 1: Filter by regions
        bcftools index -t "\${vcf_file}"
        bcftools view -R "${bed_file}" -i 'FILTER="PASS"' "\${vcf_file}" -o "\${id}.output_regions.vcf.gz" -Oz

        # Check if any variants passed the region filter
        variant_count=\$(bcftools view -H "\${id}.output_regions.vcf.gz" | wc -l)

        if [ "\${variant_count}" -eq 0 ]; then
            echo "No variants remaining after region filtering. Exiting."
            exit 0
        else
            mv "\${id}.output_regions.vcf.gz" "\${regions_filtered_vcf}"
            bcftools index -t "\${regions_filtered_vcf}"
        fi

        # Step 2: Filter by quality
        bcftools filter \\
            -i 'FORMAT/DP >= 20 && FORMAT/AD[0:1] >= 3 && FORMAT/F1R2[*] >= 1 && FORMAT/F2R1[*] >= 1 & N_ALT=1 && FORMAT/DP != "." && FORMAT/AD !="." && FORMAT/F1R2 !="." && FORMAT/F2R1 !="."' \\
            "\${regions_filtered_vcf}" \\
            -o "\${id}.output_quality.vcf.gz" -Oz

        # Check if any variants passed the quality filter
        variant_count=\$(bcftools view -H "\${id}.output_quality.vcf.gz" | wc -l)

        if [ "\${variant_count}" -eq 0 ]; then
            echo "No variants remaining after quality filtering. Exiting."
            exit 0
        else
            mv "\${id}.output_quality.vcf.gz" "\${final_filtered_vcf}"
            bcftools index -t "\${final_filtered_vcf}"
        fi
    done
    """
}