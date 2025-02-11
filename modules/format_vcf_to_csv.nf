process format_vep_and_vcf_to_csv {
    container "${params.containers.bcftools}"

    input:
    path(input_vcfs)

    output:
    path("*.annotations.csv")

    script:
    """
    # Loop over each VCF in the batch
    for vcf in ${input_vcfs}
    do
        # Define the output file name based on the sample name
        id=\$(basename \${vcf} | cut -d. -f1)

        # Step 1: Format the VEP-annotated VCF
        output_vcf="\${id}_vep_formatted.vcf.gz"

        awk 'BEGIN {OFS="\\t"} 
        /^##INFO=<.*ID=/ { gsub(/-/, "_", \$0); gsub(/\\+\\+/, "PLUS", \$0); print; next } 
        /^#/ { print; next } 
        { gsub(/-/, "_", \$8); gsub(/\\+\\+/, "PLUS", \$8); print }' "\${vcf}" > \${id}_vep.fixed.vcf

        awk -v sample="\${id}" 'BEGIN {OFS="\t"} {if (\$0 ~ /^#CHROM/) { \$10=sample"_vep"; print } 
        else {print}}' \${id}_vep.fixed.vcf > \${id}_vep.renamed.vcf

        # Compress the formatted VCF
        bcftools view \${id}_vep.renamed.vcf -o "\$output_vcf" -Oz
        bcftools index -t "\$output_vcf"

        # Step 2: Convert the formatted VCF to CSV
        output_csv="\${id}.annotations.csv"

        # Write the column names and extract data from the VCF
        (
        echo -e "Participant ID\\tCHROM\\tPOS\\tREF\\tALT\\tINFO\\tFORMAT\\tFORMAT_INFO"

        bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO\\t%FORMAT\n" "\$output_vcf" | \
        awk -v pid="\${id}" '
            {
            format_header = \$(NF-1);
            format_info = \$NF;
            print pid "\\t" \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" format_header "\\t" format_info;
            }
        '
        ) > "\$output_csv"
    done
    """
}