process format_vcf {
    container "${params.containers.bcftools}"
    input:
    path input_vep
    path input_annovar 

    output:
    path "*.vep_formatted.vcf.gz"
    path "*.annovar_formatted.vcf.gz"
    
    script:
    """
    # Extract sample name from input VEP file
    sample=\$(basename "${input_vep}" | awk -F'.' '{print \$1}')

    # Define the output file name based on the sample name
    output_vep="\${sample}.vep_formatted.vcf.gz"
    output_annovar="\${sample}.annovar_formatted.vcf.gz"

    # Process VEP file
    awk 'BEGIN {OFS="\\t"} 
    /^##INFO=<.*ID=/ { gsub(/-/, "_", \$0); gsub(/\\+\\+/, "PLUS", \$0); print; next } 
    /^#/ { print; next } 
    { gsub(/-/, "_", \$8); gsub(/\\+\\+/, "PLUS", \$8); print }' "${input_vep}" > vep_fixed.vcf

    awk -v sample="\$sample" 'BEGIN {OFS="\t"} {if (\$0 ~ /^#CHROM/) { \$10=sample"_vep"; print } 
    else {print}}' vep_fixed.vcf > vep_renamed.vcf

    # Process ANNOVAR file
    awk 'BEGIN {OFS="\\t"} 
    /^##INFO=<.*ID=/ { gsub(/-/, "_", \$0); gsub(/\\+\\+/, "PLUS", \$0); print; next } 
    /^#/ { print; next } 
    { gsub(/-/, "_", \$8); gsub(/\\+\\+/, "PLUS", \$8); print }' "${input_annovar}" > annovar_fixed.vcf

    awk -v sample="\$sample" 'BEGIN {OFS="\\t"} {if (\$0 ~ /^#CHROM/) { \$10=sample"_annovar"; print } 
    else {print}}' annovar_fixed.vcf > annovar_renamed.vcf

    # Compress files with bcftools
    bcftools view vep_renamed.vcf -o "\$output_vep" -Oz
    bcftools view annovar_renamed.vcf -o "\$output_annovar" -Oz
    """
}

process merge_vcf {
    container "${params.containers.bcftools}"
    publishDir params.results_dir, mode: 'copy'
    input:
    path input_vep
    path input_annovar 

    output:
    path "*_annotations.csv"
    script:
    """
    # Extract the sample name from the VEP input
    sample=\$(basename "${input_vep}" | awk -F'.' '{print \$1}')
    # Extract the Participant ID (second part of the sample name)
    participant_id=\$(echo "\$sample" | awk -F'_' '{print \$2}')

    # Index the VCF files 
    bcftools index "${input_vep}"
    bcftools index "${input_annovar}"

    # Perform the merge using bcftools
    bcftools merge "${input_vep}" "${input_annovar}" -o "\${sample}.annotated_merged.vcf.gz" -Oz

    # Query the merged VCF and add column headers
    output_csv="\${sample}_annotations.csv"

    # Write the column names and extract data from the VCF
    (
    # Print the header line with Participant ID as the first column
    echo -e "Participant ID\tCHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tFORMAT_INFO"

    # Extract the VCF data and process it
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO\t%FORMAT\n" "\${sample}.annotated_merged.vcf.gz" | \
    awk -v pid="\${participant_id}" '
        {
        # FORMAT column is the header (GT:AD:AF:DP:F1R2:F2R1:SB)
        format_header = \$(NF-2);
        format_info = \$NF;
        # Print the row with Participant ID and the split FORMAT and FORMAT_INFO columns
        print pid "\t" \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 "\t" format_header "\t" format_info;
        }
    '
    ) > "\${output_csv}"
    """
}





    
    