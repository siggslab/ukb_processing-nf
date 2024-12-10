
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
    path "*.annotated_merged.vcf.gz"

    script:
    """
    # Extract the sample name from the VEP input
    sample=\$(basename "${input_vep}" | awk -F'.' '{print \$1}')
    
    # Index the VCF files 
    bcftools index "${input_vep}"
    bcftools index "${input_annovar}"

    # Perform the merge using bcftools
    bcftools merge "${input_vep}" "${input_annovar}" -o "\${sample}.annotated_merged.vcf.gz" -Oz
    """
}