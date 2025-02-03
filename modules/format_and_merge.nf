process format_annovar {
    container "${params.containers.bcftools}"
    //publishDir params.results_dir, mode: 'copy'
    
    input:
    tuple val(id), path(input_vcf)

    output:
    tuple val(id), path("*_annovar_formatted.vcf.gz")
    
    script:
    """
    # Define the output file name based on the sample name
    output_file="${id}_annovar_formatted.vcf.gz"

    awk 'BEGIN {OFS="\\t"} 
    /^##INFO=<.*ID=/ { gsub(/-/, "_", \$0); gsub(/\\+\\+/, "PLUS", \$0); print; next } 
    /^#/ { print; next } 
    { gsub(/-/, "_", \$8); gsub(/\\+\\+/, "PLUS", \$8); print }' "${input_vcf}" > ${id}_annovar.fixed.vcf

    awk -v sample="${id}" 'BEGIN {OFS="\t"} {if (\$0 ~ /^#CHROM/) { \$10=sample"_annovar"; print } 
    else {print}}' ${id}_annovar.fixed.vcf > ${id}_annovar.renamed.vcf

    # Compress files with bcftools
    bcftools view ${id}_annovar.renamed.vcf -o "\$output_file" -Oz

    """
}

process format_vep {
    container "${params.containers.bcftools}"
    //publishDir params.results_dir, mode: 'copy'
    
    input:
    tuple val(id), path(input_vcf)

    output:
    tuple val(id), path("*_vep_formatted.vcf.gz")
    
    script:
    """
    output_file="${id}_vep_formatted.vcf.gz"

    awk 'BEGIN {OFS="\\t"} 
    /^##INFO=<.*ID=/ { gsub(/-/, "_", \$0); gsub(/\\+\\+/, "PLUS", \$0); print; next } 
    /^#/ { print; next } 
    { gsub(/-/, "_", \$8); gsub(/\\+\\+/, "PLUS", \$8); print }' "${input_vcf}" > ${id}_vep.fixed.vcf

    awk -v sample="${id}" 'BEGIN {OFS="\t"} {if (\$0 ~ /^#CHROM/) { \$10=sample"_vep"; print } 
    else {print}}' ${id}_vep.fixed.vcf > ${id}_vep.renamed.vcf

    # Compress files with bcftools
    bcftools view ${id}_vep.renamed.vcf -o "\$output_file" -Oz

    """
}

process vcf_to_csv {
    container "${params.containers.bcftools}"
    //publishDir params.results_dir, mode: 'copy'
    input:
    //tuple val(id), path(input_vep), path(input_annovar)
    tuple val(id), path(input_vep)

    output:
    tuple val(id), path("*_annotations.csv")
    script:
    """
    # Index the VCF files 
    bcftools index "${input_vep}"

    # Query the merged VCF and add column headers
    output_csv="${id}_annotations.csv"

    # Write the column names and extract data from the VCF
    (
    # Print the header line with Participant ID as the first column
    echo -e "Participant ID\\tCHROM\\tPOS\\tREF\\tALT\\tINFO\\tFORMAT\\tFORMAT_INFO"

    # Extract the VCF data and process it
    bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO\\t%FORMAT\n" "${id}_vep_formatted.vcf.gz" | \
    awk -v pid="${id}" '
        {
        # FORMAT column is the header (GT:AD:AF:DP:F1R2:F2R1:SB)
        format_header = \$(NF-1);
        format_info = \$NF;
        # Print the row with Participant ID and the split FORMAT and FORMAT_INFO columns
        print pid "\\t" \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" format_header "\\t" format_info;
        }
    '
    ) > "\${output_csv}"
    """
}





    
    