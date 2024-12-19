include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'
include { annotate_vep } from './modules/annotate_vep'
include { annotate_annovar } from './modules/annotate_annovar'
include { format_vcf } from './modules/format_and_merge'
include { merge_vcf } from './modules/format_and_merge'

workflow {
    // Load input VCF
    vcf_ch = Channel.fromPath(params.vcf_file, checkIfExists: true)
    // Filter VCF by regions defined in BED file
    filtered_regions_ch = filter_regions(vcf_ch, file(params.bed_file)).view()
    // Filter VCF by quality
    quality_filtered_ch = filter_quality(filtered_regions_ch).view()
    // Annotate VEP
    vep_ch = annotate_vep(quality_filtered_ch).view()
    // Annotate Annovar
    annovar_ch = annotate_annovar(quality_filtered_ch).view()
    // Format the vcf for merge 
    (formatted_vep, formatted_annovar) = format_vcf(vep_ch, annovar_ch)
    // Merge vep and Annovar annotations
    merge_vcf(formatted_vep, formatted_annovar)
}