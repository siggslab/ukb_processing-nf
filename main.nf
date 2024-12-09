include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'
include { annotate_vep } from './modules/annotate_vep'
include { annotate_annovar } from './modules/annotate_annovar'

workflow {
    // Load input VCF
    vcf_ch = Channel.fromPath(params.vcf_file, checkIfExists: true)
    // Filter VCF by regions defined in BED file
    filtered_regions_ch = filter_regions(vcf_ch, file(params.bed_file)).view()
    // Filter VCF by quality
    quality_filtered_ch = filter_quality(filtered_regions_ch).view()
    // Annotate VEP
    annotate_vep(quality_filtered_ch)
    // Annotate Annovar
    annotate_annovar(quality_filtered_ch)
}