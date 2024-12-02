include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'


workflow {
    vcf_ch = Channel.fromPath(params.vcf_file, checkIfExists: true )
    filter_regions(vcf_ch, file(params.bed_file)).view().set{ valid_vcfs } 
    filter_quality(valid_vcfs).view()
}