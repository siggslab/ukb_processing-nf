include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'
include { annotate_vep } from './modules/annotate_vep'
include { annotate_annovar } from './modules/annotate_annovar'
include { format_vcf } from './modules/format_and_merge'
include { merge_vcf } from './modules/format_and_merge'
include { merge_phenotypes } from './modules/merge_phenotypes'
include { combine_csvs } from './modules/combine_csvs'

workflow process_vcf {
    take:
    vcf_file
    main:
    // Filter VCF by regions defined in BED file
    filtered_regions_ch = filter_regions(vcf_file, file(params.bed_file))
    // Filter VCF by quality
    quality_filtered_ch = filter_quality(filtered_regions_ch)
    // Annotate VEP
    vep_ch = annotate_vep(quality_filtered_ch)
    // Annotate Annovar
    annovar_ch = annotate_annovar(quality_filtered_ch)
    // Format the vcf for merge 
    (formatted_vep, formatted_annovar) = format_vcf(vep_ch, annovar_ch)
    // Merge vep and Annovar annotations
    merge_vcf_ch = merge_vcf(formatted_vep, formatted_annovar)
    // Integrate Phenotype data
    phenotypic_csv_ch = merge_phenotypes(
        merge_vcf_ch,
        file(params.phenotype_csv),
        file(params.merge_phenotype_script)
    )
    emit:
    phenotypic_csv_ch
}

workflow {
    // Find all VCF files in the batch_29 directory and its subdirectories
    vcf_ch = Channel.fromPath(params.vcf_file, checkIfExists: true)
    // Process all the vcf individually in parallel
    csvs_ch = process_vcf(vcf_ch).view()
    // Combine the output csv 
    collected_csvs = csvs_ch.collectFile(name: 'participant_csvs', tempDir: params.results_dir)
    combine_csvs(collected_csvs)
}