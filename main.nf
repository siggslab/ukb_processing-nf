include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'
include { annotate_vep } from './modules/annotate_vep'
include { annotate_annovar } from './modules/annotate_annovar'
include { format_vcf } from './modules/format_and_merge'
include { merge_vcf } from './modules/format_and_merge'
include { merge_phenotypes } from './modules/merge_phenotypes'
include { combine_csvs } from './modules/combine_csvs'
include { format_vep } from './modules/format_vep'
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
    //Reformat VEP Anno
    new_vep_ch = format_vep(vep_ch, file(params.transform_vep_anno_script))
    // Annotate Annovar
    annovar_ch = annotate_annovar(quality_filtered_ch)
    // Format the vcf for merge 
    (formatted_vep, formatted_annovar) = format_vcf(new_vep_ch, annovar_ch)
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

process make_results_dir {
    input:
    val results_dir

    script:
    """
    mkdir -p ${results_dir}
    echo ${results_dir}
    """
}

workflow {
    // Derive the batch name from params.vcf_file
    def batch_name = params.vcf_file
        .replaceFirst("^${baseDir}/", '') // Remove $baseDir prefix
        .replaceFirst("/\\*\\*.*", '')   // Remove trailing pattern
        .replaceAll("/", "_")            // Replace / with _
    // Create results directory if it doesn't exist
    make_results_dir(params.results_dir)
    // Find all VCF files in the batch_29 directory and its subdirectories
    vcf_ch = Channel.fromPath(params.vcf_file, checkIfExists: true)
    // Process all the VCFs individually in parallel
    csvs_ch = process_vcf(vcf_ch)
    // Combine the output CSV files
    collected_csvs = csvs_ch.collectFile(name: 'participant_csvs', tempDir: params.results_dir)
    combine_csvs(collected_csvs, batch_name).view()

}