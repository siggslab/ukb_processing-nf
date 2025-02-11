include { annotate_vep } from './modules/annotate_vep'
include { merge_phenotypes } from './modules/merge_phenotypes'
include { combine_csvs } from './modules/combine_csvs'
include { format_csq_vep } from './modules/format_csq_vep'
include { filter_regions_and_quality } from './modules/filter_regions_and_quality'
include { format_vep_and_vcf_to_csv } from './modules/format_vcf_to_csv'
workflow process_vcf {
    take:
    vcf_files
    main:
    // Filter VCF by regions defined in BED file and by quality
    quality_filtered_ch = filter_regions_and_quality(vcf_files, file(params.bed_file))
    // Annotate VEP
    vep_ch = annotate_vep(quality_filtered_ch)
    // Reformat VEP Annotations and turn it into a csv
    new_vep_ch = format_csq_vep(vep_ch, file(params.transform_vep_anno_script))
    merge_vcf_ch = format_vep_and_vcf_to_csv(new_vep_ch)
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
    // Create results directory if it doesn't exist
    make_results_dir(params.results_dir)
    csv = Channel.fromPath(params.vcf_csv, checkIfExists: true)
        .splitCsv(header: true)
    vcf_ch = csv.map { row ->
        [ row.sample_id, file(row.sample_vcf, checkIfExists: true) ]
    }
    batched_vcf_ch = vcf_ch.collate(9).map{ it -> it.transpose() }
    // Process all the VCFs in batches
    tsvs_ch = process_vcf(batched_vcf_ch)
    // Combine the output CSV files
    collected_tsvs = tsvs_ch.flatten().collect()
    // Extract batch number from the CSV file name (e.g., '10_samplesheet.csv' -> '10')
    file_name = params.vcf_csv.startsWith("/") ? params.vcf_csv.tokenize('/').last() : params.vcf_csv
    batch_number = file_name.split('_')[0]  // This assumes the batch number is the first part of the file name
    combine_csvs(collected_tsvs, batch_number)
}