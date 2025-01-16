include { filter_regions } from './modules/filter_regions'
include { filter_quality } from './modules/filter_quality'
include { annotate_vep } from './modules/annotate_vep'
include { annotate_annovar } from './modules/annotate_annovar'
include { format_annovar } from './modules/format_and_merge'
include { format_vep } from './modules/format_and_merge'
include { merge_vcf } from './modules/format_and_merge'
include { merge_phenotypes } from './modules/merge_phenotypes'
include { combine_csvs } from './modules/combine_csvs'
include { format_csq_vep } from './modules/format_csq_vep'
include { format_csvs } from './modules/format_csvs'

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
    // Reformat VEP Anno
    new_vep_ch = format_csq_vep(vep_ch, file(params.transform_vep_anno_script))
    // Annotate Annovar
    annovar_ch = annotate_annovar(quality_filtered_ch)
    // Format vep and annovar
    formatted_vep = format_vep(new_vep_ch)
    formatted_annovar = format_annovar(annovar_ch)
    // Merge the Annotations
    combined_ch = formatted_vep.join(formatted_annovar).view()
    merge_vcf_ch = merge_vcf(combined_ch)

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

    // Process all the VCFs individually in parallel
    tsvs_ch = process_vcf(vcf_ch)
    // Combine the output CSV files
    collected_tsvs = tsvs_ch.collect()
    combine_csvs(collected_tsvs).view()
}