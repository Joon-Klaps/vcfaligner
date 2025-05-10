/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { softwareVersionsToYAML   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MAFFT_ALIGN as MAFFT_ALL } from '../modules/nf-core/mafft/align/main.nf'
include { MAFFT_ALIGN as MAFFT_ADD } from '../modules/nf-core/mafft/align/main.nf'
include { PATCHVCF                 } from '../modules/local/patchvcf/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCFALIGNER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    ch_samplesheet
        .map{ meta, fasta, vcf -> [[id:fasta] + meta.subMap('group'),fasta] }
        .unique()
        .collectFile{ meta, fasta ->
            ["${meta.group}.fa", fasta]
        }
        .map { file ->
            def group = file.simpleName
            def id = "${group}_references"
            [[id:id, group:group], file]}
        .set{fastas_all}

    MAFFT_ALL( fastas_all,[[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)
        .fas
        .join(fastas_all)
        .multiMap{ meta, aligned, fasta_all  ->
            reference: [meta + [id: "${meta.group}_aligned"], aligned]
            add: [meta, fasta_all]
        }
        .set{ch_aligned}

    MAFFT_ADD( ch_aligned.reference, ch_aligned.add, [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)
    added_maps = MAFFT_ADD.out.map.map{ meta, map  -> [meta.group, meta, map] }
    ch_versions = MAFFT_ADD.out.versions

    ch_map_vcf_fasta =ch_samplesheet
        .map { meta, fasta, vcf -> [meta.group, meta, fasta, vcf] }
        .combine(added_maps, by: 0)
        .map{ _group,  meta, fasta, vcf, _meta2, map ->
            [meta, map, vcf, fasta]
        }

    // Custom script to update the vcf or tsv/csv file with the new map coordinates of the alignment
    PATCHVCF(ch_map_vcf_fasta)
    ch_versions = ch_versions.mix(PATCHVCF.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'vcfaligner_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
