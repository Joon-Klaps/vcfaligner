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
        .map{ meta, fasta, vcf -> [[id:fasta],fasta] }
        .unique()
        .tap{ch_fasta}
        .collectFile(name: "all_fasta.fa"){it[1]}
        .map{it -> [[id:"references"], it]}
        .set{fastas_all}

    MAFFT_ALL( fastas_all,[[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)
        .fas
        .map{ meta, fasta  -> [[id:"ref_aligned"], fasta] }
        .set{ch_aligned}

    MAFFT_ADD( ch_aligned, fastas_all, [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)
    ch_versions = MAFFT_ADD.out.versions

    ch_map_vcf_fasta = ch_samplesheet
        .combine(MAFFT_ADD.out.map)
        .map{ meta, fasta, vcf, mafft_id, map ->
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
