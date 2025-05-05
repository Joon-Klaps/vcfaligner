/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MAFFT_ALIGN as MAFFT_ALL      } from '../modules/nf-core/mafft/align/main.nf'
include { MAFFT_ALIGN as MAFFT_INDIV    } from '../modules/nf-core/mafft/align/main.nf'
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
        .map{it -> [[id:"all_fasta"], it]}
        .set{fastas_all}

    MAFFT_ALL( fastas_all,[[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)
        .fas
        .set{ch_aligned}

    MAFFT_INDIV( ch_aligned, ch_fasta, [[:],[]], [[:],[]], [[:],[]], [[:],[]], true)
    ch_versions = MAFFT_INDIV.out.versions

    ch_map_vcf = ch_samplesheet
        .map{ meta, fasta, vcf -> [[id:fasta], meta, vcf] }
        .join(MAFFT_INDIV.out.map)
        .map{ fasta_id, meta, vcf, map ->
            [meta, map, vcf]
        }

    // Custom script to update the vcf or tsv/csv file with the new map coordinates of the alignment


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
