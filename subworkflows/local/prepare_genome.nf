//
// Index the reference genome and align reads
//

include { BWAMEM2_INDEX      } from '../modules/nf-core/bwamem2/index/main'

workflow PREPARE_GENOME {
    take:
    fasta         // [[id: 'identifier'], path(genome)]

    main:

    ch_versions = Channel.empty()
    ch_bwamem2_index = Channel.empty()
    ch_fasta = fasta

    if (!params.bwamem2_index){
        BWAMEM2_INDEX ( ch_fasta )
        ch_bwamem2_index = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else{
        ch_bwamem2_index = Channel.of([[id: 'input_genome_index'],params.bwamem2_index])
    }


    emit:
    fasta         = ch_fasta                // [val(meta), path(genome.fasta)]
    bwamem2_index = ch_bwamem2_index        // [val(meta), index_data ]
    versions      = ch_versions             // channel: [ versions.yml ]
}