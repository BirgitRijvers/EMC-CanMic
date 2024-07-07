// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { QIIME2_IMPORT                                       } from '../../modules/local/qiime2/import/main'
include { QIIME2_SUMMARIZE                                    } from '../../modules/local/qiime2/summarize/main'
include { QIIME2_BARPLOT as QIIME2_BARPLOT_RAW                } from '../../modules/local/qiime2/barplot/main'
include { QIIME2_BARPLOT as QIIME2_BARPLOT_FILTERED_WHITE     } from '../../modules/local/qiime2/barplot/main'
include { QIIME2_BARPLOT as QIIME2_BARPLOT_FILTERED_BLACK     } from '../../modules/local/qiime2/barplot/main'
include { QIIME2_FILTER_WHITELIST                             } from '../../modules/local/qiime2/decontaminate/whitelist/main'
include { QIIME2_FILTER_BLACKLIST                             } from '../../modules/local/qiime2/decontaminate/blacklist/main'
include { QIIME2_COLLAPSE as QIIME2_COLLAPSE_RAW              } from '../../modules/local/qiime2/collapse/main'
include { QIIME2_COLLAPSE as QIIME2_COLLAPSE_WHITE            } from '../../modules/local/qiime2/collapse/main'
include { QIIME2_COLLAPSE as QIIME2_COLLAPSE_BLACK            } from '../../modules/local/qiime2/collapse/main'
include { QIIME2_HEATMAP as QIIME2_HEATMAP_RAW                } from '../../modules/local/qiime2/heatmap/main'
include { QIIME2_HEATMAP as QIIME2_HEATMAP_WHITE              } from '../../modules/local/qiime2/heatmap/main'
include { QIIME2_HEATMAP as QIIME2_HEATMAP_BLACK              } from '../../modules/local/qiime2/heatmap/main'
include { QIIME2_COREDIVERSITY as QIIME2_COREDIVERSITY_RAW    } from '../../modules/local/qiime2/corediversity/main'
include { QIIME2_COREDIVERSITY as QIIME2_COREDIVERSITY_WHITE  } from '../../modules/local/qiime2/corediversity/main'
include { QIIME2_COREDIVERSITY as QIIME2_COREDIVERSITY_BLACK  } from '../../modules/local/qiime2/corediversity/main'

workflow QIIME2 {

    take:
    // TODO nf-core: edit input (take) channels
    // ch_bam // channel: [ val(meta), [ bam ] ]
    KRAKENBIOM_COM_KR.out.biom

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    // SAMTOOLS_SORT ( ch_bam )
    // ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Run QIIME2 import
    //
    QIIME2_IMPORT (
        KRAKENBIOM_COM_KR.out.biom
    )
    ch_versions = ch_versions.mix(QIIME2_IMPORT.out.versions.first())

    //
    // MODULE: Run QIIME2 feature-table summarize on raw data
    //
    QIIME2_SUMMARIZE (
        QIIME2_IMPORT.out.frequency
    )
    ch_versions = ch_versions.mix(QIIME2_SUMMARIZE.out.versions.first())

    //
    // MODULE: Run QIIME2 barplot on raw data
    //
    QIIME2_BARPLOT_RAW (
        QIIME2_IMPORT.out.frequency,
        QIIME2_IMPORT.out.taxonomy,
        "raw"
    )
    ch_versions = ch_versions.mix(QIIME2_BARPLOT_RAW.out.versions.first())

    ch_whitelist = Channel.fromPath(params.whitelist).map { it.text.trim() }

    //
    // MODULE: Run QIIME2 filter whitelist
    //
    QIIME2_FILTER_WHITELIST (
        QIIME2_IMPORT.out.frequency,
        QIIME2_IMPORT.out.taxonomy,
        ch_whitelist
    )
    ch_versions = ch_versions.mix(QIIME2_FILTER_WHITELIST.out.versions.first())

    ch_blacklist = Channel.fromPath(params.blacklist).map { it.text.trim() }

    //
    // MODULE: Run QIIME2 filter blacklist
    //
    QIIME2_FILTER_BLACKLIST (
        QIIME2_IMPORT.out.frequency,
        QIIME2_IMPORT.out.taxonomy,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(QIIME2_FILTER_BLACKLIST.out.versions.first())
    
    //
    // MODULE: Run QIIME2 barplot on whitelist filtered data
    //
    QIIME2_BARPLOT_FILTERED_WHITE (
        QIIME2_FILTER_WHITELIST.out.filtered_table,
        QIIME2_IMPORT.out.taxonomy,
        "filtered_white"
    )
    ch_versions = ch_versions.mix(QIIME2_BARPLOT_FILTERED_WHITE.out.versions.first())

    //
    // MODULE: Run QIIME2 barplot on blacklist filtered data
    //
    QIIME2_BARPLOT_FILTERED_BLACK (
        QIIME2_FILTER_BLACKLIST.out.filtered_table,
        QIIME2_IMPORT.out.taxonomy,
        "filtered_black"
    )
    ch_versions = ch_versions.mix(QIIME2_BARPLOT_FILTERED_BLACK.out.versions.first())

    // 
    // MODULE: Run QIIME2 collapse on raw data
    //
    QIIME2_COLLAPSE_RAW (
        QIIME2_IMPORT.out.frequency,
        QIIME2_IMPORT.out.taxonomy
    )
    ch_versions = ch_versions.mix(QIIME2_COLLAPSE_RAW.out.versions.first())
    
    // 
    // MODULE: Run QIIME2 collapse on whitelist filtered data
    //
    QIIME2_COLLAPSE_WHITE (
        QIIME2_FILTER_WHITELIST.out.filtered_table,
        QIIME2_IMPORT.out.taxonomy,
    )
    ch_versions = ch_versions.mix(QIIME2_COLLAPSE_WHITE.out.versions.first())

    // 
    // MODULE: Run QIIME2 collapse on blacklist filtered data
    //
    QIIME2_COLLAPSE_BLACK (
        QIIME2_FILTER_WHITELIST.out.filtered_table,
        QIIME2_IMPORT.out.taxonomy,
    )
    ch_versions = ch_versions.mix(QIIME2_COLLAPSE_BLACK.out.versions.first())

    //
    // MODULE: Run QIIME2 core diversity on raw data
    //
    QIIME2_COREDIVERSITY_RAW (
        QIIME2_COLLAPSE_RAW.out.collapsed_table,
        params.sampling_depth,
        params.metadata,
        "raw"
    )
    ch_versions = ch_versions.mix(QIIME2_COREDIVERSITY_RAW.out.versions.first())

    //
    // MODULE: Run QIIME2 core diversity on whitelist filtered data
    //
    QIIME2_COREDIVERSITY_WHITE (
        QIIME2_COLLAPSE_WHITE.out.collapsed_table,
        params.sampling_depth,
        params.metadata,
        "white"
    )
    ch_versions = ch_versions.mix(QIIME2_COREDIVERSITY_WHITE.out.versions.first())

    //
    // MODULE: Run QIIME2 core diversity on blacklist filtered data
    //
    QIIME2_COREDIVERSITY_BLACK (
        QIIME2_COLLAPSE_BLACK.out.collapsed_table,
        params.sampling_depth,
        params.metadata,
        "black"
    )
    ch_versions = ch_versions.mix(QIIME2_COREDIVERSITY_BLACK.out.versions.first())

    //
    // MODULE: Run QIIME2 heatmap on raw data
    //
    QIIME2_HEATMAP_RAW (
        QIIME2_COLLAPSE_RAW.out.collapsed_table,
        QIIME2_IMPORT.out.taxonomy,
        "raw"
    )
    ch_versions = ch_versions.mix(QIIME2_HEATMAP_RAW.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

