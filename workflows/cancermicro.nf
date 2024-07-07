/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                               } from '../modules/nf-core/fastp/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { BWAMEM2_MEM                                         } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX                                       } from '../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FASTQ                                      } from '../modules/nf-core/samtools/fastq/main'
include { KRAKEN2_KRAKEN2                                     } from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN                                     } from '../modules/nf-core/bracken/bracken/main'
include { KRAKENBIOM_COMBINED   as KRAKENBIOM_COM_KR          } from '../modules/local/krakenbiom/krakenbiom_com/main.nf'
include { KRAKENBIOM_COMBINED   as KRAKENBIOM_COM_BR          } from '../modules/local/krakenbiom/krakenbiom_com/main.nf'
include { paramsSummaryMap                                    } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_cancermicro_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if(params.fasta){
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).collect()
        .map{ it -> [[id:it[0].getSimpleName()], it[0]]}
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CANCERMICRO {

    take:
    // channel: samplesheet read in from --input
    ch_samplesheet

    main:

    // Define channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_reports = Channel.empty()
    ch_fasta = Channel.empty()
    ch_kraken2_db = Channel.empty()
    ch_bracken_index = Channel.empty()

    //
    // MODULE: Run Fastp
    //
    FASTP (
        ch_samplesheet,
        [],
        false,
        false
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_reports = ch_reports.mix(FASTP.out.json.collect{ meta, json -> json })
    ch_reports = ch_reports.mix(FASTP.out.html.collect{ meta, html -> html })

    // Make tuple for fasta input with id and path from parameters
    ch_fasta = Channel.value([[id:'input_genome_fasta'], params.fasta])

    // Check if BWA-MEM2 index is provided
    if (!params.bwamem2_index) {
        //
        // MODULE: Run BWA-MEM2 index on provided fasta
        //
        BWAMEM2_INDEX (
            ch_fasta
        )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        ch_bwamem2_index = BWAMEM2_INDEX.out.index
    }
    else {
        // Use provided BWA-MEM2 index
        ch_bwamem2_index = Channel.value([[id:'input_genome_index'], params.bwamem2_index])
    }
    
    //
    // MODULE: Run BWA-MEM2 (extra args for samtools passed in 'modules.config')
    //
    BWAMEM2_MEM (
        FASTP.out.reads,
        ch_bwamem2_index,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    //
    // MODULE: Run Samtools fastq
    //
    SAMTOOLS_FASTQ (
        BWAMEM2_MEM.out.bam,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    ch_kraken2_db = Channel.value([params.kraken2_db])

    //
    // MODULE: Run Kraken2 (Confidence set to --confidence 0.05 in modules.config)
    //
    KRAKEN2_KRAKEN2 (
        SAMTOOLS_FASTQ.out.fastq,
        ch_kraken2_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    ch_bracken_index = Channel.value([[id:'kr_db_br'], params.kraken2_db])

    //
    // MODULE: Run Bracken
    //
    BRACKEN_BRACKEN (
        KRAKEN2_KRAKEN2.out.report,
        ch_kraken2_db
    )
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions)

    // Create channel with kreports
    ch_kreports = KRAKEN2_KRAKEN2.out.report.map { it[1] }.toList()

    // Create channel with bracken kreports
    ch_br_kreports = BRACKEN_BRACKEN.out.txt.map { it[1] }.toList()
    
    // MODULE: Run Kraken-biom on Kraken2 kreports COMBINED
    KRAKENBIOM_COM_KR (
        ch_kreports,
        "Kraken2"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_COM_KR.out.versions)

    // MODULE: Run Kraken-biom on bracken kreports COMBINED
    KRAKENBIOM_COM_BR (
        ch_br_kreports,
        "Bracken"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_COM_BR.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
