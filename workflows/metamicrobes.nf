/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                  } from '../modules/nf-core/fastp/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BWAMEM2_MEM            } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX          } from '../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FASTQ         } from '../modules/nf-core/samtools/fastq/main'
include { KRAKEN2_BUILDSTANDARD  } from '../modules/nf-core/kraken2/buildstandard/main'
include { KRAKEN2_KRAKEN2        } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA              } from '../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENBIOM_KRAKENBIOM as KRAKENBIOM_KR } from '../modules/local/krakenbiom/main'
include { KRAKENBIOM_KRAKENBIOM as KRAKENBIOM_BR } from '../modules/local/krakenbiom/main'
include { BRACKEN_BUILD          } from '../modules/nf-core/bracken/build/main'
include { BRACKEN_BRACKEN        } from '../modules/nf-core/bracken/bracken/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_metamicrobes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QIIME2 as QIIME2_KRAKEN2                            } from '../subworkflows/local/qiime2'
include { QIIME2 as QIIME2_BRACKEN                            } from '../subworkflows/local/qiime2'

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

workflow METAMICROBES {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Fastp on raw reads
    //
    FASTP (
        ch_samplesheet,
        [],
        false,
        false,
        false
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions)

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
    // MODULE: Run BWA-MEM2
    //
    // Extra args for BWA-MEM2 are passed in 'modules.config' to save only unmapped readpairs
    BWAMEM2_MEM (
        FASTP.out.reads,
        ch_bwamem2_index,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    //
    // MODULE: Run Samtools fastq on BWA-MEM2 output
    //
    SAMTOOLS_FASTQ (
        BWAMEM2_MEM.out.bam,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    
    // Check if Kraken2 database is provided
    if (!params.kraken2_db) {
        //
        // MODULE: Run Kraken2 build standard database
        //
        KRAKEN2_BUILDSTANDARD (
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2_BUILDSTANDARD.out.versions)
        // Add built Kraken2 database to channel
        ch_kraken2_db = KRAKEN2_BUILDSTANDARD.out.db
    } else {
            // Use provided Kraken2 database
            ch_kraken2_db = Channel.value([params.kraken2_db])
    }

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
    
    //
    // MODULE: Run KrakenTools kreport2krona
    //
    KRAKENTOOLS_KREPORT2KRONA (
        KRAKEN2_KRAKEN2.out.report
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)

    // Create channel with Kraken2 reports list
    ch_kreports = KRAKEN2_KRAKEN2.out.report.map {it[1]}.toList()

// Check if Bracken database is provided
    if (!params.bracken_db) {
        // Create channel with built Kraken2 database and meta id
        ch_kraken2_db_path = ch_kraken2_db.map{it}
        ch_bracken_index = ch_kraken2_db_path.map { db_path -> tuple([id: 'kraken2_db_for_bracken'], db_path) }
        //
        // MODULE: Run Bracken build 
        //
        BRACKEN_BUILD (
            ch_bracken_index
        )
        ch_versions = ch_versions.mix(BRACKEN_BUILD.out.versions)
        // Extract second item from tuple channel
        ch_bracken_db = BRACKEN_BUILD.out.db.map {it[1]}
    } else {
            // Use provided Bracken database
            ch_bracken_db = Channel.value([params.bracken_db])
    }
    //
    // MODULE: Run Bracken
    //
    BRACKEN_BRACKEN (
        KRAKEN2_KRAKEN2.out.report,
        ch_bracken_db
    )
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions)
    
    // Create channel with Bracken reports list
    ch_br_kreports = BRACKEN_BRACKEN.out.txt.map { it[1] }.toList()

    //
    // MODULE: Run Kraken-biom on Kraken2 reports
    //
    KRAKENBIOM_KR (
        ch_kreports,
        "kraken2"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_KR.out.versions)

    //
    // MODULE: Run Kraken-biom on Bracken reports
    //
    KRAKENBIOM_BR (
        ch_br_kreports,
        "bracken"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_BR.out.versions)

    // Run Q2 subworkflow on Kraken2 and Bracken outputs if specified
    if (params.qiime2) {
        QIIME2_KRAKEN2 (
            KRAKENBIOM_KR.out.biom
        )
        ch_versions = ch_versions.mix(QIIME2_KRAKEN2.out.versions)
        QIIME2_BRACKEN (
            KRAKENBIOM_BR.out.biom
        )
        ch_versions = ch_versions.mix(QIIME2_BRACKEN.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

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
