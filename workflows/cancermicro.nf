/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                                               } from '../modules/nf-core/fastp/main'
include { FASTQC                                              } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { BWAMEM2_MEM                                         } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX                                       } from '../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FASTQ                                      } from '../modules/local/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_GZIP               } from '../modules/nf-core/samtools/fastq/main'
// include { KRAKEN2_SBUILD                                      } from '../modules/local/kraken2/sbuild/main'
include { KRAKEN2_KRAKEN2                                     } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA                           } from '../modules/nf-core/krakentools/kreport2krona/main'
// include { BRACKEN_BUILD                                       } from '../modules/nf-core/bracken/build/main'
include { BRACKEN_BRACKEN                                     } from '../modules/nf-core/bracken/bracken/main'
include { FARGENE                                             } from '../modules/nf-core/fargene/main'
include { UNTAR                                               } from '../modules/nf-core/untar/main'
include { SEQKIT_FQ2FA                                        } from '../modules/nf-core/seqkit/fq2fa/main'
include { RGI_CARDANNOTATION                                  } from '../modules/nf-core/rgi/cardannotation/main.nf'
include { RGI_MAIN                                            } from '../modules/nf-core/rgi/main/main.nf'
include { KRAKENBIOM_INDIVIDUAL as KRAKENBIOM_IND_KR          } from '../modules/local/krakenbiom/krakenbiom_ind/main.nf'
include { KRAKENBIOM_COMBINED   as KRAKENBIOM_COM_KR          } from '../modules/local/krakenbiom/krakenbiom_com/main.nf'
include { KRAKENBIOM_INDIVIDUAL as KRAKENBIOM_IND_BR          } from '../modules/local/krakenbiom/krakenbiom_ind/main.nf'
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
include { QIIME2                                              } from '../subworkflows/local/qiime2'

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
    // MODULE: Run FastQC raw data
    //
    FASTQC (
        ch_samplesheet
    )  
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

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
    // ch_reports = ch_reports.mix(FASTP.out.json.collect{ meta, json -> json })
    // ch_reports = ch_reports.mix(FASTP.out.html.collect{ meta, html -> html })

    //
    // MODULE: Run Seqkit FQ2FA
    //
    SEQKIT_FQ2FA (
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix(SEQKIT_FQ2FA.out.versions)

    // 
    // MODULE: Download and untar CARD database
    //
    UNTAR( 
        [ [], file('https://card.mcmaster.ca/latest/data', checkIfExists: true) ] 
    )
    ch_versions = ch_versions.mix( UNTAR.out.versions )
    rgi_db = UNTAR.out.untar.map{ it[1] }

    // 
    // MODULE: Run RGI Card Annotation
    //
    RGI_CARDANNOTATION (
        rgi_db
    )
    ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)

    //
    // MODULE: Run RGI Main
    //
    RGI_MAIN (
        SEQKIT_FQ2FA.out.fasta,
        RGI_CARDANNOTATION.out.db,
        []
    )
    ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

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
    
    ch_bam_copy1 = BWAMEM2_MEM.out.bam.map { it }
    ch_bam_copy2 = BWAMEM2_MEM.out.bam.map { it }

    //
    // MODULE: Run Samtools fastq GZIP
    //
    SAMTOOLS_FASTQ_GZIP (
        ch_bam_copy1,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_GZIP.out.versions)

    //
    // MODULE: Run Samtools fastq uncompressed
    //
    SAMTOOLS_FASTQ (
        ch_bam_copy2,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    // copy of SAMTOOLS_FASTQ.out.fastq for fastas:
    fastas = SAMTOOLS_FASTQ.out.fastq 

    // fARGene run
    if ( !params.skip_fargene ) {
        ch_fargene_classes = Channel.fromList( params.fargene_hmmmodel.tokenize(',') )

        ch_fargene_input = fastas
                            .combine( ch_fargene_classes )
                            .map {
                                meta, fastas, hmm_class ->
                                    def meta_new = meta.clone()
                                    meta_new['hmm_class'] = hmm_class
                                [ meta_new, fastas, hmm_class ]
                            }
                            .multiMap {
                                fastas: [ it[0], it[1] ]
                                hmmclass: it[2]
                            }

        FARGENE (
            ch_fargene_input.fastas, 
            ch_fargene_input.hmmclass
        )
        ch_versions = ch_versions.mix(FARGENE.out.versions)
    }

    // KRAKEN2_SBUILD does not yet work as expected so commented out for now

    // ch_kr_db = Channel.empty()
    // ch_kr_db = Channel.value([[id:'kraken2_standard_db'], '/this_is_a_name'])

    // // Check if Kraken2 database is provided
    // if (!params.kraken2_db) {
    //     //
    //     // MODULE: Build Kraken2 standard database
    //     //
    //     KRAKEN2_SBUILD (
    //         ch_kr_db,
    //         'standard_tester_db',
    //         true
    //     )
    //     ch_versions = ch_versions.mix(KRAKEN2_SBUILD.out.versions.first())
    //     ch_kraken2_db = KRAKEN2_SBUILD.out.database
    //     }
    //     else {
    //         // Use provided Kraken2 database
    //         ch_kraken2_db = Channel.value([params.kraken2_db])
    //     }

    ch_kraken2_db = Channel.value([params.kraken2_db])

    //
    // MODULE: Run Kraken2 (Confidence set to --confidence 0.05 in modules.config)
    //
    KRAKEN2_KRAKEN2 (
        SAMTOOLS_FASTQ_GZIP.out.fastq,
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

    // Bracken_build commented out as long as Kraken2_SBUILD does not work

    // //
    // // MODULE: Build Bracken database
    // //
    // BRACKEN_BUILD (
    //     ch_bracken_index
    // )

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
        "kraken2"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_COM_KR.out.versions)

    // MODULE: Run Kraken-biom on bracken kreports COMBINED
    KRAKENBIOM_COM_BR (
        ch_br_kreports,
        "bracken"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_COM_BR.out.versions)

    //
    // MODULE: Run Kraken-biom on Kraken2 kreports SEPARATELY
    //
    KRAKENBIOM_IND_KR (
        KRAKEN2_KRAKEN2.out.report,
        "kraken2"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_IND_KR.out.versions)
    
    // MODULE: Run Kraken-biom on Bracken kreports SEPARATELY
    KRAKENBIOM_IND_BR (
        BRACKEN_BRACKEN.out.txt,
        "bracken"
    )
    ch_versions = ch_versions.mix(KRAKENBIOM_IND_BR.out.versions)

    // Run subworkflow if specified
    if (params.QIIME2) {
    QIIME2 (
        KRAKENBIOM_COM_KR.out.biom
    )
    ch_versions = ch_versions.mix(QIIME2.out.versions)
    }

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
