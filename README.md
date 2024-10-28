## Introduction

**EMC/metamicrobes** is a bioinformatics pipeline that analyzes microbial signatures and AMR genes in metagenomic or metatranscriptomic data. 

As input it requires a samplesheet with paths to paired-end short-read (compressed) FASTQ files, it performs quality control and trimming on the reads, filters out reads mapping to a specified host reference genome and taxonomically classifies the remaining reads. In addition it also detects antimicrobial resistance genes with tow different approaches. As output, you receive all intermediate outputs as well as a BIOM file with the classifications and a MultiQC report of the QC metrics and tools used.

An overview of the steps implemented in MetaMicrobes are shown in the figure below:

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

And include the following:

1. Quality control ([`Fastp`](https://github.com/OpenGene/fastp))
2. Filter out reads mapping to a reference genome ([`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2) and [`Samtools`](https://www.htslib.org/doc/samtools-view.html))
3. Summarize mapping statistics ([`Samtools`](https://www.htslib.org/doc/samtools-flagstat.html))
4. Convert SAM file to FASTQ ([`Samtools`](https://www.htslib.org/doc/samtools-fasta.html))
5. Taxonomic classification ([`Kraken2`](https://github.com/DerrickWood/kraken2))
6. Visualize Kraken2 output with Krona ([`KrakenTools`](https://github.com/jenniferlu717/KrakenTools))
7. Re-estimation of microbial abundances ([`Bracken`](https://github.com/jenniferlu717/Bracken))
8. Convert Kreport to BIOM ([`Kraken-biom`](https://github.com/smdabdoub/kraken-biom))
9. Decontaminate based on a blacklist and whitelist ([`QIIME2`](https://qiime2.org/))
10. Visualize microbial profiles with barcharts and heatmaps ([`QIIME2`](https://qiime2.org/))
11. Assess microbial alpha and beta diversity ([`QIIME2`](https://qiime2.org/))
12. Generate report with quality metrics and used tools ([`MultiQC`](https://github.com/MultiQC/MultiQC))

## Usage
To use MetaMicrobes on your machine, follow the steps below:
1. Make sure you have correctly set-up Nextflow and it's dependencies
> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

2. Clone this GitHub repository
3. [Test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.
4. Prepare a samplesheet like the example below:
    `samplesheet.csv`:
    ```csv
    sample,fastq_1,fastq_2
    CONTROL_1, BR_PVP_0705_R1.fastq.gz, BR_PVP_0705_R2.fastq.gz
    ```
    Each row represents a pair of fastq files.
5. Now, you can run the MetaMicrobes pipeline using:
    ```bash
    nextflow run emc-cancermicro \
       -profile <docker/singularity/conda/.../institute> \
       --input samplesheet.csv \
       --outdir <OUTDIR> \
       --fasta <path/to/reference_genome_fasta> \
       --kraken2_db <path/to/kraken2_database/directory/>
    ```

    Change the default "null" values in "nextflow.config" to the paths you will use most often to save time. Values in this file will be overwritten by the values specified in the command.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

6.  

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run EMC/metamicrobes \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

EMC/metamicrobes was originally written by Birgit Rijvers.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use EMC/metamicrobes for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
