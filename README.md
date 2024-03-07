## Introduction

**emc/cancermicro** is a bioinformatics pipeline that can be used to analyse bacterial reads obtained from RNA sequencing of human tumour material. It takes a samplesheet and FASTQ files as input, performs quality control, removes human reads and taxonomically classifies the detected bacterial reads. It produces a quality control report and a bacterial abundance matrix.

<!-- TODO nf-core:
   Check if output is still up to date.
-->


<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline
Check if tools are still up to date. -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read QC([`fastp`](https://github.com/OpenGene/fastp))
3. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
4. Initial host-depletion ([`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
5. Second host-depletion ([`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2))
6. Classification ([`kraken2`](https://ccb.jhu.edu/software/kraken2/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

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
nextflow run emc/cancermicro \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

emc/cancermicro was originally written by Birgit Rijvers.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed. Remove if not -->


## Citations

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline to CITATIONS.MD file -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

**The nf-core framework for community-curated bioinformatics pipelines.**
Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
