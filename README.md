## Introduction

**emc/cancermicro** is a ***WORK IN PROGRESS*** bioinformatics pipeline that helps with analyzing the microbiome within RNA sequencing data, obtained from humans. As input it requires a samplesheet and paired-end FASTQ files, it performs quality control and trimming on the reads, filters out reads mapping to a human reference genome and taxonomically classifies the remaining reads. Additionally, AMR genes are detected with RGI. As output, you get a BIOM file with the classifications, multiple visualizations and a MultiQC report of the QC metrics and tools used.

<!-- TODO nf-core:
Update introduction as pipeline changes!
-->
The steps performed by the pipeline are shown in the figure below:

<!-- TODO nf-core: Update metrochart to final version!   -->
![Metrochart_pipeline_no_Q2](https://github.com/BirgitRijvers/emc-cancermicro/assets/126883391/a72e2de0-ac23-4bf6-a13b-ad0b1031c663)

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
And include the following:

1. Quality control ([`Fastp`](https://github.com/OpenGene/fastp))
2. Filter out reads mapping to GRCh38 ([`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2) and [`Samtools`](https://www.htslib.org/doc/samtools-view.html))
3. Summarize mapping statistics (`Samtools`[https://www.htslib.org/doc/samtools-flagstat.html])
4. Convert SAM file to FASTQ ([`Samtools`](https://www.htslib.org/doc/samtools-fasta.html))
5. Taxonomic classification ([`Kraken2`](https://github.com/DerrickWood/kraken2))
6. Re-estimation of abundances ([`Bracken`](https://github.com/jenniferlu717/Bracken))
7. Convert Kreport to BIOM ([`Kraken-biom`](https://github.com/smdabdoub/kraken-biom))
8. Generate report with quality metrics and used tools ([`MultiQC`](https://github.com/MultiQC/MultiQC))
9. Filter classification results based on provided contaminant lists ([`QIIME2`](https://qiime2.org/))
10. Visualize classifications ([`QIIME2`](https://qiime2.org/))
11. Calculate alpha and beta diversity metrics ([`QIIME2`](https://qiime2.org/))

The final design of the fully developed pipeline is visualized in the figure below:
![Metrochart_pipeline_with_MQC](https://github.com/user-attachments/assets/5df4639e-7ad2-48cc-a57e-c63af03df651)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_1, BR_PVP_0705_R1.fastq.gz, BR_PVP_0705_R2.fastq.gz
```

Each row represents a pair of fastq files.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run emc/cancermicro \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --fasta <path/to/reference_genome_fasta> \
   --kraken2_db <path/to/kraken2_database/directory/>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

emc/cancermicro was originally written by Birgit Rijvers.

We thank the following people for their assistance in the development of this pipeline:
- Willem de Koning

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Support
If you have questions or issues while using this pipeline, please create an issue in this repository.


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use emc/cancermicro for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
