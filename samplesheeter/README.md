# Samplesheeter
`Samplesheeter.py` is a Python script that creates a CSV samplesheet from a directory of (zipped) FASTQ files. 
The script groups the FASTQ files by sample and writes the sample information to the CSV samplesheet file. 
It supports both paired-end and single-end FASTQ files.

## Features
- Groups FASTQ files by sample
- Supports both paired-end and single-end FASTQ files
- Skips UMI files with `_R3` suffix
- Allows specifying the file extension of FASTQ files, so can also be used on compressed files

## Requirements
- Python 3.x
- `argparse` module (included in the Python standard library)
- `csv` module (included in the Python standard library)
- `os` module (included in the Python standard library)

## Usage
```bash
python samplesheeter.py --samplesheet <path_to_samplesheet> --fastq_dir <path_to_fastq_directory> [--paired] [--file_extension <file_extension>]
```
or

```bash
python samplesheeter.py -s <path_to_samplesheet> -f <path_to_fastq_directory> [-p] [-e <file_extension>]
```
### Command Line Arguments
- `--samplesheet`, `-s`: Path to the samplesheet to be created (required)
- `--fastq_dir`, `-f`: Path to the directory containing FASTQ files (required)
- `--paired`, `-p`: Flag to indicate if the FASTQ files are paired-end (optional, **default is single-end**)
- `--file_extension`, `-e`: File extension of FASTQ files (optional, **default is .fastq.gz**)

### Example
To create a samplesheet for paired-end FASTQ files in the directory "/path/to/fastq_dir" with the default file extension ".fastq.gz":

```bash
python samplesheeter.py --samplesheet samplesheet.csv --fastq_dir /path/to/fastq_dir --paired
```
or
```bash
python samplesheeter.py -s samplesheet.csv -f /path/to/fastq_dir -p
```

To create a samplesheet for single-end FASTQ files in the directory "/path/to/fastq_dir" with a custom file extension ".fq":

```bash
python samplesheeter.py --samplesheet samplesheet.csv --fastq_dir /path/to/fastq_dir --file_extension .fq
```
or
```bash
python samplesheeter.py -s samplesheet.csv -f /path/to/fastq_dir -e .fq
```

## Output
The script generates a CSV samplesheet with the following columns:

**sample**: Sample name

**fastq_1**: Path to the first FASTQ file (for single-end and paired-end)

**fastq_2**: Path to the second FASTQ file (for paired-end only)