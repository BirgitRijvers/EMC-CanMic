import os
import csv
import argparse

def create_samplesheet():
    """
    Create a CSV samplesheet from a directory of (zipped) FASTQ files.

    This script takes in command line arguments to specify the path to the samplesheet to be created,
    the path to the directory containing the FASTQ files, and whether the FASTQ files are paired-end or single-end.

    The script groups the FASTQ files by sample and writes the sample information to the CSV samplesheet file.

    Command line arguments:
    --samplesheet, -s: Path to samplesheet to be created
    --fastq_dir, -f: Path to directory containing FASTQ files
    --paired, -p: Flag to indicate if the FASTQ files are paired-end
    --file_extension, -e: File extension of FASTQ files
    """

    # Create parser
    parser = argparse.ArgumentParser(
        description=('Create a CSV samplesheet from a directory of FASTQ files.')
    )
    # Add arguments
    parser.add_argument('--samplesheet', '-s', type=str, required=True,
                        help='Path to samplesheet to be created')
    parser.add_argument('--fastq_dir', '-f', type=str, required=True,
                        help='Path to directory containing FASTQ files')
    parser.add_argument('--paired', '-p', action='store_true', default=False,
                        help='Flag to indicate if the FASTQ files are paired-end')
    parser.add_argument('-file_extension', '-e', type=str, default=".fastq.gz",
                        help='File extension of FASTQ files')
    # Parse!
    args = parser.parse_args()

    # Directory containing the FASTQ files
    data_dir = args.fastq_dir

    # Output CSV file
    output_file = args.samplesheet

    # Paired or single
    paired = args.paired

    # File extension
    file_extension = args.file_extension

    # Get a list of all FASTQ files in the directory
    fastq_files = [f for f in os.listdir(data_dir) if f.endswith(f"{file_extension}")]
    # Group the FASTQ files by sample
    samples = {}

    # If the FASTQ files are paired-end
    if paired:
        # Loop through files and group by sample
        for filename in fastq_files:
            parts = filename.split("_")
            sample = "_".join(parts[:-1])  # Join all parts except the last one (R1 or R2)
            if sample not in samples:
                samples[sample] = []
            samples[sample].append(os.path.abspath(os.path.join(data_dir, filename)))

        # Write the samples to the CSV file
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            # Write header
            writer.writerow(["sample", "fastq_1", "fastq_2"])
            # Loop through samples and write name + paths to file
            for sample, filenames in samples.items():
                # Ensure R1 is before R2
                filenames.sort()
                row = [sample]
                row.extend(filenames)
                writer.writerow(row)
    
    # If the FASTQ files are single-end
    else:
        # Loop through files and group by sample
        for filename in fastq_files:
            sample = filename.split(".")[0]
            if sample not in samples:
                samples[sample] = []
            samples[sample].append(os.path.abspath(os.path.join(data_dir, filename)))

        # Write the samples to the CSV file
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            # Write header
            writer.writerow(["sample", "fastq_1"])
            for sample, filenames in samples.items():
                row = [sample]
                row.extend(filenames)
                writer.writerow(row)

# Call the function to create the samplesheet
if __name__ == "__main__":
    create_samplesheet()
