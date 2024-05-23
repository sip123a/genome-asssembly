from Bio import SeqIO
import argparse
import subprocess
from string import Template

# command  parameters
parser = argparse.ArgumentParser(description="Generate softmasked genome file for masked hard_regions in genome.fa")
parser.add_argument("-i", "--input", help="Input hard genome fasta file", required=True)
parser.add_argument("-o", "--output", help="Output soft genome fasta file", required=True)
parser.add_argument("-r", "--ref", help="Reference genome fasta file", required=True)
args = parser.parse_args()

# Function to write a BED record
def write_bed_record(file, chrom, start, end, name=".", score="."):
    file.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\n")

# Open the input FASTA file and output BED file
bed_file_path = "hard_masked_regions.bed"
with open(bed_file_path, "w") as bed_file:
    for record in SeqIO.parse(args.input, "fasta"):
        chrom = record.id  # Chromosome or contig name
        seq = str(record.seq)

        mask_start = None
        mask_end = None

        for i, base in enumerate(seq):
            if base == "N":
                # "N" indicates masked regions
                if mask_start is None:
                    mask_start = i
                mask_end = i
            elif mask_start is not None:
                write_bed_record(bed_file, chrom, mask_start, mask_end + 1)
                mask_start = None
                mask_end = None

        # Check if a masked region ends at the end of the sequence
        if mask_start is not None:
            write_bed_record(bed_file, chrom, mask_start, mask_end + 1)

def run_command(command):
    """Run a shell command."""
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def main(bed_file_path, ref, output):
    # bedtools maskfasta
    bedtools_maskfasta_command = Template("""
    bedtools maskfasta -fi $ref -bed $bed -fo $output -soft
    """)
    run_command(bedtools_maskfasta_command.substitute(ref=ref, bed=bed_file_path, output=output).strip())

    # sed
    sed_command = Template("""
    sed -i 's/n/N/g' $output
    """)
    run_command(sed_command.substitute(output=output).strip())

if __name__ == "__main__":
    main(bed_file_path, args.ref, args.output)
