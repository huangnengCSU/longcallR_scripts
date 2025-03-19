import sys
import gzip
from Bio import SeqIO

# Check if the input file is provided
if len(sys.argv) != 4:
    sys.exit("Usage: python modify_fake_baseq.py input.fastq output.fastq fake_quality")

# Input and output file paths
input_fastq = sys.argv[1]
output_fastq = sys.argv[2]
fake_quality = int(sys.argv[3])

# Open input file (handle compressed files)
if input_fastq.endswith(".gz"):
    input_handle = gzip.open(input_fastq, "rt")
else:
    input_handle = open(input_fastq, "r")

# Open output file
output_handle = open(output_fastq, "w")

# Read the FASTQ file and modify the quality scores
with input_handle, output_handle:
    for record in SeqIO.parse(input_handle, "fastq"):
        # Set fake quality scores
        record.letter_annotations["phred_quality"] = [fake_quality] * len(record.seq)
        # Write the modified record to the output file
        SeqIO.write(record, output_handle, "fastq")

print(f"FASTQ file with modified quality scores saved to {output_fastq}")