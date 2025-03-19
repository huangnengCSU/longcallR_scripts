import sys
import gzip
from Bio import SeqIO

if len(sys.argv) != 4:
    sys.exit("Usage: python set_pseudo_baseq.py input.fastq output.fastq base_quality_value")

input_fastq = sys.argv[1]
output_fastq = sys.argv[2]
fake_quality = int(sys.argv[3])

if input_fastq.endswith(".gz"):
    input_handle = gzip.open(input_fastq, "rt")
else:
    input_handle = open(input_fastq, "r")
output_handle = open(output_fastq, "w")
with input_handle, output_handle:
    for record in SeqIO.parse(input_handle, "fastq"):
        record.letter_annotations["phred_quality"] = [fake_quality] * len(record.seq)
        SeqIO.write(record, output_handle, "fastq")
print(f"FASTQ file with modified quality scores saved to {output_fastq}")