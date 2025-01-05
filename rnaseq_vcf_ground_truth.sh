#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 -b <bam_file> -d <depth_file> -r <reference_file> -f <high_conf_bed> -a <annotation_file> -o <output_dir>"
    exit 1
}

# Initialize variables
bam_file=""
depth_file=""
reference_file=""
high_confidence_bed_file=""
annotation_file=""
output_dir=""

# Parse options using getopts
while getopts "b:d:r:f:a:o:" opt; do
    case "$opt" in
        b) bam_file=$OPTARG ;;
        d) depth_file=$OPTARG ;;
        r) reference_file=$OPTARG ;;
        f) high_confidence_bed_file=$OPTARG ;;
        a) annotation_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        *) usage ;;  # Display usage on invalid option
    esac
done

# Ensure output directory is provided
if [ -z "$output_dir" ]; then
    echo "Error: Output directory (-o) is required."
    usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Validate parameters
if [ -z "$depth_file" ]; then
    if [ -z "$bam_file" ] || [ -z "$reference_file" ]; then
        echo "Error: Either -d <depth_file> must be provided OR both -b <bam_file> and -r <reference_file> are required."
        usage
    fi
fi

if [ -z "$high_confidence_bed_file" ] || [ -z "$annotation_file" ]; then
    echo "Error: -f <high_conf_bed>, and -a <annotation_file> are required."
    usage
fi

# Generate depth file if not provided
if [ ! -f "$depth_file" ]; then
    depth_file="$output_dir/all_chr_depth.tsv"
    samtools depth -a -d 100000 "$bam_file" --reference "$reference_file" > "$depth_file"
fi

# Generate BED files for coverage ranges and store in output directory
for range in "10-20" "20-40" "40-100" "100-inf" "20-inf"; do
    awk_cmd=""
    case $range in
        "10-20") awk_cmd='(int($3) >= 10 && int($3) < 20)' ;;
        "20-40") awk_cmd='(int($3) >= 20 && int($3) < 40)' ;;
        "40-100") awk_cmd='(int($3) >= 40 && int($3) < 100)' ;;
        "100-inf") awk_cmd='(int($3) >= 100)' ;;
        "20-inf") awk_cmd='(int($3) >= 20)' ;;
    esac

    awk -F '\t' "$awk_cmd {printf(\"%s\t%d\t%s\n\", \$1, int(\$2) - 1, \$2);}" "$depth_file" | \
        sort -t $'\t' -k1,1 -k2,2n | bedtools merge > "$output_dir/all_chr.filter.${range}x.bed"
done

# Intersect with high confidence regions
for range in "10-20" "20-40" "40-100" "100-inf" "20-inf"; do
    bedtools intersect -a "$output_dir/all_chr.filter.${range}x.bed" -b "$high_confidence_bed_file" > \
        "$output_dir/all_chr.filter.${range}x.high_conf.bed"
done

# Convert annotation to BED and extract exons
gtf2bed < "$annotation_file" | grep -wF exon > "$output_dir/exons.bed"

# Intersect high confidence regions with exons
for range in "10-20" "20-40" "40-100" "100-inf" "20-inf"; do
    bedtools intersect -a "$output_dir/all_chr.filter.${range}x.high_conf.bed" -b "$output_dir/exons.bed" > \
        "$output_dir/all_chr.filter.${range}x.high_conf.exons.bed"
done

# Subtract exons from high confidence regions
for range in "10-20" "20-40" "40-100" "100-inf" "20-inf"; do
    bedtools subtract -a "$output_dir/all_chr.filter.${range}x.high_conf.bed" -b "$output_dir/exons.bed" > \
        "$output_dir/all_chr.filter.${range}x.high_conf.non_exons.bed"
done

# Display completion message
echo "Processing completed. All output files are stored in: $output_dir"
