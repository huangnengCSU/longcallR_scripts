#!/bin/bash

# Set up arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -i)
      input_phased_vcf="$2"
      shift
      shift
      ;;
    -t)
      truth_vcf="$2"
      shift
      shift
      ;;
    -o)
      output_prefix="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check if required arguments are provided
if [[ -z "$input_phased_vcf" || -z "$truth_vcf" || -z "$output_prefix" ]]; then
  echo "Usage: $0 -i <input_phased_vcf> -t <truth_vcf> -o <output_prefix>"
  exit 1
fi

# Check if input VCF files are uncompressed and end with .vcf
if [[ ! "$input_phased_vcf" =~ \.vcf$ ]]; then
  echo "Error: input_phased_vcf must be an uncompressed VCF file ending with .vcf"
  exit 1
fi
if [[ ! "$truth_vcf" =~ \.vcf$ ]]; then
  echo "Error: truth_vcf must be an uncompressed VCF file ending with .vcf"
  exit 1
fi

# Run Whatshap stats
whatshap stats --gtf=${output_prefix}.gtf ${input_phased_vcf}

# Run Whatshap compare
whatshap compare --ignore-sample-name --tsv-pairwise=${output_prefix}.tsv ${truth_vcf} ${input_phased_vcf}
