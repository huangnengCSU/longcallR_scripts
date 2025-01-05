import gzip
import argparse


def load_hap_result(hap_vcf):
    """
    Parse the haplotype VCF file to extract true positives (TP), false positives (FP), and false negatives (FN).

    Args:
        hap_vcf (str): Path to the VCF file (can be .gz compressed).

    Returns:
        tuple: Sets of TP, FP, and FN SNP identifiers (chrom:pos).
    """

    def open_file(filepath):
        return gzip.open(filepath, 'rb') if filepath.endswith(".gz") else open(filepath, 'r')

    tp_snps, fp_snps, fn_snps = set(), set(), set()
    tp_cnt, fp_cnt, fn_cnt = 0, 0, 0

    with open_file(hap_vcf) as f:
        for line in f:
            if hap_vcf.endswith(".gz"):
                line = line.decode('utf-8')
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            truth, query = fields[9], fields[10]

            snp_id = f"{chrom}:{pos}"
            if "TP" in truth and "SNP" in truth and "TP" in query:
                tp_cnt += 1
                tp_snps.add(snp_id)
            else:
                if "FN" in truth and "SNP" in truth:
                    fn_cnt += 1
                    fn_snps.add(snp_id)
                if "FP" in query and "SNP" in query:
                    fp_cnt += 1
                    fp_snps.add(snp_id)

    return tp_snps, fp_snps, fn_snps


def write_to_file(query_vcf_file, output_prefix, tp_snps, fp_snps, fn_snps):
    """
    Write TP, FP, and FN variants from the query VCF file into separate output files.

    Args:
        query_vcf_file (str): Path to the query VCF file (can be .gz compressed).
        output_prefix (str): Prefix for the output files.
        tp_snps (set): Set of true positive SNP identifiers (chrom:pos).
        fp_snps (set): Set of false positive SNP identifiers (chrom:pos).
        fn_snps (set): Set of false negative SNP identifiers (chrom:pos).
    """

    def open_file(filepath):
        return gzip.open(filepath, 'rb') if filepath.endswith(".gz") else open(filepath, 'r')

    def write_to_outputs(header, fp_output, fn_output, tp_output):
        fp_output.write(header)
        fn_output.write(header)
        tp_output.write(header)

    out_fp = open(f"{output_prefix}.fp.vcf", 'w')
    out_fn = open(f"{output_prefix}.fn.vcf", 'w')
    out_tp = open(f"{output_prefix}.tp.vcf", 'w')

    with open_file(query_vcf_file) as f:
        for line in f:
            if query_vcf_file.endswith('.gz'):
                line = line.decode('utf-8')
            if line.startswith('#'):
                write_to_outputs(line, out_fp, out_fn, out_tp)
                continue

            fields = line.strip().split('\t')
            chrom, pos = fields[0], fields[1]
            snp_id = f"{chrom}:{pos}"

            if snp_id in tp_snps:
                out_tp.write(line)
            elif snp_id in fp_snps:
                out_fp.write(line)
            elif snp_id in fn_snps:
                out_fn.write(line)

    out_fp.close()
    out_fn.close()
    out_tp.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract TP, FP, and FN variants from a query VCF file.")
    parser.add_argument("-hap", "--hap_vcf", required=True, help="Path to the haplotype VCF file.")
    parser.add_argument("-query", "--query_vcf", required=True, help="Path to the query VCF file.")
    parser.add_argument("-out", "--output_prefix", required=True, help="Prefix for the output files.")
    args = parser.parse_args()

    tp_snps, fp_snps, fn_snps = load_hap_result(args.hap_vcf)
    write_to_file(args.query_vcf, args.output_prefix, tp_snps, fp_snps, fn_snps)