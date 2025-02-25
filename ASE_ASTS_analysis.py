import time
import gzip
import pysam
import numpy as np
from scipy.stats import binomtest
from scipy.stats import chi2_contingency, chi2
from statsmodels.stats.multitest import multipletests
import concurrent.futures
from multiprocessing import Manager
import argparse


def parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type):
    start_time = time.time()  # Start time measurement
    records = {}
    with open(assignment_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[5] in assignment_type:
                read_name, isoform_id, gene_id, additional_info = parts[0], parts[3], parts[4], parts[8]
                additional_info_dict = {key.strip(): value.strip() for key, value in
                                        (pair.split('=') for pair in additional_info.split(';') if pair.strip())}
                try:
                    classification = additional_info_dict["Classification"]
                    if classification in classification_type:
                        records.setdefault(gene_id, {}).setdefault(isoform_id, []).append(read_name)
                except KeyError:
                    continue
    end_time = time.time()  # End time measurement
    elapsed_time = end_time - start_time
    print(f"parse_isoquant_read_assignment runtime: {elapsed_time:.2f} seconds")
    return records


def get_gene_regions(annotation_file, gene_types):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.
    :param annotation_file: Path to the annotation file
    :return: Gene regions, exon regions, and intron regions
    """
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unsupported annotation format"

    gene_regions = {}
    gene_names = {}

    def process_gene(parts, gene_id, gene_name):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}  # 1-based, start-inclusive, end-inclusive
        gene_names[gene_id] = gene_name

    def parse_attributes_gff3(attributes):
        return {key_value.split("=")[0]: key_value.split("=")[1] for key_value in attributes.split(";")}

    def parse_attributes_gtf(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            if attr:
                key, value = attr.strip().split(" ")
                attr_dict[key] = value.replace('"', '')
        return attr_dict

    def parse_file(file_handle, file_type):
        for line in file_handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            feature_type = parts[2]
            attributes = parts[8]

            if file_type == "gff3":
                attr_dict = parse_attributes_gff3(attributes)
            elif file_type == "gtf":
                attr_dict = parse_attributes_gtf(attributes)

            if feature_type == "gene":
                gene_id = attr_dict["gene_id"]
                gene_type = attr_dict["gene_type"]
                try:
                    gene_name = attr_dict["gene_name"]
                except KeyError:
                    gene_name = "."  # Use a placeholder if gene name is not available
                if gene_type in gene_types:
                    process_gene(parts, gene_id, gene_name)

    open_func = gzip.open if annotation_file.endswith(".gz") else open
    file_type = "gff3" if ".gff3" in annotation_file else "gtf"

    with open_func(annotation_file, "rt") as f:
        parse_file(f, file_type)

    return gene_regions, gene_names


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    fetch_start = start_pos - 1  # pysam fetch is 0-based, start-inclusive
    fetch_end = end_pos  # pysam fetch is 0-based, end-exclusive
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, fetch_start, fetch_end):
            PS = read.get_tag("PS") if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag


def get_support_reads(bam_file, vcf_file, chr, start_pos, end_pos):
    """Get reads supporting the variants in the VCF file."""
    fetch_start = start_pos - 1  # pysam fetch is 0-based, start-inclusive
    fetch_end = end_pos  # pysam fetch is 0-based, end-exclusive
    gene_snps = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch(chr, fetch_start, fetch_end):
            gt = record.samples[0]['GT']

            # Skip indels by checking if any alternate allele differs in length from the reference allele
            if any(len(record.ref) != len(alt) for alt in record.alts):
                continue

            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ref_allele = record.ref
                alt_allele = record.alts[0]
                gene_snps[record.pos] = (ref_allele, alt_allele)  # 1-based position

    support_reads = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for pileupcolumn in bam.pileup(chr, fetch_start, fetch_end):
            pos = pileupcolumn.pos + 1
            if pos not in gene_snps:
                continue
            ref_allele, alt_allele = gene_snps[pos]
            allele_counts = {'ref': 0, 'alt': 0, 'other': 0}
            ref_read_names = []
            alt_read_names = []
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_refskip:
                    continue
                if pileupread.is_del:
                    allele_counts['other'] += 1
                    continue
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base == ref_allele:
                    allele_counts['ref'] += 1
                    ref_read_names.append(pileupread.alignment.query_name)
                elif base == alt_allele:
                    allele_counts['alt'] += 1
                    alt_read_names.append(pileupread.alignment.query_name)
                else:
                    allele_counts['other'] += 1
            total_allele_cnt = sum(allele_counts.values())
            if total_allele_cnt > 0 and (allele_counts['ref'] + allele_counts['alt']) / total_allele_cnt >= 0.8:
                support_reads[pos] = (ref_read_names, alt_read_names)  # 1-based position
    return support_reads


# def calculate_ase_pvalue(isoquant_read_assignments, support_reads, gene_id):
#     candidates = support_reads.keys()
#     p_value = 1.0
#     ref_reads_num, alt_reads_num = 0, 0
#     most_significant_pos = None
#     for pos in candidates:
#         ref_reads, alt_reads = support_reads[pos]
#         assigned_reads = set()
#         for isoform_id, reads in isoquant_read_assignments[gene_id].items():
#             assigned_reads.update(reads)
#         ## calculate ASE
#         shared_ref_reads = set(ref_reads).intersection(assigned_reads)
#         shared_alt_reads = set(alt_reads).intersection(assigned_reads)
#         total_ase_counts = len(shared_ref_reads) + len(shared_alt_reads)
#         # print(f"gene_id: {gene_id}, pos: {pos}, ref_reads: {len(ref_reads)}, alt_reads: {len(alt_reads)}, shared_ref_reads: {len(shared_ref_reads)}, shared_alt_reads: {len(shared_alt_reads)}")
#         if total_ase_counts >= 10:
#             p_value_ase = binomtest(len(shared_ref_reads), total_ase_counts, 0.5, alternative='two-sided').pvalue
#             # print(f"gene_id: {gene_id}, pos: {pos}, ref_reads: {len(ref_reads)}, alt_reads: {len(alt_reads)}, shared_ref_reads: {len(shared_ref_reads)}, shared_alt_reads: {len(shared_alt_reads)}, p_value_ase: {p_value_ase}")
#             if p_value_ase < p_value:
#                 p_value = p_value_ase
#                 ref_reads_num = len(shared_ref_reads)
#                 alt_reads_num = len(shared_alt_reads)
#                 most_significant_pos = pos
#     return p_value, ref_reads_num, alt_reads_num, most_significant_pos


def calculate_ase_pvalue(isoquant_read_assignments, support_reads, gene_id):
    candidates = support_reads.keys()
    p_values = []
    pos_list = []
    ref_alt_counts = []
    max_ref_alt_counts = (0, 0)

    # Step 1: Collect p-values
    for pos in candidates:
        ref_reads, alt_reads = support_reads[pos]
        assigned_reads = set()
        for isoform_id, reads in isoquant_read_assignments[gene_id].items():
            assigned_reads.update(reads)

        # Calculate ASE
        shared_ref_reads = set(ref_reads).intersection(assigned_reads)
        shared_alt_reads = set(alt_reads).intersection(assigned_reads)
        total_ase_counts = len(shared_ref_reads) + len(shared_alt_reads)

        if total_ase_counts >= 10:
            p_value_ase = binomtest(len(shared_ref_reads), total_ase_counts, 0.5, alternative='two-sided').pvalue
            p_values.append(p_value_ase)
            pos_list.append(pos)
            ref_alt_counts.append((len(shared_ref_reads), len(shared_alt_reads)))
        else:
            if len(shared_ref_reads) + len(shared_alt_reads) > sum(max_ref_alt_counts):
                max_ref_alt_counts = (len(shared_ref_reads), len(shared_alt_reads))

    # Step 2: Apply Benjamini–Hochberg correction
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

        # Step 3: Find the most significant position with the lowest adjusted p-value
        min_p_value_index = adjusted_p_values.argmin()
        most_significant_pos = pos_list[min_p_value_index]
        ref_reads_num, alt_reads_num = ref_alt_counts[min_p_value_index]
        adjusted_p_value = adjusted_p_values[min_p_value_index]

        return adjusted_p_value, ref_reads_num, alt_reads_num, most_significant_pos

    # # step 2: fisher's method
    # if p_values:
    #     min_p_value_index = np.argmin(p_values)
    #     most_significant_pos = pos_list[min_p_value_index]
    #     ref_reads_num, alt_reads_num = ref_alt_counts[min_p_value_index]
    #     X2 = -2 * sum([np.log(p + 1e-300) for p in p_values])  # avoid log(0)
    #     dof = 2 * len(p_values)
    #     combined_p_value = chi2.sf(X2, dof)
    #     return combined_p_value, ref_reads_num, alt_reads_num, most_significant_pos

    # If no significant p-values, return default values, and the coverage of ref and alt alleles
    ref_reads_num, alt_reads_num = max_ref_alt_counts
    return 1.0, ref_reads_num, alt_reads_num, None


def calculate_asts_pvalue(isoquant_read_assignments, support_reads, gene_id):
    candidates = support_reads.keys()
    p_values = []
    cand_variant_positions = []
    cand_isoforms = []
    cand_contingency_table = []

    isoforms = isoquant_read_assignments[gene_id].keys()
    if len(isoforms) >= 2:
        for pos in candidates:
            ref_reads, alt_reads = support_reads[pos]
            ref_isoform_reads = {isoform_id: 0 for isoform_id in isoforms}
            alt_isoform_reads = {isoform_id: 0 for isoform_id in isoforms}

            for isoform_id, reads in isoquant_read_assignments[gene_id].items():
                shared_ref_reads = set(ref_reads).intersection(reads)
                shared_alt_reads = set(alt_reads).intersection(reads)
                ref_isoform_reads[isoform_id] = len(shared_ref_reads)
                alt_isoform_reads[isoform_id] = len(shared_alt_reads)

            # Ensure consistent ordering using a list of keys
            isoform_ids = list(ref_isoform_reads.keys())
            ref_counts = [ref_isoform_reads[isoform_id] + 1 for isoform_id in isoform_ids]
            alt_counts = [alt_isoform_reads[isoform_id] + 1 for isoform_id in isoform_ids]
            contingency_table = [ref_counts, alt_counts]  # size: [2, num_isoforms]
            _, p_value_asts, _, _ = chi2_contingency(contingency_table)
            p_values.append(p_value_asts)
            cand_variant_positions.append(pos)
            cand_isoforms.append(isoform_ids)
            cand_contingency_table.append(contingency_table)

    # Apply Benjamini–Hochberg correction
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        min_p_value_index = adjusted_p_values.argmin()
        most_significant_pos = cand_variant_positions[min_p_value_index]
        most_significant_isoforms = cand_isoforms[min_p_value_index]
        isoforms_str = ",".join(most_significant_isoforms)
        most_significant_contingency_table = cand_contingency_table[min_p_value_index]
        contingency_table_str = ",".join(
            [f"{most_significant_contingency_table[0][i]}:{most_significant_contingency_table[1][i]}" for i in
             range(len(most_significant_isoforms))])
        adjusted_p_value = adjusted_p_values[min_p_value_index]
        return adjusted_p_value, most_significant_pos, isoforms_str, contingency_table_str

    # # step 2: fisher's method
    # if p_values:
    #     min_p_value_index = np.argmin(p_values)
    #     most_significant_pos = cand_variant_positions[min_p_value_index]
    #     most_significant_isoforms = cand_isoforms[min_p_value_index]
    #     isoforms_str = ",".join(most_significant_isoforms)
    #     most_significant_contingency_table = cand_contingency_table[min_p_value_index]
    #     contingency_table_str = ",".join(
    #         [f"{most_significant_contingency_table[0][i]}:{most_significant_contingency_table[1][i]}" for i in
    #          range(len(most_significant_isoforms))])
    #     X2 = -2 * sum([np.log(p + 1e-300) for p in p_values])  # avoid log(0)
    #     dof = 2 * len(p_values)
    #     combined_p_value = chi2.sf(X2, dof)
    #     return combined_p_value, most_significant_pos, isoforms_str, contingency_table_str

    return 1.0, None, None, None


def calc_ase_asts(isoquant_read_assignments, reads_tag, support_reads, gene_id, p_value_threshold=0.05):
    candidates = support_reads.keys()
    for pos in candidates:
        ref_reads, alt_reads = support_reads[pos]
        assigned_reads = set()
        for isoform_id, reads in isoquant_read_assignments[gene_id].items():
            assigned_reads.update(reads)

        ## calculate ASE
        shared_ref_reads = set(ref_reads).intersection(assigned_reads)
        shared_alt_reads = set(alt_reads).intersection(assigned_reads)
        total_ase_counts = len(shared_ref_reads) + len(shared_alt_reads)
        if total_ase_counts >= 10:
            p_value_ase = binomtest(len(shared_ref_reads), total_ase_counts, 0.5, alternative='two-sided').pvalue
            if p_value_ase <= p_value_threshold:
                print(
                    f"Allele-specific: {gene_id}\t{pos}\t{len(shared_ref_reads)}\t{len(shared_alt_reads)}\t{p_value_ase}")

        ## calculate ASTS
        isoforms = isoquant_read_assignments[gene_id].keys()
        if len(isoforms) >= 2:
            ref_isoform_reads = {isoform_id: 0 for isoform_id in isoforms}
            alt_isoform_reads = {isoform_id: 0 for isoform_id in isoforms}

            for isoform_id, reads in isoquant_read_assignments[gene_id].items():
                shared_ref_reads = set(ref_reads).intersection(reads)
                shared_alt_reads = set(alt_reads).intersection(reads)
                ref_isoform_reads[isoform_id] = len(shared_ref_reads)
                alt_isoform_reads[isoform_id] = len(shared_alt_reads)

            # Ensure consistent ordering using a list of keys
            isoform_ids = list(ref_isoform_reads.keys())
            ref_counts = [ref_isoform_reads[isoform_id] for isoform_id in isoform_ids]
            alt_counts = [alt_isoform_reads[isoform_id] for isoform_id in isoform_ids]
            contingency_table = [ref_counts, alt_counts]

            # Ensure at least two isoforms have more than 10 reads
            if sum(count >= 10 for count in ref_isoform_reads.values()) >= 2:
                _, p_value_asts, _, _ = chi2_contingency(contingency_table)
                if p_value_asts <= p_value_threshold:
                    print(
                        f"Allele-specific transcript structure: {gene_id}\t{pos}\t{p_value_asts}\t{contingency_table}")

        ## calculate haplotype ASE
        containing_reads = set()
        containing_reads.update(ref_reads)
        containing_reads.update(alt_reads)
        shared_reads = containing_reads.intersection(assigned_reads)
        if len(shared_reads) == 0:
            continue
        haplotype_counts = {1: 0, 2: 0}
        for read_name in shared_reads:
            read_info = reads_tag[read_name]
            if read_info['PS'] and read_info['HP']:
                haplotype_counts[read_info['HP']] += 1
        total_counts = sum(haplotype_counts.values())
        if total_counts < 10:
            continue
        p_value = binomtest(haplotype_counts[1], total_counts, 0.5, alternative='two-sided').pvalue
        if p_value <= p_value_threshold:
            print(f"Haplotype-specific: {gene_id}\t{pos}\t{haplotype_counts[1]}\t{haplotype_counts[2]}\t{p_value}")

        ## calculate haplotype ASTS


def analyze_ase_asts(bam_file, vcf_file, gene_id, gene_region, gene_name, isoquant_read_assignments):
    chr = gene_region["chr"]
    start_pos = gene_region["start"]  # 1-based, start-inclusive
    end_pos = gene_region["end"]  # 1-based, end-inclusive
    support_reads = get_support_reads(bam_file, vcf_file, chr, start_pos, end_pos)
    ase_p_value, ref_reads, alt_reads, ase_var_pos = calculate_ase_pvalue(isoquant_read_assignments,
                                                                          support_reads,
                                                                          gene_id)
    asts_p_value, asts_var_pos, isoforms_str, contingency_table_str = calculate_asts_pvalue(isoquant_read_assignments,
                                                                                            support_reads,
                                                                                            gene_id)
    return gene_id, gene_name, chr, ase_var_pos, ase_p_value, ref_reads, alt_reads, \
        asts_var_pos, asts_p_value, isoforms_str, contingency_table_str


# def analyze_asts(bam_file, vcf_file, gene_id, gene_region, gene_name, isoquant_read_assignments):
#     chr = gene_region["chr"]
#     start_pos = gene_region["start"]  # 1-based, start-inclusive
#     end_pos = gene_region["end"]  # 1-based, end-inclusive
#     support_reads = get_support_reads(bam_file, vcf_file, chr, start_pos, end_pos)
#     p_value, variant_pos, isoforms, contingency_table = calculate_asts_pvalue(isoquant_read_assignments, support_reads,
#                                                                               gene_id)
#     return gene_id, gene_name, chr, variant_pos, p_value, isoforms, contingency_table


def multiple_processes_run(bam_file, vcf_file, annotation_file, assignment_file, output_prefix, threads, gene_types,
                           assignment_type, classification_type):
    isoquant_read_assignments = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions, gene_names = get_gene_regions(annotation_file, gene_types)
    gene_args = [(bam_file, vcf_file, gene_id, gene_regions[gene_id], gene_names[gene_id])
                 for gene_id in gene_regions.keys()
                 if gene_id in isoquant_read_assignments]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(isoquant_read_assignments)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(analyze_ase_asts, *gene_data, shared_assignments) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    fase = open(output_prefix + ".ase.tsv", "w")
    fase.write("Gene\tChr\tPos\tRef\tAlt\tP-value\n")
    fasts = open(output_prefix + ".asts.tsv", "w")
    fasts.write("Gene\tChr\tPos\tP-value\tIsoforms\tTable\n")
    # # Apply Benjamini–Hochberg correction
    # ase_pvalues = np.array([result[4] for result in results])
    # asts_p_values = np.array([result[8] for result in results])
    # reject_ase, adjusted_ase_p_values, _, _ = multipletests(ase_pvalues, alpha=0.05, method='fdr_bh')
    # reject_asts, adjusted_asts_p_values, _, _ = multipletests(asts_p_values, alpha=0.05, method='fdr_bh')
    for idx, (gene_id, gene_name, chr, ase_var_pos, ase_p_value, ref_reads, alt_reads, asts_var_pos, asts_p_value,
              isoforms_str, contingency_table_str) in enumerate(results):
        # ase_p_value = adjusted_ase_p_values[idx]
        # asts_p_value = adjusted_asts_p_values[idx]
        fase.write(f"{gene_name}\t{chr}\t{ase_var_pos}\t{ref_reads}\t{alt_reads}\t{ase_p_value}\n")
        fasts.write(f"{gene_name}\t{chr}\t{asts_var_pos}\t{asts_p_value}\t{isoforms_str}\t{contingency_table_str}\n")
    fase.close()
    fasts.close()


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="Analyze allele-specific expression and transcript structure.")
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("-v", "--vcf_file", help="VCF file", required=True)
    parse.add_argument("-a", "--annotation_file", help="Annotation file", required=True)
    parse.add_argument("-i", "--assignment_file", help="Assignment file", required=True)
    parse.add_argument("-o", "--output_prefix", help="Prefix of output file", required=True)
    parse.add_argument("-t", "--threads", help="Number of threads", type=int, default=1)
    parse.add_argument("--gene_types", help="Gene types", nargs="+", default=["protein_coding", "lncRNA"])
    parse.add_argument("--assignment_type", help="Assignment type", nargs="+",
                       default={"unique", "unique_minor_difference"})
    parse.add_argument("--classification_type", help="Classification type", nargs="+",
                       default={"full_splice_match", "incomplete_splice_match", "mono_exon_match"})
    args = parse.parse_args()
    multiple_processes_run(args.bam_file, args.vcf_file, args.annotation_file, args.assignment_file, args.output_prefix,
                           args.threads, args.gene_types, args.assignment_type, args.classification_type)
