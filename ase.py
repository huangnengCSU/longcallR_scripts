import math
import pysam
from collections import defaultdict
import gzip
from intervaltree import Interval, IntervalTree
from multiprocessing import Manager
import concurrent.futures
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import argparse
import time


def get_gene_regions(annotation_file, gene_types):
    """Parse gene, exon, and intron regions from a GFF3 or GTF file.
    :param annotation_file: Path to the annotation file
    :return: Gene regions, exon regions, and intron regions
    """
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unknown annotation file format"

    gene_regions = {}
    gene_names = {}
    gene_strands = {}
    exon_regions = defaultdict(lambda: defaultdict(list))
    intron_regions = defaultdict(lambda: defaultdict(list))

    def process_gene(parts, gene_id, gene_name):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}  # 1-based, start-inclusive, end-inclusive
        gene_names[gene_id] = gene_name
        strand = parts[6]
        gene_strands[gene_id] = strand

    def process_exon(parts, gene_id, transcript_id):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        exon_regions[gene_id][transcript_id].append((chr, start, end))  # 1-based, start-inclusive, end-inclusive

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
            elif feature_type == "exon":
                gene_type = attr_dict["gene_type"]
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
                if gene_type in gene_types:
                    process_exon(parts, gene_id, transcript_id)

    open_func = gzip.open if annotation_file.endswith(".gz") else open
    file_type = "gff3" if ".gff3" in annotation_file else "gtf"

    with open_func(annotation_file, "rt") as f:
        parse_file(f, file_type)

    # Calculate intron regions based on exons
    for gene_id, transcripts in exon_regions.items():
        for transcript_id, exons in transcripts.items():
            if len(exons) == 1:
                continue
            exons_sorted = sorted(exons, key=lambda x: x[1])
            for i in range(1, len(exons_sorted)):
                intron_start = exons_sorted[i - 1][2] + 1
                intron_end = exons_sorted[i][1] - 1
                if intron_start < intron_end:
                    intron_regions[gene_id][transcript_id].append(
                        (exons_sorted[i - 1][0], intron_start, intron_end))  # 1-based, start-inclusive, end-inclusive

    return gene_regions, gene_names, gene_strands, exon_regions, intron_regions


def merge_gene_exon_regions(exon_regions):
    """Merge transcript exons into gene regions."""
    # merged_genes_exons_sorted_by_start = dict()  # key: chr, value: list of sorted (collapsed_exons, gene_id, gene_name)

    # merged_genes_exons, key: chr, value: dict of gene_id: [(start, end), ..., (start, end)]
    merged_genes_exons = defaultdict(lambda: defaultdict(list))

    for gene_id, transcripts in exon_regions.items():
        collapsed_exons = IntervalTree()
        chromosome = None

        # Iterate over transcripts and collect intervals
        for transcript_id, exons in transcripts.items():
            for (chr, start, end) in exons:
                if chromosome is None:
                    chromosome = chr
                else:
                    assert chromosome == chr, f"Error: Inconsistent chromosome in gene {gene_id}"
                collapsed_exons.add(Interval(start, end + 1))  # Interval is left-inclusive, right-exclusive

        # Merge overlapping intervals and adjust to 1-based closed intervals
        collapsed_exons.merge_overlaps()
        collapsed_exons = sorted((interval.begin, interval.end - 1) for interval in collapsed_exons)
        merged_genes_exons[chromosome][gene_id].extend(collapsed_exons)

        # # Insert collapsed exons into merged_genes_exons_sorted_by_start using binary search
        # if chromosome not in merged_genes_exons_sorted_by_start:
        #     merged_genes_exons_sorted_by_start[chromosome] = [(collapsed_exons, gene_id, gene_names[gene_id])]
        # else:
        #     # Extract start positions for binary search
        #     start_positions = [elem[0][0][0] for elem in merged_genes_exons_sorted_by_start[chromosome]]
        #     insert_position = bisect.bisect_left(start_positions, collapsed_exons[0][0])
        #     merged_genes_exons_sorted_by_start[chromosome].insert(insert_position,
        #                                                           (collapsed_exons, gene_id, gene_names[gene_id]))
    # return merged_genes_exons_sorted_by_start
    return merged_genes_exons


def assign_reads_to_gene(bam_file, merged_genes_exons):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    trees = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals = defaultdict(lambda: defaultdict(IntervalTree))  # key: chr, value: dict of gene_id: IntervalTree
    for chrom in merged_genes_exons.keys():
        for gene_id, merged_exons in merged_genes_exons[chrom].items():
            gene_region = (merged_exons[0][0], merged_exons[-1][1])  # 1-based, start-inclusive, end-inclusive
            trees[chrom].add(Interval(gene_region[0], gene_region[1] + 1, gene_id))

            # Build IntervalTree for exon regions within the gene
            for exon_start, exon_end in merged_exons:
                gene_intervals[chrom][gene_id].add(Interval(exon_start, exon_end + 1))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            chromosome = read.reference_name
            start_pos = read.reference_start  # 0-based, inclusive
            end_pos = read.reference_end  # 0-based, exclusive
            if chromosome not in trees:
                continue
            # query should be 1-based, left-inclusive, right-exclusive
            overlapping_intervals = trees[chromosome].overlap(start_pos + 1, end_pos + 1)
            candidate_gene_ids = [interval.data for interval in overlapping_intervals]

            # parse cigar string to get the splice alignment regions
            cigar = read.cigartuples
            splice_regions = []  # list of splice alignment positions of a read, 1-based and start/end inclusive
            current_pos = start_pos + 1  # 1-based
            shift = 0
            for operation, length in cigar:
                if operation in {0, 2, 7, 8}:
                    shift += length
                elif operation == 3:
                    if shift > 0:
                        splice_regions.append((current_pos, current_pos + shift - 1))
                    current_pos += (shift + length)  # Move past the skipped region
                    shift = 0
            if shift > 0:
                splice_regions.append((current_pos, current_pos + shift - 1))

            read_overlap_length = {}  # key: gene_id, value: overlap_length
            # calculate the overlap of splice regions with gene exons to determine assignment of read to gene
            for gene_id in candidate_gene_ids:
                if gene_id not in gene_intervals[chromosome]:
                    continue
                overlap_length = 0
                for splice_region in splice_regions:
                    overlap_length += sum(
                        max(0, min(splice_region[1], interval.end - 1) - max(splice_region[0], interval.begin) + 1)
                        for interval in gene_intervals[chromosome][gene_id].overlap(*splice_region)
                    )
                read_overlap_length[gene_id] = overlap_length
            if read_overlap_length:
                best_gene_id = max(read_overlap_length, key=read_overlap_length.get)
                read_assignment[read.query_name] = best_gene_id
    return read_assignment


def process_chunk(bam_file, chromosome, start, end, shared_tree, shared_gene_intervals):
    read_assignment = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_unmapped:
                continue
            start_pos = read.reference_start
            end_pos = read.reference_end

            # query should be 1-based, left-inclusive, right-exclusive
            overlapping_intervals = shared_tree.overlap(start_pos + 1, end_pos + 1)
            candidate_gene_ids = [interval.data for interval in overlapping_intervals]

            # parse cigar string to get the splice alignment regions
            cigar = read.cigartuples
            splice_regions = []  # list of splice alignment positions of a read, 1-based and start/end inclusive
            current_pos = start_pos + 1  # 1-based
            shift = 0
            for operation, length in cigar:
                if operation in {0, 2, 7, 8}:
                    shift += length
                elif operation == 3:
                    if shift > 0:
                        splice_regions.append((current_pos, current_pos + shift - 1))
                    current_pos += (shift + length)  # Move past the skipped region
                    shift = 0
            if shift > 0:
                splice_regions.append((current_pos, current_pos + shift - 1))

            read_overlap_length = {}  # key: gene_id, value: overlap_length
            # calculate the overlap of splice regions with gene exons to determine assignment of read to gene
            for gene_id in candidate_gene_ids:
                if gene_id not in shared_gene_intervals:
                    continue
                overlap_length = 0
                for splice_region in splice_regions:
                    overlap_length += sum(
                        max(0, min(splice_region[1], interval.end - 1) - max(splice_region[0], interval.begin) + 1)
                        for interval in shared_gene_intervals[gene_id].overlap(*splice_region)
                    )
                read_overlap_length[gene_id] = overlap_length
            if read_overlap_length:
                best_gene_id = max(read_overlap_length, key=read_overlap_length.get)
                read_assignment[read.query_name] = best_gene_id
    return read_assignment


def assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    trees_by_chr = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals_by_chr = defaultdict(
        lambda: defaultdict(IntervalTree))  # key: chr, value: dict of gene_id: IntervalTree
    for chrom in merged_genes_exons.keys():
        for gene_id, merged_exons in merged_genes_exons[chrom].items():
            gene_region = (merged_exons[0][0], merged_exons[-1][1])  # 1-based, start-inclusive, end-inclusive
            trees_by_chr[chrom].add(Interval(gene_region[0], gene_region[1] + 1, gene_id))

            # Build IntervalTree for exon regions within the gene
            for exon_start, exon_end in merged_exons:
                gene_intervals_by_chr[chrom][gene_id].add(Interval(exon_start, exon_end + 1))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = bam.references
        chromosome_lengths = dict(zip(bam.references, bam.lengths))

    chunks = []
    for chromosome in chromosomes:
        total_length = chromosome_lengths[chromosome]
        chunk_size = max(1, math.ceil(total_length / threads))
        for i in range(threads):
            start = i * chunk_size  # 0-based, inclusive
            end = min((i + 1) * chunk_size, total_length)  # 0-based, exclusive
            tree = trees_by_chr[chromosome]
            gene_intervals = gene_intervals_by_chr[chromosome]
            chunks.append((bam_file, chromosome, start, end, tree, gene_intervals))

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_chunk, *chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            read_assignment.update(result)

    return read_assignment


def transform_read_assignment(read_assignment):
    """select the best gene assignment for each read"""
    gene_assigned_reads = defaultdict(list)  # key: gene_id, value: list of read_name
    for read_name, gene_id in read_assignment.items():
        gene_assigned_reads[gene_id].append(read_name)
    return gene_assigned_reads


def load_whole_genome_phased_vcf(vcf_file):
    """
    Load phased variants from a whole-genome VCF file.
    :param vcf_file:
    :return: wg_vcfs: key is chr:pos, value is directory{genotype, paternal, maternal}
    """
    wg_vcfs = {}  # key is chr:pos, value is 0|1 or 1|0
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]['GT']
            # Skip indels by checking if any alternate allele differs in length from the reference allele
            if any(len(record.ref) != len(alt) for alt in record.alts):
                continue
            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ref_allele = record.ref
                alt_allele = record.alts[0]
                if gt == (0, 1):
                    wg_vcfs[f"{record.contig}:{record.pos}"] = {"gt": gt,
                                                                "pat": alt_allele,
                                                                "mat": ref_allele}
                else:
                    wg_vcfs[f"{record.contig}:{record.pos}"] = {"gt": gt,
                                                                "pat": ref_allele,
                                                                "mat": alt_allele}
    return wg_vcfs


def load_longcallR_phased_vcf(vcf_file):
    """
    get the longcallR phased vcf
    :param vcf_file:
    :return: rna_vcfs: key is ps, value is variants
    """
    import pysam
    from collections import defaultdict
    rna_vcfs = defaultdict(list)  # key is ps, value is positions
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]['GT']
            # Skip indels by checking if any alternate allele differs in length from the reference allele
            if any(len(record.ref) != len(alt) for alt in record.alts):
                continue
            # Check for only phased heterozygous variants (0|1 or 1|0)
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ps = record.samples[0].get('PS', None)
                if ps:
                    rna_vcfs[ps].append(f"{record.contig}:{record.pos}")  # 1-based
    return rna_vcfs


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, start_pos, end_pos):
            PS = read.get_tag("PS") if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag


def calculate_ase_pvalue(bam_file, gene_id, gene_name, gene_region, min_count, gene_assigned_reads):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    # assigned_reads = set()
    #
    # # Collect assigned reads from isoquant_read_assignments
    # for isoform_id, reads in isoquant_read_assignments[gene_id].items():
    #     assigned_reads.update(reads)
    assigned_reads = set(gene_assigned_reads[gene_id])

    # Track haplotype counts for each phase set
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        if rname in reads_tag:
            ps = reads_tag[rname]["PS"]
            hp = reads_tag[rname]["HP"]
            if ps and hp:
                phase_set_hap_count[ps][hp] += 1

    # Step 1: Collect p-values for each phase set
    p_values = []
    phase_sets = []
    hap_counts = []

    for ps, hap_count in phase_set_hap_count.items():
        if hap_count[1] + hap_count[2] < min_count:
            continue
        # Binomial test p-value
        total_reads = hap_count[1] + hap_count[2]
        p_value_ase = binomtest(hap_count[1], total_reads, 0.5, alternative='two-sided').pvalue
        p_values.append(p_value_ase)
        phase_sets.append(ps)
        hap_counts.append((hap_count[1], hap_count[2]))

    # Step 2: Apply Benjamini–Hochberg correction
    if p_values:
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

        # Step 3: Find the most significant adjusted p-value
        min_p_value_index = adjusted_p_values.argmin()
        most_significant_phase_set = phase_sets[min_p_value_index]
        h1_count, h2_count = hap_counts[min_p_value_index]
        adjusted_p_value = adjusted_p_values[min_p_value_index]

        return (gene_name, gene_region["chr"], adjusted_p_value, most_significant_phase_set, h1_count, h2_count)

    # # fisher's method
    # if p_values:
    #     min_p_value_index = np.argmin(p_values)
    #     most_significant_phase_set = phase_sets[min_p_value_index]
    #     h1_count, h2_count = hap_counts[min_p_value_index]
    #     X2 = -2 * sum([np.log(p + 1e-300) for p in p_values])  # avoid log(0)
    #     dof = 2 * len(p_values)
    #     combined_p_value = chi2.sf(X2, dof)
    #     return (gene_name, gene_region["chr"], combined_p_value, most_significant_phase_set, h1_count, h2_count)

    # If no valid p-values are found, return defaults
    return (gene_name, gene_region["chr"], 1.0, ".", 0, 0)


def calculate_ase_pvalue_pat_mat(bam_file, gene_id, gene_name, gene_region, min_count, gene_assigned_reads,
                                 rna_vcfs, wg_vcfs):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    # assigned_reads = set()
    #
    # # Collect assigned reads from isoquant_read_assignments
    # for isoform_id, reads in isoquant_read_assignments[gene_id].items():
    #     assigned_reads.update(reads)
    assigned_reads = set(gene_assigned_reads[gene_id])

    # Track haplotype counts for each phase set
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
        if rname in reads_tag:
            ps = reads_tag[rname]["PS"]
            hp = reads_tag[rname]["HP"]
            if ps and hp:
                phase_set_hap_count[ps][hp] += 1

    # Step 1: Collect p-values for each phase set
    p_values = []
    phase_sets = []
    hap_counts = []

    for ps, hap_count in phase_set_hap_count.items():
        if hap_count[1] + hap_count[2] < min_count:
            continue
        # Binomial test p-value
        total_reads = hap_count[1] + hap_count[2]
        p_value_ase = binomtest(hap_count[1], total_reads, 0.5, alternative='two-sided').pvalue
        p_values.append(p_value_ase)
        phase_sets.append(ps)
        hap_counts.append((hap_count[1], hap_count[2]))

    if p_values:
        # # Step 2: fisher's method
        # min_p_value_index = np.argmin(p_values)
        # most_significant_phase_set = phase_sets[min_p_value_index]
        # h1_count, h2_count = hap_counts[min_p_value_index]
        # X2 = -2 * sum([np.log(p + 1e-300) for p in p_values])  # avoid log(0)
        # dof = 2 * len(p_values)
        # combined_p_value = chi2.sf(X2, dof)

        # apply Benjamini–Hochberg correction
        reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
        min_p_value_index = adjusted_p_values.argmin()
        most_significant_phase_set = phase_sets[min_p_value_index]
        h1_count, h2_count = hap_counts[min_p_value_index]
        adjusted_p_value = adjusted_p_values[min_p_value_index]

        # determine paternal and maternal alleles for H1 and H2
        ps_variants = rna_vcfs.get(most_significant_phase_set, [])
        ps_reads = [rname for rname in assigned_reads if
                    rname in reads_tag and reads_tag[rname]["PS"] == most_significant_phase_set]
        h1_reads = [rname for rname in ps_reads if reads_tag[rname]["HP"] == 1]
        h2_reads = [rname for rname in ps_reads if reads_tag[rname]["HP"] == 2]

        ps_variant_pos = [int(pos.split(":")[1]) - 1 for pos in ps_variants]  # 0-based
        reads_pat_mat_cnt = defaultdict(
            lambda: {"pat": 0, "mat": 0})  # key: read name, value: {paternal: count, maternal: count}
        for pileupcolumn in pysam.AlignmentFile(bam_file, "rb").pileup(gene_region["chr"], gene_region["start"] - 1,
                                                                       gene_region["end"], max_depth=100000):
            if pileupcolumn.pos not in ps_variant_pos or f"{gene_region['chr']}:{pileupcolumn.pos + 1}" not in wg_vcfs:
                continue
            for pileup_read in pileupcolumn.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    read_name = pileup_read.alignment.query_name
                    if read_name in ps_reads:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                        # wg_vcfs is 1-based, pileupcolumn.pos is 0-based
                        if base in wg_vcfs[f"{gene_region['chr']}:{pileupcolumn.pos + 1}"]["pat"]:
                            reads_pat_mat_cnt[read_name]["pat"] += 1
                        elif base in wg_vcfs[f"{gene_region['chr']}:{pileupcolumn.pos + 1}"]["mat"]:
                            reads_pat_mat_cnt[read_name]["mat"] += 1
                        else:
                            continue
        h1_pat_cnt, h1_mat_cnt = 0, 0
        for reads in h1_reads:
            if reads in reads_pat_mat_cnt:
                if reads_pat_mat_cnt[reads]["pat"] > reads_pat_mat_cnt[reads]["mat"]:
                    h1_pat_cnt += 1
                elif reads_pat_mat_cnt[reads]["pat"] < reads_pat_mat_cnt[reads]["mat"]:
                    h1_mat_cnt += 1
                else:
                    continue
        h2_pat_cnt, h2_mat_cnt = 0, 0
        for reads in h2_reads:
            if reads in reads_pat_mat_cnt:
                if reads_pat_mat_cnt[reads]["pat"] > reads_pat_mat_cnt[reads]["mat"]:
                    h2_pat_cnt += 1
                elif reads_pat_mat_cnt[reads]["pat"] < reads_pat_mat_cnt[reads]["mat"]:
                    h2_mat_cnt += 1
                else:
                    continue
        return (gene_name, gene_region["chr"], adjusted_p_value, most_significant_phase_set,
                h1_count, h2_count, h1_pat_cnt, h1_mat_cnt, h2_pat_cnt, h2_mat_cnt)

    # If no valid p-values are found, return defaults
    return (gene_name, gene_region["chr"], 1.0, ".", 0, 0, 0, 0, 0, 0)


def analyze_ase_genes(annotation_file, bam_file, out_file, threads, gene_types, min_support):
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    read_assignment = assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads)
    gene_assigned_reads = transform_read_assignment(read_assignment)

    gene_args = [(bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support)
                 for gene_id in gene_regions.keys()
                 if gene_id in gene_assigned_reads]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(gene_assigned_reads)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(calculate_ase_pvalue, *gene_data, shared_assignments) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    # # apply Benjamini–Hochberg correction
    # p_values = [result[2] for result in results]
    # reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    with open(out_file, "w") as f:
        f.write("Gene\tChr\tPS\tH1\tH2\tP-value\n")
        for idx, (gene_name, chrom, p_value, ps, h1, h2) in enumerate(results):
            # p_value = adjusted_p_values[idx]
            f.write(f"{gene_name}\t{chrom}\t{ps}\t{h1}\t{h2}\t{p_value}\n")


def analyze_ase_genes_pat_mat(annotation_file, bam_file, vcf_file1, vcf_file2, out_file, threads, gene_types,
                              min_support):
    rna_vcfs = load_longcallR_phased_vcf(vcf_file1)
    wg_vcfs = load_whole_genome_phased_vcf(vcf_file2)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    merged_genes_exons = merge_gene_exon_regions(exon_regions)
    read_assignment = assign_reads_to_gene_parallel(bam_file, merged_genes_exons, threads)
    gene_assigned_reads = transform_read_assignment(read_assignment)

    gene_args = [(bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support)
                 for gene_id in gene_regions.keys()
                 if gene_id in gene_assigned_reads]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(gene_assigned_reads)
        shared_rna_vcfs = manager.dict(rna_vcfs)
        shared_wg_vcfs = manager.dict(wg_vcfs)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(calculate_ase_pvalue_pat_mat, *gene_data, shared_assignments, shared_rna_vcfs,
                                       shared_wg_vcfs) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    # # apply Benjamini–Hochberg correction
    # p_values = [result[2] for result in results]
    # reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    with open(out_file, "w") as f:
        f.write("Gene\tChr\tPS\tH1\tH2\tP-value\tH1_Paternal\tH1_Maternal\tH2_Paternal\tH2_Maternal\n")
        for idx, (gene_name, chrom, p_value, ps, h1, h2, h1_pat, h1_mat, h2_pat, h2_mat) in enumerate(results):
            # p_value = adjusted_p_values[idx]
            f.write(f"{gene_name}\t{chrom}\t{ps}\t{h1}\t{h2}\t{p_value}\t{h1_pat}\t{h1_mat}\t{h2_pat}\t{h2_mat}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="phased BAM file")
    parser.add_argument("--vcf1", required=False, help="longcallR phased vcf file", default=None)
    parser.add_argument("--vcf2", required=False, help="whole genome haplotype phased vcf file", default=None)
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-o", "--output", required=True, help="prefix of output file")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Number of process to run")
    parser.add_argument("--gene_types", type=str, nargs="+", default=["protein_coding", "lncRNA"],
                        help='Gene types to be analyzed. Default is ["protein_coding", "lncRNA"]', )
    parser.add_argument("--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")

    args = parser.parse_args()

    gene_types = set(args.gene_types)

    if args.vcf1 and args.vcf2:
        analyze_ase_genes_pat_mat(args.annotation, args.bam, args.vcf1, args.vcf2, args.output + ".patmat.ase.tsv",
                                  args.processes, gene_types, args.min_support)
    else:
        analyze_ase_genes(args.annotation, args.bam, args.output + ".ase.tsv", args.processes, gene_types,
                          args.min_support)
