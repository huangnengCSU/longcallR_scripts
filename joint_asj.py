import argparse
import concurrent.futures
from multiprocessing import Manager
import gzip
import math
from collections import defaultdict
import networkx as nx
import csv
import numpy as np
import pysam
from scipy.stats import fisher_exact, power_divergence, chi2, chi2_contingency
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
from intervaltree import IntervalTree, Interval
import time


def get_tissue_readnames(tissue_readnames):
    tissue_readname_map = {}
    with open(tissue_readnames, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            read_name, tissue_name = row
            tissue_readname_map[read_name] = tissue_name
    return tissue_readname_map


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
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            key, value = attr.strip().split("=")
            attr_dict[key] = value.replace('"', '')
        return attr_dict

    def parse_attributes_gtf(attributes):
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            if attr:
                key, value = attr.strip().split(" ")
                if key == "tag":
                    attr_dict[key] = attr_dict.get(key, []) + [value.replace('"', '')]
                else:
                    attr_dict[key] = value.replace('"', '')
        attr_dict["tag"] = ",".join(attr_dict.get("tag", []))
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
                tag = attr_dict.get("tag", "")
                try:
                    gene_name = attr_dict["gene_name"]
                except KeyError:
                    gene_name = "."  # Use a placeholder if gene name is not available
                if gene_type in gene_types and "readthrough" not in tag:
                    process_gene(parts, gene_id, gene_name)
            elif feature_type == "exon":
                gene_type = attr_dict["gene_type"]
                transcript_id = attr_dict["transcript_id"]
                gene_id = attr_dict["gene_id"]
                tag = attr_dict.get("tag", "")
                if gene_type in gene_types and "readthrough" not in tag:
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


def get_exon_intron_regions(read, ref_seq, no_gtag):
    exon_regions = []  # 1-based, start-inclusive, end-inclusive
    intron_regions = []  # 1-based, start-inclusive, end-inclusive, gt-ag tag
    reference_start = read.reference_start + 1  # 1-based
    current_position = reference_start
    for cigartuple in read.cigartuples:
        operation, length = cigartuple
        if operation in {0, 7, 8}:
            if exon_regions and exon_regions[-1][1] + 1 == current_position:
                # Extend the last exon region if it is contiguous
                exon_regions[-1] = (exon_regions[-1][0], exon_regions[-1][1] + length)
            else:
                # Start a new exon region
                exon_start = current_position
                exon_end = current_position + length - 1
                exon_regions.append((exon_start, exon_end))
            current_position += length
        elif operation == 2:  # 'D' operation represents deletions (still part of exon on the reference)
            if exon_regions and exon_regions[-1][1] + 1 == current_position:
                # Extend the last exon region if it is contiguous
                exon_regions[-1] = (exon_regions[-1][0], exon_regions[-1][1] + length)
            else:
                # Start a new exon region, case: 100M20N10D100M
                exon_start = current_position
                exon_end = current_position + length - 1
                exon_regions.append((exon_start, exon_end))
            current_position += length
        elif operation == 3:  # 'N' operation represents skipped region (intron)
            intron_start = current_position  # 1-based, start-inclusive
            intron_end = current_position + length - 1  # 1-based, end-inclusive
            intron_left_seq = ref_seq[intron_start - 1: intron_start + 1].upper()
            intron_right_seq = ref_seq[intron_end - 2: intron_end].upper()
            if no_gtag:
                intron_regions.append((intron_start, intron_end, False))
            else:
                if (intron_left_seq == "GT" and intron_right_seq == "AG") or (
                        intron_left_seq == "CT" and intron_right_seq == "AC"):
                    intron_regions.append((intron_start, intron_end, True))
                else:
                    intron_regions.append((intron_start, intron_end, False))
            current_position += length
        else:
            pass
    return exon_regions, intron_regions


def merge_gene_exon_regions(exon_regions):
    """Merge transcript exons into gene regions."""
    # merged_genes_exons_sorted_by_start = dict()  # key: chr, value: list of sorted (collapsed_exons, gene_id, gene_name)

    # merged_genes_exons, key: chr, value: dict of gene_id: [(start, end), ..., (start, end)]
    merged_genes_exons = defaultdict(lambda: defaultdict(list))
    for gene_id, transcripts in exon_regions.items():
        collapsed_exons = IntervalTree()
        chromosome = None
        chr_set = set()
        for transcript_id, exons in transcripts.items():
            for (chr, start, end) in exons:
                chr_set.add(chr)
        if len(chr_set) > 1:
            # this gene has exons on multiple chromosomes
            continue
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
    return merged_genes_exons


def process_chunk(bam_file, chromosome, start, end, ref_seq, no_gtag, shared_tree, shared_gene_intervals):
    read_assignment = {}
    reads_positions = {}
    reads_tags = {}
    reads_exons = {}
    reads_junctions = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_unmapped:
                continue
            start_pos = read.reference_start
            end_pos = read.reference_end

            # get read positions and read tags
            if not read.has_tag("HP"):
                continue
            HP_tag = read.get_tag("HP")
            if read.has_tag("PS"):
                PS_tag = read.get_tag("PS")
            else:
                PS_tag = "."
            reads_tags[read.query_name] = {"PS": PS_tag, "HP": HP_tag}
            # get read start position and end position, 1-based, start-inclusive, end-inclusive
            reads_positions[read.query_name] = (read.reference_start + 1, read.reference_end)

            # get all exons and introns
            exon_regions, intron_regions = get_exon_intron_regions(read, ref_seq, no_gtag)
            reads_exons[read.query_name] = exon_regions
            reads_junctions[read.query_name] = intron_regions

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
    return read_assignment, reads_positions, reads_tags, reads_exons, reads_junctions


def load_reads(bam_file, genome_dict, merged_genes_exons, threads, no_gtag):
    """Assign reads to genes based on their alignment positions."""

    # read_assignment, key: read_name, value: gene_id
    read_assignment = {}

    reads_positions = {}  # key: read_name, value: (start, end)
    reads_tags = {}  # key: read_name, value: {"PS": phase set, "HP": haplotype}

    # key: read_name, value: exons(list)/introns(list)
    reads_exons = {}
    reads_junctions = {}

    trees_by_chr = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals_by_chr = defaultdict(lambda: defaultdict(IntervalTree))  # key1: chr, key2: gene_id
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
        if chromosome not in genome_dict:
            continue
        total_length = chromosome_lengths[chromosome]  # 1-based
        chunk_size = max(1, math.ceil(total_length / threads))
        for i in range(threads):
            start = i * chunk_size  # 0-based, inclusive
            end = min((i + 1) * chunk_size, total_length)  # 0-based, exclusive
            if start >= end:
                continue
            tree = trees_by_chr[chromosome]
            gene_intervals = gene_intervals_by_chr[chromosome]
            chunks.append((bam_file, chromosome, start, end, genome_dict[chromosome], no_gtag, tree, gene_intervals))

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_chunk, *chunk) for chunk in chunks]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            read_assignment.update(result[0])
            reads_positions.update(result[1])
            reads_tags.update(result[2])
            reads_exons.update(result[3])
            reads_junctions.update(result[4])

    return read_assignment, reads_positions, reads_tags, reads_exons, reads_junctions


def transform_read_assignment(read_assignment):
    """select the best gene assignment for each read"""
    gene_assigned_reads = defaultdict(list)  # key: gene_id, value: list of read_name
    for read_name, gene_id in read_assignment.items():
        gene_assigned_reads[gene_id].append(read_name)
    return gene_assigned_reads


def parse_reads_from_alignment(bam_file, chr, start_pos, end_pos):
    """Parse reads from a BAM file that overlap a specific region.
    :param bam_file: Path to the BAM file
    :param chr: Chromosome
    :param start_pos: Start position, 1-based inclusive
    :param end_pos: End position, 1-based inclusive
    :return:
    """
    reads_exons = {}
    reads_junctions = {}
    reads_positions = {}  # 1-based, start-inclusive, end-inclusive
    reads_tags = {}  # key: read name, value: {"PS": phase set, "HP": haplotype}
    samfile = pysam.AlignmentFile(bam_file, "rb")
    fetch_start = start_pos - 1  # 0-based, start-inclusive
    fetch_end = end_pos  # 0-based, end-exclusive
    contigs = samfile.references
    if chr not in contigs:
        return reads_positions, reads_exons, reads_junctions, reads_tags
    for read in samfile.fetch(chr, fetch_start, fetch_end):
        # read is not phased, ignore
        # if not read.has_tag("PS") or not read.has_tag("HP"):
        #     continue
        if not read.has_tag("HP"):
            continue
        if read.has_tag("PS"):
            PS_tag = read.get_tag("PS")
        else:
            PS_tag = "."
        HP_tag = read.get_tag("HP")
        reads_tags[read.query_name] = {"PS": PS_tag, "HP": HP_tag}
        # get read start position and end position, 1-based, start-inclusive, end-inclusive
        read_start = read.reference_start + 1
        read_end = read.reference_end
        reads_positions[read.query_name] = (read_start, read_end)
        # get all exons and introns
        exon_regions, intron_regions = get_exon_intron_regions(read)
        reads_exons[read.query_name] = exon_regions
        reads_junctions[read.query_name] = intron_regions
    return reads_positions, reads_exons, reads_junctions, reads_tags


def cluster_exons(reads_exons, min_count=10):
    exons = {}  # key: (start, end), value: count
    for read_name, exon_regions in reads_exons.items():
        # single exon
        if len(exon_regions) == 1:
            pass
        # two exons
        if len(exon_regions) == 2:
            pass
        # more than two exons
        if len(exon_regions) > 2:
            for i, exon_region in enumerate(exon_regions):
                if i == 0 or i == len(exon_regions) - 1:
                    continue
                exons[exon_region] = exons.get(exon_region, 0) + 1
    # delete exons with count less than 10
    exons = {k: v for k, v in exons.items() if v >= min_count}
    return exons


def cluster_junctions(reads_junctions, min_count=10):
    junctions = {}  # key: (start, end), value: count
    for read_name, junction_regions in reads_junctions.items():
        for junction_region in junction_regions:
            junctions[junction_region] = junctions.get(junction_region, 0) + 1
    # delete junctions with count less than 10
    junctions = {k: v for k, v in junctions.items() if v >= min_count}
    return junctions


def cluster_junctions_connected_components(reads_junctions, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count

    # Count occurrences of each junction region
    for read_name, junction_regions in reads_junctions.items():
        for junction_region in junction_regions:
            junctions[junction_region] = junctions.get(junction_region, 0) + 1

    # Remove junctions with count less than min_count
    junctions = {k: v for k, v in junctions.items() if v >= min_count}

    # Create a graph
    G = nx.Graph()

    # Add nodes for each junction
    for junction in junctions.keys():
        G.add_node(junction)  # "junction = (intron_start, intron_end)"

    # Add edges between overlapping junctions (sharing start or end positions)
    junction_list = list(junctions.keys())
    for i in range(len(junction_list)):
        for j in range(i + 1, len(junction_list)):
            start1, end1 = junction_list[i]
            start2, end2 = junction_list[j]

            # Check if they share a start or end position
            if start1 == start2 or end1 == end2:
                G.add_edge(junction_list[i], junction_list[j])

    # Find connected components
    connected_components = list(nx.connected_components(G))

    # Collect the clusters of junctions
    for component in connected_components:
        junctions_clusters.append(list(component))

    return junctions_clusters, junctions


def cluster_junctions_exons_connected_components(reads_junctions, reads_exons, min_count=10):
    junctions_clusters = []
    junctions = {}  # key: (start, end), value: count
    gt_ag_dict = {}  # key: (start, end), value: gt-ag tag
    for read_name, list_of_junctions in reads_junctions.items():
        for (start, end, tag) in list_of_junctions:
            junction_region = (start, end)
            junctions[junction_region] = junctions.get(junction_region, 0) + 1
            gt_ag_dict[junction_region] = tag
    junctions = {k: v for k, v in junctions.items() if v >= min_count}
    exons = {}  # key: (start, end), value: count
    for read_name, exon_regions in reads_exons.items():
        if len(exon_regions) == 0:
            pass
        # single exon
        if len(exon_regions) == 1:
            pass
        # two exons
        if len(exon_regions) == 2:
            pass
        # more than two exons
        if len(exon_regions) > 2:
            for i, exon_region in enumerate(exon_regions):
                if i == 0 or i == len(exon_regions) - 1:
                    continue
                exons[exon_region] = exons.get(exon_region, 0) + 1
    exons = {k: v for k, v in exons.items() if v >= min_count}
    # Create a graph
    G = nx.Graph()
    # Add nodes for each junction
    for junction in junctions.keys():
        G.add_node((junction[0], junction[1], "junction"))  # "junction = (intron_start, intron_end)"
    # Add nodes for each exon
    for exon in exons.keys():
        G.add_node((exon[0] - 1, exon[1] + 1, "exon"))  # "exon = (exon_start-1, exon_end+1)"
    # Add edges between overlapping junctions (sharing start or end positions)
    junction_list = [(junction[0], junction[1], "junction") for junction in list(junctions.keys())]
    # need to connected exon and junction, but exon and junction have one base difference
    exon_list = [(exon[0] - 1, exon[1] + 1, "exon") for exon in list(exons.keys())]
    merged_list = junction_list + exon_list
    for i in range(len(merged_list)):
        for j in range(i + 1, len(merged_list)):
            start1, end1, type1 = merged_list[i]
            start2, end2, type2 = merged_list[j]
            # Check if they share a start or end position
            if type1 == type2:
                # junction-junction or exon-exon should share the donor or acceptor site
                if start1 == start2 or end1 == end2:
                    G.add_edge(merged_list[i], merged_list[j])
            else:
                # junction-exon or exon-junction should connect to each other
                if start1 == end2 or end1 == start2:
                    G.add_edge(merged_list[i], merged_list[j])
    connected_components = list(nx.connected_components(G))
    for component in connected_components:
        clu = []
        for node in component:
            if node[2] == "junction":
                gtag_tag = gt_ag_dict[(node[0], node[1])]
                clu.append((node[0], node[1], gtag_tag))
        if len(clu) > 0:
            junctions_clusters.append(clu)
    return junctions_clusters, junctions


def check_absent_present(start_pos, end_pos, reads_positions, reads_junctions):
    """
    Find the reads where an exon or junction is absent or present.
    :param start_pos: 1-based, start-inclusive
    :param end_pos: 1-based, end-inclusive
    :param reads_positions:
    :param reads_junctions:
    :return:
    """
    absent_reads = []
    present_reads = []
    for read_name, (read_start, read_end) in reads_positions.items():
        # if read_start > start_pos or read_end < end_pos:
        #     continue
        ## Based on Max Marin case, we need to use the following condition
        if read_start > end_pos or read_end < start_pos:
            continue
        present = False
        for junction_start, junction_end, _ in reads_junctions[read_name]:
            if junction_start == start_pos and junction_end == end_pos:
                present_reads.append(read_name)
                present = True
                break
        if not present:
            absent_reads.append(read_name)
    return absent_reads, present_reads


class AseEvent:
    def __init__(self, chr, start, end, novel, gt_ag_tag, gene_name, strand, junction_set, phase_set, hap1_absent,
                 hap1_present,
                 hap2_absent, hap2_present, h1_tissues_str, h2_tissues_str, tissues_ases_str, p_value_tsase,
                 p_value_ts):
        self.chr = chr
        self.start = start  # 1-based, inclusive
        self.end = end  # 1-based, inclusive
        self.strand = strand
        self.junction_set = junction_set
        self.phase_set = phase_set
        self.hap1_absent = hap1_absent
        self.hap1_present = hap1_present
        self.hap2_absent = hap2_absent
        self.hap2_present = hap2_present
        self.h1_tissues = h1_tissues_str
        self.h2_tissues = h2_tissues_str
        self.tissues_ases = tissues_ases_str
        self.p_value_tsase = p_value_tsase
        self.p_value_ts = p_value_ts
        self.novel = novel
        self.gt_ag_tag = gt_ag_tag
        self.gene_name = gene_name

    @staticmethod
    def __header__():
        return ("#Junction\tStrand\tJunction_set\tPhase_set\tHap1_tissues\tHap2_tissues\tTissues_ASEs\t"
                "P_value_tsase\tP_value_ts\tNovel\tGT_AG\tGene_name")

    def __str__(self):
        return (f"{self.chr}:{self.start}-{self.end}\t{self.strand}\t{self.junction_set}\t{self.phase_set}\t"
                f"{self.h1_tissues}\t{self.h2_tissues}\t{self.tissues_ases}\t"
                f"{self.p_value_tsase}\t{self.p_value_ts}\t{self.novel}\t{self.gt_ag_tag}\t{self.gene_name}")


def calc_sor(hap1_absent, hap1_present, hap2_absent, hap2_present):
    R = ((hap1_absent + 1) * (hap2_present + 1)) / ((hap1_present + 1) * (hap2_absent + 1))
    R_inverse = 1 / R
    sum = R + R_inverse
    SOR = math.log(sum)
    return SOR


def haplotype_event_test(absent_reads, present_reads, reads_tags, tissue_readname_map):
    """
    Perform Fisher's exact test to determine if the haplotype distribution is significantly different between absent and present reads.
    :param absent_reads:
    :param present_reads:
    :param reads_tags:
    :return:
    """
    hap_absent_counts = defaultdict(lambda: {1: [], 2: []})  # key: phase set, value: {hap1: readnames, hap2: readnames}
    hap_present_counts = defaultdict(lambda: {1: [], 2: []})
    for read_name in absent_reads:
        hap = reads_tags[read_name]["HP"]
        phase_set = reads_tags[read_name]["PS"]
        hap_absent_counts[phase_set][hap].append(read_name)
    for read_name in present_reads:
        hap = reads_tags[read_name]["HP"]
        phase_set = reads_tags[read_name]["PS"]
        hap_present_counts[phase_set][hap].append(read_name)
    all_phase_sets = set(hap_absent_counts.keys()).union(set(hap_present_counts.keys()))
    # get ps with the highest read coverage
    ps_read_count = {}
    for ps in all_phase_sets:
        h1_a, h2_a = len(hap_absent_counts[ps][1]), len(hap_absent_counts[ps][2])
        h1_p, h2_p = len(hap_present_counts[ps][1]), len(hap_present_counts[ps][2])
        ps_read_count[ps] = h1_a + h1_p + h2_a + h2_p
    if ps_read_count:
        most_reads_ps = sorted(ps_read_count.items(), key=lambda x: x[1], reverse=True)[0][0]
    else:
        return None
    phase_set = most_reads_ps
    # key: tissue, value: {"absent": count, "present": count}
    h1_count_table = defaultdict(lambda: {"absent": 0, "present": 0})
    h2_count_table = defaultdict(lambda: {"absent": 0, "present": 0})
    for readname in hap_absent_counts[phase_set][1]:
        tissue = tissue_readname_map[readname]
        h1_count_table[tissue]["absent"] += 1
        h1_count_table[tissue]["present"] += 0
        h2_count_table[tissue]["absent"] += 0
        h2_count_table[tissue]["present"] += 0
    for readname in hap_absent_counts[phase_set][2]:
        tissue = tissue_readname_map[readname]
        h2_count_table[tissue]["absent"] += 1
        h2_count_table[tissue]["present"] += 0
        h1_count_table[tissue]["absent"] += 0
        h1_count_table[tissue]["present"] += 0
    for readname in hap_present_counts[phase_set][1]:
        tissue = tissue_readname_map[readname]
        h1_count_table[tissue]["present"] += 1
        h1_count_table[tissue]["absent"] += 0
        h2_count_table[tissue]["present"] += 0
        h2_count_table[tissue]["absent"] += 0
    for readname in hap_present_counts[phase_set][2]:
        tissue = tissue_readname_map[readname]
        h2_count_table[tissue]["present"] += 1
        h2_count_table[tissue]["absent"] += 0
        h1_count_table[tissue]["present"] += 0
        h1_count_table[tissue]["absent"] += 0
    assert set(h1_count_table.keys()) == set(h2_count_table.keys())
    tissues = list(h1_count_table.keys())
    # table1 = []
    # for tissue in tissues:
    #     table1.append([[h1_count_table[tissue]["absent"] + 1, h1_count_table[tissue]["present"] + 1],
    #                    [h2_count_table[tissue]["absent"] + 1, h2_count_table[tissue]["present"] + 1]])
    # stratified_table = StratifiedTable(table1)
    # # whether there is an overall association between the row and column variables across all strata
    # cmh_overall = stratified_table.test_null_odds()
    # cmh_overall_p_value = cmh_overall.pvalue
    # cmh_overall_statistic = cmh_overall.statistic
    # # whether the odds ratios are consistent (homogeneous) across all strata
    # cmh_homogeneous = stratified_table.test_equal_odds()
    # cmh_homogeneous_p_value = cmh_homogeneous.pvalue
    # cmh_homogeneous_statistic = cmh_homogeneous.statistic

    # h1_cnts, h2_cnts = [], []
    # for tissue in tissues:
    #     h1_cnts.append(h1_count_table[tissue]["absent"])
    #     h1_cnts.append(h1_count_table[tissue]["present"])
    #     h2_cnts.append(h2_count_table[tissue]["absent"])
    #     h2_cnts.append(h2_count_table[tissue]["present"])
    # table = [[v + 1 for v in h1_cnts], [v + 1 for v in h2_cnts]]  # 2 x 2n table
    # _, p_value, _, _ = chi2_contingency(table)
    #
    # Tissue-specific allele-specific junction
    table2 = [[], [], [], []]
    for tissue in tissues:
        table2[0].append(h1_count_table[tissue]["absent"] + 1)
        table2[1].append(h1_count_table[tissue]["present"] + 1)
        table2[2].append(h2_count_table[tissue]["absent"] + 1)
        table2[3].append(h2_count_table[tissue]["present"] + 1)
    _, p_value_tsase, _, _ = chi2_contingency(table2)  # 4 x n table test

    # Tissue-specific junction
    table3 = [[], []]  # 2 x n table
    for tissue in tissues:
        table3[0].append(h1_count_table[tissue]["absent"] + h2_count_table[tissue]["absent"] + 1)  # absent
        table3[1].append(h1_count_table[tissue]["present"] + h2_count_table[tissue]["present"] + 1)  # present
    _, p_value_ts, _, _ = chi2_contingency(table3)  # 2 x n table test

    # separately call ASE for each tissue
    tissues_ases = {}
    for tissue in tissues:
        h1_absent = h1_count_table[tissue]["absent"]
        h1_present = h1_count_table[tissue]["present"]
        h2_absent = h2_count_table[tissue]["absent"]
        h2_present = h2_count_table[tissue]["present"]
        table = [[h1_absent, h2_absent], [h1_present, h2_present]]
        _, p_value = fisher_exact(table)
        sor = calc_sor(h1_absent, h1_present, h2_absent, h2_present)
        tissues_ases[tissue] = (p_value, sor)
    tissues_ases_str = []
    for tissue, (p_value, sor) in tissues_ases.items():
        tissues_ases_str.append(f"{tissue}:{p_value}:{sor}")
    tissues_ases_str = ",".join(tissues_ases_str)

    h1_tissue_str = []
    h2_tissue_str = []
    for tissue in tissues:
        h1_tissue_str.append(f'{tissue}:{h1_count_table[tissue]["absent"]}:{h1_count_table[tissue]["present"]}')
        h2_tissue_str.append(f'{tissue}:{h2_count_table[tissue]["absent"]}:{h2_count_table[tissue]["present"]}')
    h1_tissue_str = ",".join(h1_tissue_str)
    h2_tissue_str = ",".join(h2_tissue_str)
    event = (phase_set, p_value_tsase, p_value_ts,
             h1_tissue_str, h2_tissue_str, tissues_ases_str,
             len(hap_absent_counts[phase_set][1]), len(hap_present_counts[phase_set][1]),
             len(hap_absent_counts[phase_set][2]), len(hap_present_counts[phase_set][2]))
    return event


def analyze_gene(gene_name, gene_strand, annotation_exons, annotation_junctions, gene_region, reads_positions,
                 reads_tags, reads_exons, reads_introns, min_count, tissue_readname_map):
    chr = gene_region["chr"]
    start = gene_region["start"]
    end = gene_region["end"]
    gene_junction_set = set()
    for transcript_id, anno_junctions in annotation_junctions.items():
        for anno_junc in anno_junctions:
            gene_junction_set.add(anno_junc)
    gene_exon_set = set()
    for transcript_id, anno_exons in annotation_exons.items():
        for anno_exon in anno_exons:
            gene_exon_set.add(anno_exon)

    # # Extract relevant reads and regions
    # reads_positions, reads_exons, reads_introns, reads_tags = parse_reads_from_alignment(bam_file, chr, start, end)
    # junctions_clusters, read_junctions = cluster_junctions_connected_components(reads_introns, min_count)
    junctions_clusters, read_junctions = cluster_junctions_exons_connected_components(reads_introns, reads_exons,
                                                                                      min_count)

    # filter reads which have no overlapped exons with current gene exons
    intervalt = IntervalTree()
    for anno_exon in gene_exon_set:
        exon_start, exon_end = anno_exon[1:3]
        intervalt.addi(exon_start, exon_end + 1)  # interval is half-open, left inclusive, right exclusive

    reads_to_remove = []
    for qname, read_exons in reads_exons.items():
        overlapped = False
        for (exon_start, exon_end) in read_exons:
            if bool(intervalt.overlap(exon_start, exon_end + 1)):
                overlapped = True
                break
        if not overlapped:
            reads_to_remove.append(qname)

    # Remove reads after collecting them
    for qname in reads_to_remove:
        del reads_positions[qname]
        del reads_exons[qname]
        del reads_introns[qname]
        del reads_tags[qname]

    gene_ase_events = []
    # Analyze junction regions
    for junc_cluster in junctions_clusters:
        if len(junc_cluster) == 0:
            continue
        junction_set = f"{chr}:{junc_cluster[0][0]}-{junc_cluster[0][1]}"
        for read_junc in junc_cluster:
            junction_start = read_junc[0]
            junction_end = read_junc[1]
            gt_ag_tag = read_junc[2]
            novel = (chr, junction_start, junction_end) not in gene_junction_set
            # (extended_junction_start, extended_junction_end) = junctions_extended[(junction_start, junction_end)]
            # absences, presents = check_absent_present(extended_junction_start, extended_junction_end, reads_positions,
            #                                           reads_introns)
            absences, presents = check_absent_present(junction_start, junction_end, reads_positions, reads_introns)
            test_result = haplotype_event_test(absences, presents, reads_tags, tissue_readname_map)
            if test_result is None:
                continue
            (phase_set, pvalue_tsase, pvalue_ts, h1_tissue_str, h2_tissue_str, tissues_ases_str, h1_absent, h1_present,
             h2_absent, h2_present) = test_result
            gene_ase_events.append(AseEvent(chr, junction_start, junction_end, novel, gt_ag_tag, gene_name, gene_strand,
                                            junction_set, phase_set, h1_absent, h1_present, h2_absent, h2_present,
                                            h1_tissue_str, h2_tissue_str, tissues_ases_str,
                                            pvalue_tsase, pvalue_ts))
    return gene_ase_events


def analyze(annotation_file, bam_file, tissue_readnames, reference_file, output_prefix, min_count, gene_types, threads,
            no_gtag):
    tissue_readname_map = get_tissue_readnames(tissue_readnames)
    all_ase_events = {}  # key: (chr, start, end), value: AseEvent
    start_time = time.time()
    anno_gene_regions, anno_gene_names, anno_gene_strands, anno_exon_regions, anno_intron_regions = get_gene_regions(
        annotation_file, gene_types)
    print(f"Annotation file parsed in {time.time() - start_time:.2f} seconds")
    merged_genes_exons = merge_gene_exon_regions(anno_exon_regions)
    genome_dict = {}
    ref_genome = pysam.FastaFile(reference_file)
    for chrom in ref_genome.references:
        genome_dict[chrom] = ref_genome.fetch(chrom)
    start_time = time.time()
    read_assignment, reads_positions, reads_tags, reads_exons, reads_junctions = load_reads(bam_file,
                                                                                            genome_dict,
                                                                                            merged_genes_exons,
                                                                                            threads,
                                                                                            no_gtag)
    print(f"Reads assigned to genes in {time.time() - start_time:.2f} seconds")
    gene_assigned_reads = transform_read_assignment(read_assignment)
    gene_data_list = [(anno_gene_names[gene_id], anno_gene_strands[gene_id], anno_exon_regions[gene_id],
                       anno_intron_regions[gene_id], gene_region,
                       {name: reads_positions[name] for name in gene_assigned_reads[gene_id]},
                       {name: reads_tags[name] for name in gene_assigned_reads[gene_id]},
                       {name: reads_exons[name] for name in gene_assigned_reads[gene_id]},
                       {name: reads_junctions[name] for name in gene_assigned_reads[gene_id]},
                       min_count) for gene_id, gene_region in anno_gene_regions.items() if
                      (gene_region["chr"] in genome_dict) and (len(gene_assigned_reads[gene_id]) > 0)]
    print(f"Total genes to be analyzed: {len(gene_data_list)}")
    start_time = time.time()
    with Manager() as manager:
        shared_tissue_readname_map = manager.dict(tissue_readname_map)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(analyze_gene, *gene_data, shared_tissue_readname_map) for gene_data in
                       gene_data_list]
            for future in concurrent.futures.as_completed(futures):
                gene_ase_events = future.result()
                for event in gene_ase_events:
                    # multiple junctions in one gene_ase_events
                    key = (event.chr, event.start, event.end)
                    if key in all_ase_events.keys():
                        all_ase_events[key][event.gene_name] = event
                    else:
                        all_ase_events[key] = {event.gene_name: event}
    print(f"All Gene processed completed in {time.time() - start_time:.2f} seconds")

    # Perform multiple testing correction
    pass_idx = []  # index of junctions
    # p_values = []
    junctions = []
    for key in all_ase_events.keys():
        for gname in all_ase_events[key].keys():
            junctions.append((key, gname))  # key: (chr, start, end), gname: gene name
    print(f"Total junctions: {len(junctions)}")
    for idx in range(len(junctions)):
        junc = junctions[idx][0]
        gname = junctions[idx][1]
        event = all_ase_events[junc][gname]
        if event.hap1_absent + event.hap1_present + event.hap2_absent + event.hap2_present >= min_count:
            pass_idx.append(idx)
            # p_values.append(event.p_value)
    print(f"number of junctions with at least {min_count} reads: {len(pass_idx)}")
    # reject, adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    # asj_genes = {}
    # with open(output_prefix + ".joint_diff_splice.tsv", "w") as f:
    #     f.write(AseEvent.__header__() + "\n")
    #     for pi in range(len(pass_idx)):
    #         junc = junctions[pass_idx[pi]][0]
    #         gname = junctions[pass_idx[pi]][1]
    #         event = all_ase_events[junc][gname]
    #         event.p_value = adjusted_p_values[pi]
    #         f.write(event.__str__() + "\n")
    #         if gname not in asj_genes:
    #             asj_genes[gname] = [event.chr, event.p_value, event.p_value_tsase, event.p_value_ts, event.h1_tissues,
    #                                 event.h2_tissues]
    #         else:
    #             if event.p_value < asj_genes[gname][1]:
    #                 asj_genes[gname] = [event.chr, event.p_value, event.p_value_tsase, event.p_value_ts,
    #                                     event.h1_tissues, event.h2_tissues]
    # print(f"number of genes with allele-specific junctions: {len(asj_genes.keys())}")
    # with open(output_prefix + ".joint_asj_gene.tsv", "w") as f:
    #     f.write(f"#Gene_name\tChr\tP_value\tP_value_tsase\tP_value_ts\tHap1_tissues\tHap2_tissues\n")
    #     for gene_name in asj_genes:
    #         chr, p_value, p_value_tsase, p_value_ts, h1_tissues, h2_tissues = asj_genes[gene_name]
    #         f.write(f"{gene_name}\t{chr}\t{p_value}\t{p_value_tsase}\t{p_value_ts}\t{h1_tissues}\t{h2_tissues}\n")

    with open(output_prefix + ".joint_diff_splice.tsv", "w") as f:
        f.write(AseEvent.__header__() + "\n")
        for pi in range(len(pass_idx)):
            junc = junctions[pass_idx[pi]][0]
            gname = junctions[pass_idx[pi]][1]
            event = all_ase_events[junc][gname]
            f.write(event.__str__() + "\n")


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("-a", "--annotation_file", help="Annotation file in GFF3 or GTF format", required=True)
    parse.add_argument("-b", "--bam_file", help="BAM file", required=True)
    parse.add_argument("-f", "--reference", help="Reference genome file", required=True)
    parse.add_argument("-r", "--tissue_readnames", required=True,
                       help="File stores tissue name of each read in tsv format, read_name\ttissue_name")
    parse.add_argument("-o", "--output_prefix", help="prefix of output file", required=True)
    parse.add_argument("-t", "--threads", help="Number of threads", default=1, type=int)
    parse.add_argument("-g", "--gene_types", type=str, nargs="+", default=["protein_coding", "lncRNA"],
                       help='Gene types to be analyzed. Default is ["protein_coding", "lncRNA"]', )
    parse.add_argument("-m", "--min_sup", help="Minimum support of phased reads for exon or junction", default=10,
                       type=int)
    parse.add_argument("--no_gtag", action="store_true", help="Do not filter read junction with GT-AG signal")
    args = parse.parse_args()
    analyze(args.annotation_file, args.bam_file, args.tissue_readnames, args.reference, args.output_prefix, args.min_sup,
            args.gene_types, args.threads, args.no_gtag)
