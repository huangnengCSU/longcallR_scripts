import pysam
from collections import defaultdict
import gzip
from intervaltree import Interval, IntervalTree
import bisect


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

    # read_assignment, key: read_name, value: directory of gene_id: overlap_length
    read_assignment = defaultdict(lambda: defaultdict(int))

    trees = defaultdict(IntervalTree)  # key: chr, value: IntervalTree
    gene_intervals = defaultdict(lambda: defaultdict(IntervalTree))  # key: chr, value: dict of gene_id: IntervalTree
    for chrom in merged_genes_exons.keys():
        for gene_id, merged_exons in merged_genes_exons[chrom].items():
            gene_region = (merged_exons[0][0], merged_exons[-1][1])  # 1-based, start-inclusive, end-inclusive
            trees[chrom].add(Interval(gene_region[0], gene_region[1] + 1), gene_id)

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

            # read_overlap_length = {}  # key: gene_id, value: overlap_length
            # calculate the overlap of splice regions with gene exons to determine assignment of read to gene
            for gene_id in candidate_gene_ids:
                if gene_id not in gene_intervals[chromosome]:
                    continue
                overlap_length = 0
                for splice_region in splice_regions:
                    overlap_length += sum(
                        interval.end - interval.begin
                        for interval in gene_intervals[chromosome][gene_id].overlap(*splice_region)
                    )
                # read_overlap_length[gene_id] = overlap_length
                read_assignment[read.query_name][gene_id] = overlap_length
    return read_assignment
