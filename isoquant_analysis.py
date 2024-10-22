import argparse
import gzip
import pysam
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from multiprocessing import Manager
import concurrent.futures


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


def parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type):
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
    return records


def parse_isoquant_tpm(tpm_file):
    tpm_dict = {}
    with open(tpm_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            feature_id, tpm = parts[0], float(parts[1])
            tpm_dict[feature_id] = tpm
    return tpm_dict


# def get_gene_regions(annotation_file):
#     assert ".gff3" in annotation_file or ".gtf" in annotation_file, "Error: Unsupported annotation file format"
#     gene_regions = {}
#
#     def process_line(parts, gene_id):
#         chr, start, end = parts[0], int(parts[3]), int(parts[4])
#         gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}
#
#     def parse_attributes(attributes, key):
#         return [x for x in attributes if x.startswith(key)][0].split("=")[1]
#
#     if annotation_file.endswith(".gz"):
#         import gzip
#         with gzip.open(annotation_file, "rt") as f:
#             if ".gff3" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         gene_id = parse_attributes(attributes, "ID=")
#                         process_line(parts, gene_id)
#             elif ".gtf" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         tag, gene_id = attributes[0].split(" ")
#                         gene_id = gene_id.replace('"', '')
#                         assert tag == "gene_id"
#                         process_line(parts, gene_id)
#     else:
#         with open(annotation_file) as f:
#             if ".gff3" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         gene_id = parse_attributes(attributes, "ID=")
#                         process_line(parts, gene_id)
#             elif ".gtf" in annotation_file:
#                 for line in f:
#                     if line.startswith("#"):
#                         continue
#                     parts = line.strip().split("\t")
#                     if parts[2] == "gene":
#                         attributes = parts[8].split(";")
#                         tag, gene_id = attributes[0].split(" ")
#                         gene_id = gene_id.replace('"', '')
#                         assert tag == "gene_id"
#                         process_line(parts, gene_id)
#     return gene_regions


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


def get_reads_tag(bam_file, chr, start_pos, end_pos):
    reads_tag = {}
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for read in f.fetch(chr, start_pos, end_pos):
            PS = read.get_tag("PS") if read.has_tag("PS") else None
            HP = read.get_tag("HP") if read.has_tag("HP") else None
            reads_tag[read.query_name] = {"PS": PS, "HP": HP}
    return reads_tag


# def calculate_ase_pvalue(bam_file, gene_id, gene_name, gene_region, min_count, isoquant_read_assignments):
#     reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
#     assigned_reads = set()
#     for isoform_id, reads in isoquant_read_assignments[gene_id].items():
#         assigned_reads.update(reads)
#     phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
#     for rname in assigned_reads:
#         ps = reads_tag[rname]["PS"]
#         hp = reads_tag[rname]["HP"]
#         if ps and hp:
#             phase_set_hap_count[ps][hp] += 1
#     p_value = 1.0
#     h1_count, h2_count = 0, 0
#     phase_set = "."  # Placeholder for phase set
#     for ps, hap_count in phase_set_hap_count.items():
#         if hap_count[1] + hap_count[2] < min_count:
#             continue
#         # Calculate Binomial test p-value
#         total_reads = hap_count[1] + hap_count[2]
#         p_value_ase = binomtest(hap_count[1], total_reads, 0.5, alternative='two-sided').pvalue
#         if p_value_ase < p_value:
#             p_value = p_value_ase
#             h1_count = hap_count[1]
#             h2_count = hap_count[2]
#             phase_set = ps
#     return (gene_name, p_value, phase_set, h1_count, h2_count)


def calculate_ase_pvalue(bam_file, gene_id, gene_name, gene_region, min_count, isoquant_read_assignments):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set()

    # Collect assigned reads from isoquant_read_assignments
    for isoform_id, reads in isoquant_read_assignments[gene_id].items():
        assigned_reads.update(reads)

    # Track haplotype counts for each phase set
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
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

        return (gene_name, adjusted_p_value, most_significant_phase_set, h1_count, h2_count)

    # If no valid p-values are found, return defaults
    return (gene_name, 1.0, ".", 0, 0)


def calculate_ase_pvalue_pat_mat(bam_file, gene_id, gene_name, gene_region, min_count, isoquant_read_assignments,
                                 rna_vcfs, wg_vcfs):
    reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
    assigned_reads = set()

    # Collect assigned reads from isoquant_read_assignments
    for isoform_id, reads in isoquant_read_assignments[gene_id].items():
        assigned_reads.update(reads)

    # Track haplotype counts for each phase set
    phase_set_hap_count = defaultdict(lambda: {1: 0, 2: 0})  # key: phase set, value: {haplotype: count}
    for rname in assigned_reads:
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

        # determine paternal and maternal alleles for H1 and H2
        ps_variants = rna_vcfs.get(most_significant_phase_set, [])
        ps_reads = [rname for rname in assigned_reads if reads_tag[rname]["PS"] == most_significant_phase_set]
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
        return (gene_name, adjusted_p_value, most_significant_phase_set,
                h1_count, h2_count, h1_pat_cnt, h1_mat_cnt, h2_pat_cnt, h2_mat_cnt)

    # If no valid p-values are found, return defaults
    return (gene_name, 1.0, ".", 0, 0, 0, 0, 0, 0)


def analyze_ase_genes(assignment_file, annotation_file, bam_file, out_file, threads, gene_types, assignment_type,
                      classification_type, min_support):
    isoquant_read_assignments = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)

    gene_args = [(bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support)
                 for gene_id in gene_regions.keys()
                 if gene_id in isoquant_read_assignments]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(isoquant_read_assignments)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(calculate_ase_pvalue, *gene_data, shared_assignments) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    with open(out_file, "w") as f:
        f.write("Gene\tPS\tH1\tH2\tP-value\n")
        for gene_name, p_value, ps, h1, h2 in results:
            f.write(f"{gene_name}\t{ps}\t{h1}\t{h2}\t{p_value}\n")


def analyze_ase_genes_pat_mat(assignment_file, annotation_file, bam_file, vcf_file1, vcf_file2, out_file, threads,
                              gene_types, assignment_type, classification_type, min_support):
    isoquant_read_assignments = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    rna_vcfs = load_longcallR_phased_vcf(vcf_file1)
    wg_vcfs = load_whole_genome_phased_vcf(vcf_file2)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)

    gene_args = [(bam_file, gene_id, gene_names[gene_id], gene_regions[gene_id], min_support)
                 for gene_id in gene_regions.keys()
                 if gene_id in isoquant_read_assignments]
    results = []
    # Use a Manager to share isoquant_read_assignments across processes
    with Manager() as manager:
        shared_assignments = manager.dict(isoquant_read_assignments)
        shared_rna_vcfs = manager.dict(rna_vcfs)
        shared_wg_vcfs = manager.dict(wg_vcfs)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            # Submit all tasks at once without chunking
            futures = [executor.submit(calculate_ase_pvalue_pat_mat, *gene_data, shared_assignments, shared_rna_vcfs,
                                       shared_wg_vcfs) for gene_data in gene_args]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
    with open(out_file, "w") as f:
        f.write("Gene\tPS\tH1\tH2\tP-value\tH1_Paternal\tH1_Maternal\tH2_Paternal\tH2_Maternal\n")
        for gene_name, p_value, ps, h1, h2, h1_pat, h1_mat, h2_pat, h2_mat in results:
            f.write(f"{gene_name}\t{ps}\t{h1}\t{h2}\t{p_value}\t{h1_pat}\t{h1_mat}\t{h2_pat}\t{h2_mat}\n")


def analyze_tpm(assignment_file, annotation_file, bam_file, gene_types, assignment_type, classification_type,
                isoquant_gene_tpm,
                isoquant_transcript_tpm, output_gene_tpm, output_transcript_tpm):
    isoquant_read_assignment = parse_isoquant_read_assignment(assignment_file, assignment_type, classification_type)
    gene_regions, gene_names, gene_strands, exon_regions, intron_regions = get_gene_regions(annotation_file, gene_types)
    gene_tpm = parse_isoquant_tpm(isoquant_gene_tpm)
    transcript_tpm = parse_isoquant_tpm(isoquant_transcript_tpm)
    gene_tpm_writer = open(output_gene_tpm, "w")
    gene_tpm_writer.write("Gene\tTPM\tHap1_TPM\tHap2_TPM\n")
    transcript_tpm_writer = open(output_transcript_tpm, "w")
    transcript_tpm_writer.write("Isoform\tTPM\tHap1_TPM\tHap2_TPM\tHap1_reads\tHap2_reads\tP-value\n")
    for gene_id, isoform_records in isoquant_read_assignment.items():
        if gene_id not in gene_regions:
            continue
        gene_region = gene_regions[gene_id]
        reads_tag = get_reads_tag(bam_file, gene_region["chr"], gene_region["start"], gene_region["end"])
        gene_dict = defaultdict(int)
        for isoform_id, read_names in isoform_records.items():
            ps_dict = defaultdict(int)
            for rname in read_names:
                if rname in reads_tag:
                    read_info = reads_tag[rname]
                    if read_info['HP']:
                        ps_dict[read_info['HP']] += 1
                        gene_dict[read_info['HP']] += 1
            h1_cnt = ps_dict[1]
            h2_cnt = ps_dict[2]
            tpm = transcript_tpm.get(isoform_id, 0)
            if h1_cnt + h2_cnt == 0:
                h1_tpm, h2_tpm = 0, 0
                pvalue = 1.0
            else:
                h1_tpm = h1_cnt / (h1_cnt + h2_cnt) * tpm
                h2_tpm = h2_cnt / (h1_cnt + h2_cnt) * tpm
                pvalue = binomtest(h1_cnt, h1_cnt + h2_cnt, 0.5, alternative='two-sided').pvalue
            transcript_tpm_writer.write(f"{isoform_id}\t{tpm}\t{h1_tpm}\t{h2_tpm}\t{h1_cnt}\t{h2_cnt}\t{pvalue}\n")
        h1_cnt = gene_dict[1]
        h2_cnt = gene_dict[2]
        tpm = gene_tpm.get(gene_id, 0)
        if h1_cnt + h2_cnt == 0:
            h1_tpm, h2_tpm = 0, 0
        else:
            h1_tpm = h1_cnt / (h1_cnt + h2_cnt) * tpm
            h2_tpm = h2_cnt / (h1_cnt + h2_cnt) * tpm
        gene_tpm_writer.write(f"{gene_id}\t{tpm}\t{h1_tpm}\t{h2_tpm}\n")
    gene_tpm_writer.close()
    transcript_tpm_writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True, help="phased BAM file")
    parser.add_argument("--vcf1", required=False, help="longcallR phased vcf file", default=None)
    parser.add_argument("--vcf2", required=False, help="whole genome haplotype phased vcf file", default=None)
    parser.add_argument("-a", "--annotation", required=True, help="Annotation file")
    parser.add_argument("-i", "--assignment", required=True, help="Isoquant read assignment file")
    parser.add_argument("-g", "--gene_tpm", required=True, help="Isoquant gene TPM file")
    parser.add_argument("-t", "--transcript_tpm", required=True, help="Isoquant transcript TPM file")
    parser.add_argument("-o", "--output", required=True, help="prefix of output file")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Number of process to run")
    parser.add_argument("--gene_types", type=str, nargs="+", default=["protein_coding", "lncRNA"],
                        help='Gene types to be analyzed. Default is ["protein_coding", "lncRNA"]', )
    parser.add_argument('--assignment_type', type=str, nargs='+', default=["unique", "unique_minor_difference"],
                        help='Assignment types to include. Default is ["unique", "unique_minor_difference"].')
    parser.add_argument('--classification_type', type=str, nargs='+',
                        default=["full_splice_match", "incomplete_splice_match", "mono_exon_match"],
                        help='Classification types to include. Default is ["full_splice_match", "incomplete_splice_match", "mono_exon_match"].')
    parser.add_argument("--min_support", type=int, default=10,
                        help="Minimum support reads for counting event (default: 10)")

    args = parser.parse_args()

    gene_types = set(args.gene_types)
    assignment_type = set(args.assignment_type)
    classification_type = set(args.classification_type)

    if args.vcf1 and args.vcf2:
        analyze_ase_genes_pat_mat(args.assignment, args.annotation, args.bam, args.vcf1, args.vcf2,
                                  args.output + ".patmat.ase.tsv", args.processes, gene_types, assignment_type,
                                  classification_type, args.min_support)
    else:
        analyze_ase_genes(args.assignment, args.annotation, args.bam, args.output + ".ase.tsv", args.processes,
                          gene_types, assignment_type, classification_type, args.min_support)

    analyze_tpm(args.assignment, args.annotation, args.bam, gene_types, assignment_type, classification_type,
                args.gene_tpm, args.transcript_tpm, args.output + ".haplotype_gene_tpm.tsv",
                args.output + ".haplotype_transcript_tpm.tsv")
