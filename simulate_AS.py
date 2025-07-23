import gzip
from collections import defaultdict
from intervaltree import IntervalTree, Interval
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def generate_simulated_transcripts(annotation_file, reference_genome, output_transcripts, gene_types, asj_case_ratio,
                                   snp_rate):
    """
    1. Pick two longest transcripts with overlapped regions for a given gene.
    2. for two transcripts, randomly apply variants to the overlapped regions for ASJ case.
    3. Do nothing for non-ASJ case.
    """
    assert annotation_file.endswith((".gff3", ".gtf", ".gff3.gz", ".gtf.gz")), "Error: Unknown annotation file format"
    gene_regions = {}
    gene_names = {}
    exon_regions = defaultdict(lambda: defaultdict(list))

    def process_gene(parts, gene_id, gene_name):
        chr, start, end = parts[0], int(parts[3]), int(parts[4])
        gene_regions[gene_id] = {"chr": chr, "start": start, "end": end}  # 1-based, start-inclusive, end-inclusive
        gene_names[gene_id] = gene_name

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

    # Load reference genome
    ref_seqs = {rec.id: rec.seq for rec in SeqIO.parse(reference_genome, "fasta")}

    isolated_genes = set()
    for contig in ref_seqs.keys():
        # filter out genes with overlaps
        gene_ivtree = IntervalTree()
        for gene_id, region in gene_regions.items():
            chr = region["chr"]
            start = region["start"]
            end = region["end"]
            if chr!=contig:
                continue
            iv = Interval(start-1, end, gene_id)  # Convert to 0-based, start-inclusive, end-exclusive
            gene_ivtree.add(iv)
        for iv in gene_ivtree:
            gene_ivtree.remove(iv)
            overlaps = gene_ivtree.overlap(iv.begin, iv.end)
            gene_ivtree.add(iv)  # Re-add the original interval
            if not overlaps:
                isolated_genes.add(iv.data)

    # Filter out genes with overlaps
    print(f"Total genes before filtering: {len(gene_regions)}")
    exon_regions = {gene_id: exons for gene_id, exons in exon_regions.items() if gene_id in isolated_genes}
    print(f"Total genes after filtering: {len(exon_regions)}")

    # Simulate transcripts
    simulated_records = []
    ASJ_case_count = 0
    non_ASJ_case_count = 0
    for gene_id, exons in exon_regions.items():
        gene_name = gene_names.get(gene_id, None)
        if not gene_name:
            continue

        sorted_transcripts = sorted(
            exons.items(),
            key=lambda x: sum(e - s + 1 for _, s, e in x[1] if len(x[1]) >= 2),
            reverse=True
        )

        if len(sorted_transcripts) < 2:
            continue

        transcripts = sorted_transcripts[:2]
        t1_id, t1_exons = transcripts[0]
        t2_id, t2_exons = transcripts[1]

        # sort t1_exons and t2_exons, since the transcript may be reverse-stranded
        t1_exons = sorted(t1_exons, key=lambda x: x[1])
        t2_exons = sorted(t2_exons, key=lambda x: x[1])

        tree = IntervalTree()
        for _, s, e in t1_exons:
            if s < e:
                tree[s:e] = True
            else:
                # Invalid exon, e.g., chr7    HAVANA  exon    110511060       110511060
                tree = IntervalTree()  # Reset tree if invalid exon
                break

        overlapped_regions = []
        for _, s, e in t2_exons:
            for iv in tree.overlap(s, e):
                ov_start = max(s, iv.begin)
                ov_end = min(e, iv.end)
                if ov_start < ov_end:
                    overlapped_regions.append((ov_start, ov_end))

        if not overlapped_regions:
            continue

        # Simulate transcript sequences:
        def extract_seq(exon_list):
            seq = ''
            for chr, s, e in exon_list:
                subseq = ref_seqs[chr][s - 1:e]  # FASTA is 1-based inclusive
                seq += str(subseq)
            return seq

        seq1 = extract_seq(t1_exons)
        seq2 = extract_seq(t2_exons)

        # Randomly apply variants for ASJ case (simulate SNPs)
        def apply_random_snps(seq, exon_list, overlapped_regions, snp_rate=0.05):
            seq_list = list(seq)
            offset = 0
            variant_count = 0
            variant_positions = []
            for chr, s, e in exon_list:
                for ov_s, ov_e in overlapped_regions:
                    if s <= ov_s <= e or s <= ov_e <= e:
                        for pos in range(max(s, ov_s), min(e, ov_e) + 1):
                            idx = offset + (pos - s)
                            if random.random() < snp_rate:  # 5% chance to introduce SNP
                                ref_allele = seq_list[idx]
                                if seq_list[idx] == "A":
                                    seq_list[idx] = random.choice(['C', 'T'])
                                elif seq_list[idx] == "C":
                                    seq_list[idx] = random.choice(['A', 'G', 'T'])
                                elif seq_list[idx] == "G":
                                    seq_list[idx] = random.choice(['A', 'C', 'T'])
                                elif seq_list[idx] == "T":
                                    seq_list[idx] = random.choice(['A', 'G'])
                                alt_allele = seq_list[idx]
                                variant_positions.append((chr, pos, ref_allele, alt_allele))
                                variant_count += 1
                offset += e - s + 1
            return ''.join(seq_list), variant_count, variant_positions

        # Randomly apply variants together for both transcripts
        def apply_random_snps_dual(seq1, seq2, exon_list1, exon_list2, overlapped_regions, snp_rate=0.05):

            # Collect overlap positions
            overlap_positions = set()
            for ov_s, ov_e in overlapped_regions:
                overlap_positions.update(range(ov_s, ov_e + 1))

            # Determine SNP positions and their replacement base (same for both transcripts)
            snp_info = {}  # {genomic_pos: new_base}
            for pos in overlap_positions:
                if random.random() < snp_rate:
                    snp_info[pos] = None  # Placeholder for replacement base to determine later

            variant_count = len(snp_info)
            variant_positions = []

            # Determine replacement base for each SNP position (same alt base for both transcripts)
            # Helper to pick an alternate base
            def pick_alt_base(base):
                if base == "A":
                    return random.choice(['C', 'T'])
                elif base == "C":
                    return random.choice(['A', 'G', 'T'])
                elif base == "G":
                    return random.choice(['A', 'C', 'T'])
                elif base == "T":
                    return random.choice(['A', 'G'])
                else:
                    return base  # Leave non-ACGT unchanged just in case

            # Helper to mutate a single transcript sequence
            def mutate(seq, exon_list, snp_info):
                seq_list = list(seq)
                transcript_pos = 0
                for chr, s, e in exon_list:
                    for pos in range(s, e + 1):
                        if pos in snp_info:

                            if snp_info[pos] is None:
                                ref_allele = seq_list[transcript_pos]
                                # Determine alternate base using the original base in this transcript
                                alt_base = pick_alt_base(seq_list[transcript_pos])
                                snp_info[pos] = alt_base
                                variant_positions.append((chr, pos, ref_allele, alt_base))
                            else:
                                alt_base = snp_info[pos]

                            seq_list[transcript_pos] = alt_base
                        transcript_pos += 1
                return ''.join(seq_list)

            # Mutate both transcripts using the same snp_info
            mutated_seq1 = mutate(seq1, exon_list1, snp_info)
            mutated_seq2 = mutate(seq2, exon_list2, snp_info)

            return mutated_seq1, mutated_seq2, variant_count, variant_positions

        if random.random() < asj_case_ratio:
            # ASJ case, only apply SNPs to the first transcript
            seq1_mut, variant_count1, variant_positions1 = apply_random_snps(seq1, t1_exons, overlapped_regions,
                                                                             snp_rate)
            seq2_mut = seq2
            if variant_count1 > 0:
                rec1 = SeqRecord(Seq(seq1_mut), id=f"{gene_name}:{gene_id}:{t1_id}:ASJ",
                                 description=f"Gene: {gene_name}, Variants Count: {variant_count1}, Variant Positions: {variant_positions1}")
                rec2 = SeqRecord(Seq(seq2_mut), id=f"{gene_name}:{gene_id}:{t2_id}:ASJ",
                                 description=f"Gene: {gene_name}")
                ASJ_case_count += 1
            else:
                rec1 = SeqRecord(Seq(seq1_mut), id=f"{gene_name}:{gene_id}:{t1_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                rec2 = SeqRecord(Seq(seq2_mut), id=f"{gene_name}:{gene_id}:{t2_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                non_ASJ_case_count += 1
            simulated_records.extend([rec1, rec2])
        else:
            # Non-ASJ case, both transcripts apply the same SNPs + original transcripts (4 transcripts)
            seq1_no_mut = seq1
            seq2_no_mut = seq2
            seq1_mut, seq2_mut, variant_count1, variant_positions1 = apply_random_snps_dual(seq1, seq2, t1_exons,
                                                                                            t2_exons,
                                                                                            overlapped_regions,
                                                                                            snp_rate)
            if variant_count1 > 0:
                rec1 = SeqRecord(Seq(seq1_no_mut), id=f"{gene_name}:{gene_id}:{t1_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                rec2 = SeqRecord(Seq(seq2_no_mut), id=f"{gene_name}:{gene_id}:{t2_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                rec3 = SeqRecord(Seq(seq1_mut), id=f"{gene_name}:{gene_id}:{t1_id}:non-ASJ",
                                 description=f"Gene: {gene_name}, Variants Count: {variant_count1}, Variant Positions: {variant_positions1}")
                rec4 = SeqRecord(Seq(seq2_mut), id=f"{gene_name}:{gene_id}:{t2_id}:non-ASJ",
                                 description=f"Gene: {gene_name}, Variants Count: {variant_count1}, Variant Positions: {variant_positions1}")
                non_ASJ_case_count += 1
                simulated_records.extend([rec1, rec2, rec3, rec4])
            else:
                rec1 = SeqRecord(Seq(seq1_no_mut), id=f"{gene_name}:{gene_id}:{t1_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                rec2 = SeqRecord(Seq(seq2_no_mut), id=f"{gene_name}:{gene_id}:{t2_id}:non-ASJ",
                                 description=f"Gene: {gene_name}")
                non_ASJ_case_count += 1
                simulated_records.extend([rec1, rec2])

    # Write simulated transcripts
    with open(output_transcripts, "w") as out_f:
        SeqIO.write(simulated_records, out_f, "fasta")

    print(
        f"Total ASJ cases: {ASJ_case_count}, Non-ASJ cases: {non_ASJ_case_count}, Allele-specific ratio: {asj_case_ratio}, SNP rate: {snp_rate}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate allele-specific junctions (ASJ) in transcripts.")
    parser.add_argument("-a", "--annotation", type=str, required=True, help="Annotation file (GFF3 or GTF)")
    parser.add_argument("-r", "--reference", type=str, required=True, help="Reference genome FASTA file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output FASTA file for simulated transcripts")
    parser.add_argument("-g", "--gene_types", type=str, nargs='+', default=["protein_coding"],
                        help="List of gene types to consider (default: protein_coding)")
    parser.add_argument("-c", "--asj_case_ratio", type=float, default=0.5,
                        help="Ratio of ASJ cases to non-ASJ cases (default: 0.5)")
    parser.add_argument("-s", "--snp_rate", type=float, default=0.01,
                        help="SNP rate for ASJ case (default: 0.01)")
    args = parser.parse_args()

    generate_simulated_transcripts(args.annotation, args.reference, args.output, args.gene_types,
                                   args.asj_case_ratio, args.snp_rate)