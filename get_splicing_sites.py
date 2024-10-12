import gzip
import argparse

class Exon:
    def __init__(self, chr, start, end, strand, transcript_id) -> None:
        self.chr = chr
        self.start = start  # 1-based, inclusive
        self.end = end  # 1-based, inclusive
        self.strand = strand
        self.transcript_id = transcript_id

class DonorSite:
    def __init__(self, chr, pos, strand, transcript_id) -> None:
        self.chr = chr
        self.pos = pos
        self.strand = strand
        self.transcript_id = transcript_id

class AcceptorSite:
    def __init__(self, chr, pos, strand, transcript_id) -> None:
        self.chr = chr
        self.pos = pos
        self.strand = strand
        self.transcript_id = transcript_id
        

def parse_gff(gff_file):
    gz_flag = False
    if gff_file.endswith('.gz'):
        fopen = gzip.open(gff_file, 'rb')
        gz_flag = True
    else:
        fopen = open(gff_file, 'r')

    donor_sites = []
    acceptor_sites = []
    exon_list = []
    transcript_id = ''

    for line in fopen:
        if gz_flag:
            line = line.decode('utf-8')
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chr = fields[0]
        type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        if type == 'transcript':
            transcript_id = fields[8].strip().split(';')[0].replace('ID=', '')
            if len(exon_list) == 0:
                continue
            if exon_list[0].strand == '+':
                ## forward strand, exon is present from left to right
                for i, exon in enumerate(exon_list):
                    if len(exon_list) == 1:
                        ## single exon transcript
                        continue
                    if i==0:
                        # the first exon
                        dn = DonorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        donor_sites.append(dn)
                    if i==len(exon_list)-1:
                        # the last exon
                        ac = AcceptorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        acceptor_sites.append(ac)
                    if i>0 and i<len(exon_list)-1:
                        dn = DonorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        ac = AcceptorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        donor_sites.append(dn)
                        acceptor_sites.append(ac)
            else:
                ## reverse strand, exon is present from right to left, donor site is 1-base before the start of exon, acceptor site is 1-base after the end of exon
                for i, exon in enumerate(exon_list):
                    if len(exon_list) == 1:
                        ## single exon transcript
                        continue
                    if i==0:
                        # the first exon
                        dn = DonorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        donor_sites.append(dn)
                    if i==len(exon_list)-1:
                        # the last exon
                        ac = AcceptorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        acceptor_sites.append(ac)
                    if i>0 and i<len(exon_list)-1:
                        ac = AcceptorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        dn = DonorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        acceptor_sites.append(ac)
                        donor_sites.append(dn)

            exon_list = []

        if type == 'exon':
            exon = Exon(chr, start, end, strand, transcript_id)
            exon_list.append(exon)
    
    # process the last transcript
    if len(exon_list) > 0:
        if exon_list[0].strand == '+':
                ## forward strand, exon is present from left to right
                for i, exon in enumerate(exon_list):
                    if len(exon_list) == 1:
                        ## single exon transcript
                        continue
                    if i==0:
                        # the first exon
                        dn = DonorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        donor_sites.append(dn)
                    if i==len(exon_list)-1:
                        # the last exon
                        ac = AcceptorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        acceptor_sites.append(ac)
                    if i>0 and i<len(exon_list)-1:
                        dn = DonorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                        ac = AcceptorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                        donor_sites.append(dn)
                        acceptor_sites.append(ac)
        else:
            ## reverse strand, exon is present from right to left, donor site is 1-base before the start of exon, acceptor site is 1-base after the end of exon
            for i, exon in enumerate(exon_list):
                if len(exon_list) == 1:
                        ## single exon transcript
                        continue
                if i==0:
                    # the first exon
                    dn = DonorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                    donor_sites.append(dn)
                if i==len(exon_list)-1:
                    # the last exon
                    ac = AcceptorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                    acceptor_sites.append(ac)
                if i>0 and i<len(exon_list)-1:
                    ac = AcceptorSite(exon.chr, exon.end+1, exon.strand, exon.transcript_id)
                    dn = DonorSite(exon.chr, exon.start-1, exon.strand, exon.transcript_id)
                    acceptor_sites.append(ac)
                    donor_sites.append(dn)
    fopen.close()
    return donor_sites, acceptor_sites


def main():
    parse = argparse.ArgumentParser(description="Get donor and acceptor sites from gff file")
    parse.add_argument('-i', help='GFF file')
    args = parse.parse_args()
    gff_file = args.i
    donor_sites, acceptor_sites = parse_gff(gff_file)
    print('Donor sites:', len(donor_sites))
    print('Acceptor sites:', len(acceptor_sites))

    donor_site_set = set()
    for site in donor_sites:
        ## remove transcript id, because same donor site may belong to different transcripts
        donor_site_set.add((site.chr, site.pos, site.strand))
    acceptor_site_set = set()
    for site in acceptor_sites:
        ## remove transcript id, because same acceptor site may belong to different transcripts
        acceptor_site_set.add((site.chr, site.pos, site.strand))

    print('Remove duplicate donor sites:', len(donor_site_set))
    print('Remove duplicate acceptor sites:', len(acceptor_site_set))

    # for site in donor_site_set:
    #     print(site[0], site[1], site[2], sep='\t')

    # for site in acceptor_site_set:
    #     print(site[0], site[1], site[2], sep='\t')

if __name__=="__main__":
    main()





