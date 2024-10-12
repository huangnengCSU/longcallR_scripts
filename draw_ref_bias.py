import matplotlib.pyplot as plt
import pysam


def draw_ref_bias(vcf_file):
    alt_allele_fractions = []
    ref_alleles, alt_alleles = 0, 0
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            if record.filter.keys() != ['PASS']:
                continue
            ref = record.ref
            alt = record.alts[0]

            # get genotype information
            gt = record.samples[0]['GT']
            # Check for heterozygous variants
            if gt == (0, 1) or gt == (1, 0):
                # Get the allele fraction of the alternate allele
                af = record.samples[0]['AF'][0]
                alt_allele_fractions.append(af)
                dp = record.samples[0]['DP']
                if dp == 0:
                    continue
                # Get the allele counts
                alt_alleles += int(af * dp)
                ref_alleles += int((1 - af) * dp)
    print(f"R= {ref_alleles}, A= {alt_alleles}")
    print(f"R/(R+A)= {ref_alleles / (ref_alleles + alt_alleles)}")

    # plot the distribution of alternate allele fractions
    plt.hist(alt_allele_fractions, bins=100, alpha=0.75, color='blue', edgecolor='black')
    plt.title('Distribution of Alternate Allele Fractions')
    plt.xlabel('Alternate Allele Fraction')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
