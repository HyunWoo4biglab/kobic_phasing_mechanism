



def convert(vcf, sample_id, output):

    with open(output, 'w') as outfile:
        outfile.write('#chrom\tstart\tend\tkid_id\tvar_type\n')
        parser = md.VariantParser(vcf)
        parser.parse()
        v_l = parser.get_variant_list()
        for v in v_l:
            outfile.write(f"{v.contig}\t{v.vstart}\t{v.vend}\t{sample_id}\t{v.vtype}\n")

    return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for processing VCF and calculating performance ")
    parser.add_argument("-v", "--vcf", help = " input vcf file ", required = True)
    parser.add_argument("-i", "--sample_id", help = " input sample_id ", required = True)
    parser.add_argument("-o", "--output", help = " output bed ", required = True)

    args = parser.parse_args()

    import sys, os, pysam, pickle, re, math, resource
    from Bio.Seq import reverse_complement
    import numpy as np
    import dnsv_module as md
    import kmer_module as kmer_md

    convert(args.vcf, args.sample_id, args.output)
