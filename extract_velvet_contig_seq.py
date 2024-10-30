

def extract_contig_seq(bam, output, target_region):

    sam = pysam.AlignmentFile(bam, 'rb')
#    test_region = 'chr9:18032624-18033021'
#    contig = 'chr9'
#    start = 18027700; end = 18037947
    contig = target_region.split(':')[0]
    start = int(target_region.split(':')[-1].split('-')[0]) - 5000
    end = int(target_region.split(':')[-1].split('-')[-1]) + 5000

    for read in sam.fetch(contig, start, end):
        print(read.query_name, read.cigarstring)

        cigar = read.cigarstring
        cigar_match = re.findall(r'(\d+)([A-Z])', cigar)
        if len(cigar_match) > 2:
            bp1_pos = int(cigar_match[0][0]); bp2_pos = int(cigar_match[-1][0])
            print(cigar_match)
            try:

                bp1_seq = read.query_sequence[:bp1_pos]
                bp2_seq = read.query_sequence[-bp2_pos:]

                bp1_out = output + '.bp1.fasta'
                bp2_out = output + '.bp2.fasta'
                with open(bp1_out, 'w') as infile:
                    infile.write(f'>bp1\n{bp1_seq}\n')

                with open(bp2_out, 'w') as infile:
                    infile.write(f'>bp2\n{bp2_seq}\n')
            except Exception as e:
                print('Exception occurred', e)
                print(read.query_sequence)
                pass
            finally:
                pass

    sam.close()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for processing VCF and calculating performance ")
    parser.add_argument("-i", "--bam", help = " input child bam file ", required = True)
    parser.add_argument("-r", "--region", help = " input target region ", required = True)
    parser.add_argument("-o", "--output", help = " output file prefix ", required = True)
    args = parser.parse_args()

    import re, pysam

    extract_contig_seq(args.bam, args.output, args.region)
