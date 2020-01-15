#!/usr/bin/env python3

#this script creates a custom 'vcf' without doing genotyping for circleseq output. Takes a samtools mpileup and a header text file (which can be obtained by creating a vcf output with mpileup, -uv, and grepping the # lines.)

#import
import argparse
import sys
#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--header', help = 'File containing header text for a vcf of the given reference genome.')
    parser.add_argument('-p', '--pileup', help = 'Pileup to parse and force into a VCF format. Default is standard in', default = None)
    parser.add_argument('-o', '--output', help = 'Name of the vcf output file. Default is stdout', default = None)
    parser.add_argument('-g', '--germline', type = bool, help = 'Set to True to retain germline mutations. Default is False and ignores high frequency mutations', default = False)
    args = parser.parse_args()
    return args

def make_vcf_line(spent, germline = False):
    #convert a stripped and split mpileup line into a fake vcf line, filling in default values.
    #pileup: chrom loc ref depth vector_of_alts vector_of_quals
    #vcf: chrom loc ID(.) ref alt qual filter info(dp=...)
    #minimalist entry looking at the docs is just "DP" in info, qual is ., id is ., filter is PASS. So try those.
    chrom, loc, ref, ndepth, vector_of_alts, vector_of_quals = spent
    #real depth isn't depth, its the length of the vector of alts without other symbols or Ns, because most Ns are introduced by the consensus builder and not in the original reads.
    depth = len([v for v in vector_of_alts if v in 'ACGT.'])
    quality_alts = []
    cleaner = [b for b in vector_of_alts if b in 'ACGTN.']
    assert len(cleaner) == len(vector_of_quals)
    for i,b in enumerate(cleaner):
        if b in 'ACGT':
            #check if the qual value is high enough.
            if int(vector_of_quals[i]) > 1:
                quality_alts.append(b)
    #in quality alts, the majority or entirety of the set may all be the same base, which happens when its a germline mutation.
    #note that I can't distinguish germline from somatic at low depths, but with Wri datasets at higher ones I can.
    if not germline:
        for b in 'ACGT': #for all possible bases
            if quality_alts.count(b) > depth/4: #if that base is more than 25% of seen bases at this point
                quality_alts = [base for base in quality_alts if base != b] #remove it, it's almost certainly a germline mutation.

    fixed_alts = ','.join(list(set(quality_alts))) #ignoring any weirdness around the mpileup, at least for now. Also, removing duplicates of somatic mutations.
    # print(len(fixed_alts), fixed_alts)
    #the assumption is also that germline variants have been cleaned out by pilon.
    if depth == 0 or len(fixed_alts) == 0: #nothing but Ns here.
        return None
    else:
        vcf_line = chrom + '\t' + loc + '\t.\t' + ref + '\t' + fixed_alts + '\t.\tPASS\tDP=' + str(depth)
        return vcf_line

def main():
    args = argparser()
    #insert code
    if args.pileup == None:
        pilein = sys.stdin
    else:
        pilein = open(args.pileup)
    if args.output == None:
        outf = sys.stdout
    else:
        outf = open(args.output, 'w+')
    # with open(args.output, 'w+') as outf:
    with open(args.header) as tin:
        for entry in tin:
            print(entry.strip(), file = outf)
    for entry in pilein:
        nline = make_vcf_line(entry.strip().split(), args.germline)
        if nline != None:
            print(nline, file = outf)
    if args.pileup != None:
        pilein.close()
    if args.output != None:
        outf.close()

if __name__ == "__main__":
    main()