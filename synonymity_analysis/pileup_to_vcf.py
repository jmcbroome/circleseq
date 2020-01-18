#!/usr/bin/env python3

#this script creates a custom 'vcf' without doing genotyping for circleseq output. Takes a samtools mpileup and a header text file (which can be obtained by creating a vcf output with mpileup, -uv, and grepping the # lines.)

#import
import argparse
import sys
import numpy as np
import statistics as st

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--header', help = 'File containing header text for a vcf of the given reference genome.')
    parser.add_argument('-p', '--pileup', help = 'Pileup to parse and force into a VCF format. Default is standard in', default = None)
    parser.add_argument('-o', '--output', help = 'Name of the vcf output file. Default is stdout', default = None)
    parser.add_argument('-g', '--germline', type = bool, help = 'Set to True to retain germline mutations. Default is False and ignores high frequency mutations', default = False)
    parser.add_argument('-d', '--mind', type = int, help = 'Minimum depth for somatic mutation identification. Default 10', default = 10)
    args = parser.parse_args()
    return args
#the following functions are intended to construct a tracking structure which can identify and ignore likely PCR duplicate errors
#these errors generally manifest as a series of alternative alleles which come from adjacently mapping consensus sequences, e.g. a line of alternative alleles will appear as "......AAAAAA......"
#in this case, we only want to count the single A error rather than counting it 5 times.
#we use a permuter which generates random sequences with an equal length and number of alternatives, measures median distance between each instances of the alternative allele, and determines whether a given read has an average density of alternative alleles which falls below this threshold
def get_dindex(altstring):
    distances = {}
    for i,base in enumerate(altstring):
        if base != '.':
            if base not in distances:
                last = i
                distances[base] = []
            else:
                distances[base].append(i-last)
                last = i
    dindex = {}
    for k,v in distances.items():
        if len(v)> 0:
            dindex[k] = st.median(v)
    return dindex

def make_random(length = 100, bases_to_use = 'A', num = 5):
    string = list('.' * length)
    for b in bases_to_use:
        locs = np.random.choice(length,num,replace = False)
        for l in locs:
            string[l] = b
    return ''.join(string)

def perm_index(leng = 100, num = 3, pnum = 1000):
    indeces = []
    for p in range(pnum):
        tstr = make_random(length = leng, num = num, bases_to_use='A')
        index = get_dindex(tstr)
        indeces.append(index['A'])
    return np.percentile(indeces,5)

def make_vcf_line(spent, mind = 10, germline = False, pcr_duplicate_track = {}):
    #convert a stripped and split mpileup line into a fake vcf line, filling in default values.
    #pileup: chrom loc ref depth vector_of_alts vector_of_quals
    #vcf: chrom loc ID(.) ref alt qual filter info(dp=...)
    #minimalist entry looking at the docs is just "DP" in info, qual is ., id is ., filter is PASS. So try those.
    chrom, loc, ref, ndepth, vector_of_alts, vector_of_quals = spent
    #real depth isn't depth, its the length of the vector of alts without other symbols or Ns, because most Ns are introduced by the consensus builder and not in the original reads.
    # pcr_duplicate_track = {} #using dynamic programming to save compute cycles for this qc measure
    depth = len([v for v in vector_of_alts if v in 'ACGT.'])
    if depth >= mind:
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
            # for b in 'ACGT': #for all possible bases
                # if quality_alts.count(b) > depth/4: #if that base is more than 25% of seen bases at this point
                    # quality_alts = [base for base in quality_alts if base != b] #remove it, it's almost certainly a germline mutation.
            #if I'm removing germline I can apply an additional filter which should remove PCR duplicates
            skip = '' #record no more than one of the bases that will be included here because of pcr duplicate inflation.
            dindeces = get_dindex(quality_alts)
            pcrc = [] #bases which appear to be in a pcr duplicate cluster get a copy here and then skipped by the main iterator
            for base in 'ACGT':
                basecount = quality_alts.count(base)
                if 2 <= basecount <= depth/4: #doesn't make sense to calculate for singletons, which I intend to skip by default now.
                    key = (len(quality_alts), basecount)
                    if key not in pcr_duplicate_track:
                        pcr_duplicate_track[key] = perm_index(leng = key[0], num = key[1])
                    thresh = pcr_duplicate_track[key]
                    if dindeces[base] < thresh: #less than 5% chance of getting a cluster like this. Lock this one to 1 instance
                        #print("QC: Base is skipped for clustering")
                        skip += base
                        #still count it once though
                        pcrc.append(base)
                    #if a base exists at appropriate levels but stays above the cluster threshold, it's collected in the fixed alts statement below
                else:
                    skip += base
        #collect individual bases not counted as a pcr cluster, plus a singleton of each cluster base
        fixed_alts = ','.join(pcrc + [q for q in quality_alts if q not in skip])
        # fixed_alts = ','.join(list(set(quality_alts))) #ignoring any weirdness around the mpileup, at least for now. Also, removing duplicates of somatic mutations.
        # print(len(fixed_alts), fixed_alts)
        #the assumption is also that germline variants have been cleaned out by pilon.
        if depth == 0 or len(fixed_alts) == 0: #nothing but Ns here.
            return None, pcr_duplicate_track
        else:
            vcf_line = chrom + '\t' + loc + '\t.\t' + ref + '\t' + fixed_alts + '\t.\tPASS\tDP=' + str(depth)
            return vcf_line, pcr_duplicate_track
    else:
        return None, pcr_duplicate_track

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
    pdt = {}
    for entry in pilein:
        nline, pdt = make_vcf_line(entry.strip().split(), args.mind, args.germline, pcr_duplicate_track=pdt)
        if nline != None:
            print(nline, file = outf)
    if args.pileup != None:
        pilein.close()
    if args.output != None:
        outf.close()

if __name__ == "__main__":
    main()