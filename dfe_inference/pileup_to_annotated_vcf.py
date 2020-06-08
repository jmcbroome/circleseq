#!/usr/bin/env python3

#this script creates a custom 'vcf' without doing genotyping for circleseq output. Takes a samtools mpileup and a header text file (which can be obtained by creating a vcf output with mpileup, -uv, and grepping the # lines.)
#this alternative version specifically formats the output vcf to follow the specifications required for the software SNPGenie, using --vcfformat=2 .
#import
import argparse
import sys
import numpy as np
import statistics as st
from scipy.stats import beta, binom, betabinom
from scipy.optimize import fmin
import subprocess
import math

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--header', help = 'File containing header text for a vcf of the given reference genome.')
    parser.add_argument('-p', '--pileup', help = 'Pileup to parse and force into a VCF format. Default is standard in', default = None)
    parser.add_argument('-o', '--output', help = 'Name of the vcf output file. Default is stdout', default = None)
    parser.add_argument('-g', '--germline', type = bool, help = 'Set to True to retain germline mutations. Default is False and ignores high frequency mutations', default = False)
    #parser.add_argument('-f', '--freq_cutoff', type = float, help = 'Value between 0 and 1, above which all entries are ignored. Default 1.')
    parser.add_argument('-d', '--mind', type = int, help = 'Minimum depth for somatic mutation identification. Default 10', default = 10)
    parser.add_argument('-r', '--params', nargs = 2, type = float, help = 'Alpha and beta parameters for the betabinomial')
    parser.add_argument('-b', '--bayes', type = bool, help = 'Set to True to use the Bayesian Posterior Mean Estimator for most likely frequency. False instead uses Nelder-Mead MLE. Default True', default = True)
    args = parser.parse_args()
    return args

ldf = {} #update this with likelihoods to save on redundant calculations

def get_fp(f, a, b, d, x):
    return -1 * beta.pdf(a=a,b=b,x=f) * binom.pmf(n=d,k=x,p=f)

def add_likelihood_info_bayes(nline, a, b):
    '''
    Adds additional fields to a vcf info line describing the likelihood of producing the site and predicted frequency, with and without testing correction. Uses the bayesian beta posterior mean estimator.
    '''
    chro, pos, rid, ref, alt, qual, rfil, info = nline.strip().split()
    altset = alt.split(',')
    #need alt count and depth value from info
    iv = info.split(';')
    dp = int(iv[0].strip("DP="))
    arats = [float(v) for v in iv[1].strip('AF=').split(',')] 
    if len(altset) != len(arats):
        print('Mismatch', nline)

    ninfo = ''
    for i,ar in enumerate(arats):
        alt_count = math.floor(ar * dp) #should always be intable.
        #the structure of the new information will go as follows:
        #for each alt, there will be a series of comma delineated entries
        #that represent the most likely frequency and overall likelihood of the site
        #this is uncorrected so don't go thinking things are so super special or anything.
        #if I have one mutation that's an A, I'll add ";A=.003,.000001" for example.
        likelihood = betabinom.pmf(n=dp,a=a,b=b,k=alt_count)
        #calculate the most likely frequency. 
        #using the simple bayesian beta posterior mean equation
        #phat = (x + a) / (n + a + b)
        lf = (alt_count + a) / (dp + a + b)
        #add information to the line.
        ninfo += ';' + altset[i] + '=f:' + str(lf) + ',l:' + str(likelihood)
    return '\t'.join([chro, pos, rid, ref, alt, qual, rfil, info + ninfo])

def add_likelihood_info_mle(nline, a, b):
    '''
    Adds additional fields to a vcf info line describing the likelihood of producing the site, with and without testing correction. Uses an MLE approach.
    '''
    chro, pos, rid, ref, alt, qual, rfil, info = nline.strip().split()
    altset = alt.split(',')
    #need alt count and depth value from info
    iv = info.split(';')
    dp = int(iv[0].strip("DP="))
    arats = [float(v) for v in iv[1].strip('AF=').split(',')] 
    if len(altset) != len(arats):
        print('Mismatch', nline)

    ninfo = ''
    for i,ar in enumerate(arats):
        alt_count = math.floor(ar * dp) #should always be intable.
        #the structure of the new information will go as follows:
        #for each alt, there will be a series of comma delineated entries
        #that represent the most likely frequency and overall likelihood of the site
        #this is uncorrected so don't go thinking things are so super special or anything.
        #if I have one mutation that's an A, I'll add ";A=.003,.000001" for example.
        likelihood = betabinom.pmf(n=dp,a=a,b=b,k=alt_count)
        #calculate the most likely frequency. 
        #using dynamic programming to save on runtime.
        if (alt_count,dp) in ldf:
            lf = ldf[(alt_count,dp)]
        else:
            lf = fmin(func = get_fp, x0 = .01, args = (a,b,dp,alt_count), disp = False)
            ldf[(alt_count,dp)] = lf[0]
        #add information to the line.
        ninfo += ';' + altset[i] + '=f:' + str(lf) + ',l:' + str(likelihood)
    return '\t'.join([chro, pos, rid, ref, alt, qual, rfil, info + ninfo])


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
    basedepth = len([v for v in vector_of_alts if v in 'ACGT.'])
    if basedepth >= mind and ref != 'N':
        depth = 0
        altbase_counts = {}
        quality_alts = []
        cleaner = [b for b in vector_of_alts if b in 'ACGTN.']
        assert len(cleaner) == len(vector_of_quals)
        for i,b in enumerate(cleaner):
            if b in 'ACGT':
                #check if the qual value is high enough.
                if int(vector_of_quals[i]) > 1:
                    quality_alts.append(b)
                    depth += 1 #only count bases which have 2+ consensus representations as part of the depth
                    altbase_counts[b] = altbase_counts.get(b,0) + 1 #count the number of altbases per quality base.
        #altbase_counts[ref] = basedepth - sum(altbase_counts.values()) #the rest of basedepth are all reference and need to be accounted as such
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
                if 2 <= basecount <= basedepth/4: #doesn't make sense to calculate for singletons, which I intend to skip by default now.
                    key = (len(quality_alts), basecount)
                    if key not in pcr_duplicate_track:
                        pcr_duplicate_track[key] = perm_index(leng = key[0], num = key[1])
                    thresh = pcr_duplicate_track[key]
                    if dindeces[base] < thresh: #less than 5% chance of getting a cluster like this. Lock this one to 1 instance
                        #print("QC: Base is skipped for clustering")
                        skip += base
                        #still count it once though
                        pcrc.append(base)
                        altbase_counts[base] = 1 #set its count representation to 1. Note that all members are still counted separately for depth.
                    #if a base exists at appropriate levels but stays above the cluster threshold, it's collected in the fixed alts statement below
                else:
                    skip += base
        #collect individual bases not counted as a pcr cluster, plus a singleton of each cluster base
            #fixed_alts = ','.join(sorted(list(set(pcrc + [q for q in quality_alts if q not in skip])))) #save in order ACGT
            fixed_alts = ','.join(sorted(list(altbase_counts.keys()))) #save in order ACGT
        else: #retain higher frequencies, including PCR clusters which are indistinguishable from higher frequency mutations.
            fixed_alts = ','.join(sorted(list(set(quality_alts)))) #save in order ACGT
        aacountstr = ','.join([str(v/basedepth) for k,v in sorted(altbase_counts.items())])
        if depth == 0 or len(fixed_alts) == 0: #nothing but Ns or reference here.
            return None, pcr_duplicate_track
        else:
            if len(fixed_alts.split(',')) != len(altbase_counts):
                print("Strange Entry")
                print(spent)
                print(skip)
                print(fixed_alts, altbase_counts)
                return None, pcr_duplicate_track
            vcf_line = chrom + '\t' + loc + '\t.\t' + ref + '\t' + fixed_alts + '\t.\tPASS\tDP=' + str(basedepth) + ';AF=' + aacountstr
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
    alpha, beta = args.params
    for entry in pilein:
        nline, pdt = make_vcf_line(entry.strip().split(), args.mind, args.germline, pcr_duplicate_track=pdt)
        if nline != None:
            if args.bayes:
                nline = add_likelihood_info_bayes(nline, alpha, beta)
            else:
                nline = add_likelihood_info_mle(nline, alpha, beta)
            print(nline, file = outf)
    if args.pileup != None:
        pilein.close()
    if args.output != None:
        outf.close()

if __name__ == "__main__":
    main()
