#!/usr/bin/env python3

#script to remove Ns, bases without enough support from the original consensus reads, and PCR duplicates from the raw mpileup of consensus reads.
#additionally filters on depth to return only regions where somatic and germline mutations can be distinguished and pcr duplicates identified

#import
import argparse
import sys
import numpy as np
import statistics as st

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 2', default = 2)
    parser.add_argument('-i', '--input', help = 'path to input file.', default = None)
    parser.add_argument('-o', '--output', help = 'name of output pileup, default is stdout', default = None)
    parser.add_argument('-s', '--stats', help = 'save statistics on number of reads affected and total coverage removed')
    args = parser.parse_args()
    return args

def count_bases(path):
    counts = {k:0 for k in ['A','C','G','T']}
    with open(path) as inf:
        for line in inf:
            if line[0] != '>':
                for base in line.strip():
                    if base in counts:
                        counts[base] += 1
    return counts

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

def main():
    args = argparser()
    good_entries = []

    pcr_duplicate_track = {} #using dynamic programming to save compute cycles for this qc measure
    if args.input == None:
        inputf = sys.stdin
    else:
        inputf = open(args.input)
    for entry in inputf:
        spent = entry.strip().split()
        ref = spent[2].upper()
        if len(spent) > 4 and ref != "N": #ignore empty lines from end of files etc
            spent = entry.strip().split()
            alts = [b for b in spent[4] if b in 'ACGTN.']
            quals = spent[5]
            assert len(alts) == len(quals)
            nalts = ''
            nquals = ''
            for i, base in enumerate(alts):
                if int(quals[i]) >= args.threshold and base != 'N': #strip out Ns and low quality alleles
                    nalts += base
                    nquals += quals[i]
            #apply filters for calling mutations here.
            #first, the depth must be at least five in order to differentiate between germline and somatic mutations.
            #depth being the non-N content of the alternative allele string.
            #second, any mutations which exist at higher than a 25% frequency in the string are probably germline and should be ignored for somatic mutation analysis.
            if len(nalts) > 5 and all([nalts.count(b) < len(nalts)/4 for b in 'ACGT']):
                #now, apply the pcr duplicate permutation filter structure using functions above.
                skip = '' #record no more than one of the bases that will be included here because of pcr duplicate inflation.
                dindeces = get_dindex(nalts)
                for base in 'ACGT':
                    basecount = nalts.count(base)
                    if 2 <= basecount <= len(nalts)/4: #doesn't make sense to calculate for singletons, which I intend to skip by default now.
                        key = (len(nalts), nalts.count(base))
                        if key not in pcr_duplicate_track:
                            pcr_duplicate_track[key] = perm_index(leng = key[0], num = key[1])
                        thresh = pcr_duplicate_track[key]
                        if dindeces[base] < thresh: #less than 5% chance of getting a cluster like this. Lock this one to 1 instance
                            #print("QC: Base is skipped for clustering")
                            skip += base
                    else:
                        #print("QC: Base is singleton or too high frequency in pileup")
                        skip += base
                #now actually go through and count the scattered mutations.
                recorded = []
                dnalts = ''
                dnquals = ''
                for i,base in enumerate(nalts):
                    if base in skip and base in recorded:
                        continue
                    recorded.append(base) #it can only go in once if it's in skip from the pcr duplicate remover.
                    dnalts += base
                    dnquals += str(nquals[i])
                #reconstruct the entry and append it to the output.
                nent = spent
                nent[4] = dnalts
                nent[5] = dnquals
                nent[3] = len(dnalts)
                good_entries.append('\t'.join([str(v) for v in nent]))
            else:
                # print("QC: Read is skipped for having no alts")
                continue
    if args.output == None:
        outf = sys.stdout
    else:
        outf = open(args.output, 'w+')
    for nen in good_entries:
        print(nen, file = outf)   
    if args.input != None:
        inputf.close()
    elif args.output != None:
        outf.close()
if __name__ == "__main__":
    main()
