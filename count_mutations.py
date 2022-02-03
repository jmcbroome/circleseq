#!/usr/bin/env python3

import argparse
import numpy as np
import statistics as st

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 2', default = 2)
    parser.add_argument('-m', '--mutations', help = 'path to input file.')
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
    sumerrors = {}
    for a in "ACGT":
        for b in "ACGT":
            if a != b:
                sumerrors[(a,b)] = 0
    counts = {k:0 for k in 'ACGT'}
    pcr_duplicate_track = {} #using dynamic programming to save compute cycles for this qc measure
    with open(args.mutations) as ef:
        for entry in ef:
            spent = entry.strip().split()
            ref = spent[2].upper()
            if len(spent) > 4 and ref != "N": #ignore empty lines from end of files etc
                spent = entry.strip().split()
                alts = [b for b in spent[4] if b in 'ACGTN.']
                quals = spent[5]
                try:
                    assert len(alts) == len(quals)
                except:
                    #print("Alt and Qual Mismatch:", entry)
                    continue
                nalts = ''
                nquals = ''
                for i, base in enumerate(alts):
                    if int(quals[i]) >= args.threshold and base != 'N': #strip out Ns and low quality alleles
                        nalts += base
                        nquals += quals[i]
                counts[ref] += len(nalts) #number of bases I have sufficient information about (not Ns and not 1x circles)
                #apply filters for calling mutations here.
                #first, the depth must be at least ten in order to differentiate between germline and somatic mutations.
                #depth being the non-N content of the alternative allele string.
                #second, any mutations which exist at higher than a 25% frequency in the string are probably germline and should be ignored for somatic mutation analysis.
                if len(nalts) > 12: #try setting this to 12 for now to best distinguish?
                    # [nalts.count(base) < len(nalts)/4 for base in 'ACGT']):
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
                                sumerrors[(ref,base)] += 1 #still record it once.
                        else:
                            #print("QC: Base is singleton or too high frequency in pileup")
                            skip += base
                    #now actually go through and count the scattered mutations.
                    for base in nalts:
                        if base != '.' and base not in skip:
                            sumerrors[(ref,base)] += 1
                else:
                    # print("QC: Read is skipped for having no alts")
                    continue
    print("Error Counts")
    if counts == None:
        for k,v in sumerrors.items():
            print(k, '\t', v)
    else:
        for k,v in sumerrors.items():
            print(k, '\t', v, '\t', v/counts[k[0]]) #divide by the number of bases in the reference matching the original type

if __name__ == "__main__":
    main()
