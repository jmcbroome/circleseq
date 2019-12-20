#!/usr/bin/env python3

#import
import argparse
import sys
#define functions/classes
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 2', default = 2)
    parser.add_argument('-e', '--errors', help = 'path to input file.')
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

def main():
    args = argparser()
    #insert code
    sumerrors = {}
    for a in "ACGT":
        for b in "ACGT":
            if a != b:
                sumerrors[(a,b)] = 0
    counts = {k:0 for k in 'ACGT'}
    with open(args.errors) as ef:
        for entry in ef:
            spent = entry.strip().split()
            if len(spent) > 4: #ignore empty lines from end of files etc
                ref = spent[2].upper() #don't particularly care whether its forward or reverse
                cir = ''.join(c for c in spent[4].upper() if c in 'ACGTN.') #remove characters indicating read start or end or mapping scores or whatever.
                if ref != 'N':
                    #depth = int(spent[3])
                    assert len(cir) == int(spent[3])
                    #counts[ref] += depth
                    if int(spent[3]) == len(spent[5]):
                        counts[ref] += len([v for v in cir if v != 'N']) #don't want to count bases the read has no information about as either reference or nonreference
                    #if depth == len(spent[5]): #ignore entries with indels or other complications for this iteration of the pipeline
                        for i, base in enumerate(cir): #may be length 1.
                            if base in "ACGT" and base != ref and spent[5][i] in '0123456789':
                                if int(spent[5][i]) >= args.threshold:
                                    sumerrors[(ref,base)] += 1 

    print("Error Counts")
    if counts == None:
        for k,v in sumerrors.items():
            print(k, '\t', v)
    else:
        for k,v in sumerrors.items():
            print(k, '\t', v, '\t', v/counts[k[0]]) #divide by the number of bases in the reference matching the original type

if __name__ == "__main__":
    main()
