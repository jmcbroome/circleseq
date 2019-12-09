#!/usr/bin/env python3

#import
import argparse
import sys

#define functions/classes
def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    #add args
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 1', default = 1)
    parser.add_argument('-e', '--errors', help = 'path to input file. default is stdin', default = sys.stdin)
    parser.add_argument('-s', '--seqs', help = 'path to file containing all covered reference sequences for rate calculations. Exclude to print only hard counts', default = None)
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
    if args.seqs != None:
        counts = count_bases(args.seqs)
    else:
        counts = None
    with open(args.errors) as ef:
        for entry in ef:
            spent = entry.strip().split()
            if len(spent) > 4:
                ref = spent[2].upper()
                cir = spent[4].upper()
                if cir in 'ACGT' and ref in "ACGT" and spent[5] != '/' and spent[5] != '.' and len(cir) == 1:
                    if cir != ref and int(spent[5]) >= args.threshold:
                        # if ord(spent[5]) + 33 > 80:
                        sumerrors[(ref, cir)] += 1
    print("Error Counts")
    if counts == None:
        for k,v in sumerrors.items():
            print(k, '\t', v)
    else:
        for k,v in sumerrors.items():
            print(k, '\t', v, '\t', v/counts[k[0]]) #divide by the number of bases in the reference matching the original type

if __name__ == "__main__":
    main()
