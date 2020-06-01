#!/usr/bin/python3
#the concept behind this script is that germline mutations will be shared between separate samples from the same strain, but somatic mutations will not
#this script takes any number of VCF files, reads each, and records how many times each entry is seen as encoded by chromosome - location - alternate, regardless of frequency
#it then returns two separate files; one will contain unique mutations, and one will contain germline (seen in multiple samples) mutations
import argparse
import sys
import numpy as np

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--somatic', help = 'Filename for somatic mutation output', default = 'somatic.vcf')
    parser.add_argument('-g', '--germline', help = 'Filename for germline mutation output', default = 'germline.vcf')
    parser.add_argument('-d', '--diagnostic', help = 'Set to a filename to save the count dictionary for mutations seen to a file. Default is to not save', default = None)
    parser.add_argument('-m', '--mode', help = 'Set to VCF or TXT to parse VCF or plaintext files. Default is VCF', default = 'VCF')
    parser.add_argument('inputs', nargs = '+', help = 'Paths to any number of input vcf files')
    args = parser.parse_args()
    return args

def iterator(paths, somatic = 'somatic.vcf', germline = 'germline.vcf', mode = 'VCF', diagnostic = None):
    seen = {}
    for vcf in paths:
        if mode == 'VCF':
            header, body = parse_vcf(vcf)
        elif mode == 'TXT':
            header, body = parse_text(vcf)
        else:
           print("Mode invalid")
           raise ValueError
        for key, entry in body.items():
            if key not in seen:
                seen[key] = (1, entry, [vcf])
            else:
                seen[key] = (seen[key][0]+1, entry, seen[key][2] + [vcf])
    if diagnostic != None:
        with open(diagnostic, 'w+') as outf:
            for key, values in seen.items():
                count, entry, vcflist = values
                print(entry + '\t' + str(count) + '\t' + ','.join(vcflist), file = outf)

    with open(somatic, 'w+') as sout:
        printed = set()
        if mode == 'VCF':
            for line in header: #print the header (of the last parsed vcf, but they should all match anyways.)
                print(line, file = sout)
        #note that the output is unsorted, and may have to be passed through bcftools sort or other sorting tool options
        for key, cent in {k:v for k,v in seen.items() if v[0] < 2}.items():
            if cent[1] not in printed:
                print(cent[1], file = sout)
                printed.add(cent[1])
    with open(germline, 'w+') as gout:
        printed = set()
        if mode == 'VCF':
            for line in header:
                print(line, file = gout)
        for key, cent in {k:v for k,v in seen.items() if v[0] >= 2}.items():
            if cent[1] not in printed:
                print(cent[1], file = gout)
                printed.add(cent[1])

def parse_text(infile):
    header = None
    body = {}
    with open(infile) as inf:
        for entry in inf:
            spent = entry.strip().split()
            body[(spent[0],spent[1],spent[3])] = entry.strip()
    return header, body

def parse_vcf(infile):
    header = []
    body = {}
    with open(infile) as inf:
        for entry in inf:
            if entry[0] == '#':
                header.append(entry.strip())
            else:
                spent = entry.strip().split()
                for m in spent[4].split(','):
                    body[(spent[0],spent[1],m)] = entry.strip()
    return header, body

def main():
    args = argparser()
    iterator(args.inputs, args.somatic, args.germline, args.mode, args.diagnostic)

if __name__ == '__main__':
    main()
