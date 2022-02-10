#script that takes a set of paired reads and simulates PCR duplication of the original reads by generating additional copies and rotating their content
#it saves the expanded fasta files and a logging file containing the names of duplicated reads and how much rotation was selected

import argparse
import numpy as np
from scipy.stats import poisson
import random
import gzip
import sys

class Read:
    def __init__(self,n,seq,qual):
        self.name = n
        self.seq = seq
        self.qual = qual
    def rotate(self):
        """
        Randomly chooses a split point, and appends everything before the split point to the string after the split point, then returns the rotated string.
        This is intended to simulate taking a different section of the identical sequence of repeats resulting from a duplicated circle.
        It's possible to return an identical sequence and quality string.
        """
        split = random.choice(range(len(self.seq)))
        if split == 0 or split == len(self.seq):
            return (self.seq,self.qual)
        ns = self.seq[split:] + self.seq[:split]
        nq = self.qual[split:] + self.qual[:split]
        return (ns,nq)
    def write(self,outf=sys.stdout):
        print(self.name,file=outf)
        print(self.seq,file=outf)
        print("+",file=outf)
        print(self.qual,file=outf)

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help = 'prefixes of input fastqs. Assumes fastq files end in "_R*.fq.gz".', required = True)
    parser.add_argument('-o', '--output', help = 'name of output fastqs. Default is the same as input with "_duplicated" inserted before the file extension', default = None)
    parser.add_argument('-s', '--stats', help = 'Record number of reads duplicated.',default='duplicate.log')
    parser.add_argument('-r', '--rate', help='Choose a proportion of reads to undergo duplication. Bounded 0-1. Default 0.1',default=0.1,type=float)
    parser.add_argument('-m', '--mu', help='Choose an amount of duplication to perform. Based on a poisson distribution. Default mu is 10.',type=float,default=10)
    args = parser.parse_args()
    return args

def generate_duplicates(r1,r2,mu):
    count = poisson.rvs(mu)
    nreads = []
    for c in range(count):
        n1,q1 = r1.rotate()
        n2,q2 = r2.rotate()
        nr1 = Read(r1.name + "_" + str(c), n1, q1)
        nr2 = Read(r2.name + "_" + str(c), n2, q2)
        nreads.append((nr1,nr2))    
    return nreads, count

def generate_reads(allrs, rate, mu):
    all_newreads = []
    total_duplicated = 0
    for r1,r2 in allrs:
        if random.random() < rate:
            nrv, c = generate_duplicates(r1,r2,mu)
            all_newreads.extend(nrv)
            total_duplicated += c
    return all_newreads, total_duplicated

def read_fasta(f):
    fd = {}
    cr = ""
    cseq = ""
    cquals = ""
    adding_seq = True
    with gzip.open(f,'rt') as inf:
        for entry in inf:
            if entry[0] == "@" and len(cquals) == len(cseq):
                if cr != "":
                    fd[cr] = (cseq,cquals)
                cseq = ""
                cquals = ""
                cr = entry.strip()
                adding_seq = True
            elif entry[0] == "+":
                adding_seq = False
            elif adding_seq:
                cseq += entry.strip()
            else:
                cquals += entry.strip()
    return fd

def gather_reads(f1,f2):
    readpairs = []
    unpaired_count = 0
    R1 = read_fasta(f1)
    R2 = read_fasta(f2)
    all_names = set(list(R1.keys()) + list(R2.keys()))
    for n in all_names:
        #only save paired reads.
        if n in R1 and n in R2:
            cr1 = Read(n,R1[n][0],R1[n][1])
            cr2 = Read(n,R2[n][0],R2[n][1])
            readpairs.append((cr1,cr2))
        else:
            unpaired_count += 1
    return readpairs, unpaired_count

def write_fastqs(reads,outpref):
    with gzip.open(outpref + "_R1.fq.gz","wt+") as outf1:
        with gzip.open(outpref + "_R2.fq.gz","wt+") as outf2:
            for r1,r2 in reads:
                r1.write(outf1)
                r2.write(outf2)

def main():
    args = argparser()
    if args.rate < 0 or args.rate > 1:
        print("ERROR: Choose a rate value between 0 and 1.")
        exit(1)
    log = open(args.stats,'w+')
    reads, upc= gather_reads(args.input + "_R1.fq.gz", args.input + "_R2.fq.gz")
    print("Total True Reads: " + str(len(reads)),file=log)
    print("Unpaired Reads Removed: " + str(upc),file=log)
    new_reads,total_created = generate_reads(reads, args.rate, args.mu)
    print("Duplicates Generated: " + str(total_created),file=log)
    reads.extend(new_reads)
    if args.output == None:
        outpref = args.input + "_duplicated"
    else:
        outpref = args.output
    write_fastqs(reads,outpref)

if __name__ == "__main__":
    main()