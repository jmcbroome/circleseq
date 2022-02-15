#!/usr/bin/env python3

#script to remove consensus sequences which are identical to or encompassed by other consensus sequences
#this should conseratively remove all pcr duplicates and possibly some weaker data if it's shorter in length

import argparse
import numpy as np
import pysam

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', type = int, help = 'Set a minimum number of times a base must be seen. default 2', default = 2)
    parser.add_argument('-i', '--input', help = 'Path to input bam.', required = True)
    parser.add_argument('-o', '--output', help = 'Path to output bam. Default is input with ".dedup" added before the file extension.', default = None)
    parser.add_argument('-s', '--stats', help = 'Save statistics on number of alignments affected and total coverage removed to the target file.', default = None)
    parser.add_argument('-m', '--min_length', help = 'Set to a minimum value for the length of a consensus sequence to be retained in any case. Default 50', type = int, default = 50)
    args = parser.parse_args()
    return args

def write_stats(reads,outfname,mean_prior,mean_post):
    outf = open(outfname,'w+')
    print("Number of Reads Removed: " + str(len(reads)),file=outf)
    print("Mean Length of Reads Removed: " + str(np.mean([r.query_length for r in reads])),file=outf)
    print("Percentage of Reads Removed: " + str(mean_prior),file=outf)
    print("Percentage of Coverage Remaining: " + str(mean_post),file=outf)

def main():
    args = argparser()
    inputf = pysam.AlignmentFile(args.input,'rb')
    mapped_chromosomes = inputf.references
    bad_reads = set()
    total_covered = 0
    total_count = 0
    post_covered = 0
    for mc in mapped_chromosomes:
        allreads = inputf.fetch(mc)
        rsd = {}
        rlens = 0
        total_reads = 0
        for r in allreads:
            if r.query_length < args.min_length:
                bad_reads.add(r)
                continue
            if r.reference_start not in rsd:
                rsd[r.reference_start] = []
            rsd[r.reference_start].append(r)
            rlens += r.query_length
            total_reads += 1
        #pick a critical value- assume that the distribution of counts of number of reads which start at a location are normally distributed across the genome
        #and identify groups of reads that share a start and are unusual...
        cv = [len(v) for v in rsd.values()]
        critical = 2*np.std(cv) + np.mean(cv)
        for p,rc in rsd.items():
            if len(rc) > critical:
                #this site is considered a candidate location for pcr duplication
                #strip out the reads and any reads downstream they contain.
                #for now, don't even retain one copy of the putative duplicate.
                maxl = max([r.query_length for r in rc])
                for i in range(p,p+maxl-args.min_length): 
                    for r in rsd.get(i,[]):
                        if i + r.query_length <= p + maxl:
                            #mark this read for removal.
                            bad_reads.add(r)
        total_covered += rlens
        total_count += total_reads
    if args.output == None:
        outputname = args.input[:-4] + ".dedup.bam"
    else:
        outputname = args.output
    outputf = pysam.AlignmentFile(outputname,'wb',template=inputf)
    #finally, go over it again and write out every read that passed muster.
    total_removed = 0
    inputf.reset()
    for read in inputf.fetch():
        if read not in bad_reads:
            outputf.write(read)
        else:
            total_removed += read.query_length
    post_covered = total_covered - total_removed
    if args.stats != None:
        write_stats(bad_reads, args.stats, len(bad_reads)/total_count*100, post_covered/total_covered*100)    
    inputf.close()
    outputf.close()
if __name__ == "__main__":
    main()
