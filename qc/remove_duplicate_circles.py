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
    args = parser.parse_args()
    return args

def is_encompassed(r1,r2):
    #we consider r2 to be encompassed by r1 if they have exactly the same length or r2 is shorter
    #and r2 ends at the same point or earlier
    #this is a conservative filtering step.
    return (r1.reference_start <= r2.reference_start and r2.query_length <= r1.query_length and r1.reference_start + r1.query_length >= r2.reference_start + r2.query_length)

def write_stats(reads,outfname,mean_prior,mean_post):
    outf = open(outfname,'w+')
    print("Number of Reads Removed: " + str(len(reads)),file=outf)
    print("Mean Length of Reads Removed: " + str(np.mean([r.query_length for r in reads])),file=outf)
    print("Mean Coverage Initial: " + str(mean_prior),file=outf)
    print("Mean Coverage After Removal: " + str(mean_post),file=outf)

def main():
    args = argparser()
    inputf = pysam.AlignmentFile(args.input,'rb')
    mapped_chromosomes = inputf.references
    bad_reads = set()
    total_covered = 0
    total_count = 0
    post_covered = 0
    for mc in mapped_chromosomes:
        #get the set of reads aligned to each base in the chromosome
        #sort by size of alignment and ask whether each is encompassed by the larger reads in their set
        #if they are, tag them and ignore them in the future and don't output them
        #if they aren't, they can be saved to the output
        for pcol in inputf.pileup(mc):
            #previously marked reads can also be ignored as candidates for encompasser, since their encompasser will always encompass anything they would encompass
            #print(dir(pcol.pileups[0].alignment))
            allreads = pcol.pileups
            total_covered += len(allreads)
            total_count += 1
            reads = sorted([p.alignment for p in allreads if p.alignment not in bad_reads], key=lambda x:x.query_length)
            icheck = 1
            removed_here = 0
            for r in reads:
                for br in reads[icheck:]:
                    if is_encompassed(br,r):
                        bad_reads.add(r)
                        removed_here += 1
                        break

                icheck += 1
            post_covered += len(reads) - removed_here
    if args.output == None:
        outputname = args.input[:-4] + ".dedup.bam"
    else:
        outputname = args.output
    outputf = pysam.AlignmentFile(outputname,'wb',template=inputf)
    #finally, go over it again and write out every read that passed muster.
    inputf.reset()
    for read in inputf.fetch():
        if read not in bad_reads:
            outputf.write(read)

    if args.stats != None:
        write_stats(bad_reads, args.stats, total_covered/total_count, post_covered/total_count)    
    inputf.close()
    outputf.close()
if __name__ == "__main__":
    main()
