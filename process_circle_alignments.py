#!/usr/bin/env python3

import argparse
import re
from Bio import pairwise2

def create_fasta(entry, name):
    '''
    Process a sam alignment back into N fasta entries by splitting based on cigar string.
    Additionally, if more bases are trimmed than belong to the primary alignment match, an additional split will be performed on the long trimmed section of the length of the primary alignment under the circle assumption.
    '''
    cigar = re.findall(r'(\d+)(\w)', entry[5])
    indeces = [0] + [int(i[0]) for i in cigar] + [-1]
    nseqs, nquals = [[entry[v][indeces[i]:indeces[i] + indeces[i+1]]
        for i in range(len(indeces)-1)] 
        for v in [9,10]]
    i = 0
    entries = []
    for seq, qual in zip(nseqs, nquals):
        if len(seq) > 0:
            subname = name + '_' + str(i)
            nent = build_fa_entry(subname, seq, qual)
            i += 1
            entries.append(nent)
    return entries    

def build_fa_entry(name, seq, qual):
    assert len(seq) == len(qual)
    return '@' + name.strip() + '\n' + seq.strip() + '\n+\n' + qual.strip()

def additional_split(entry):
    '''
    Uses a local alignment of two subsections of sequence to generate a new cigar string and returns the altered alignment with cigar.
    Returns the entry list if it doesn't fill conditions.
    '''
    #a lot of this calculation is technically redundant but I want to be able to disable this altogether if I want to skip this step
    #also the vast majority of reads won't be entering this function so overall it doesn't add much redundant computation
    cigar = re.findall(r'(\d+)(\w)', entry[5])
    if len(cigar) == 1:
        return entry #whole thing mapped, probably from a repetitive region, just throw it back
    #find the largest span in the cigar, if its an S do another splitting step based on the results of find_align
    big = max(enumerate(cigar), key = lambda x:int(x[1][0]))
    if big[1][1] == 'S': #split this big one further.
        #grab the biggest matching region. which should also be the only one in this case.
        match = max([(i,c) for i,c in enumerate(cigar) if c[1] == 'M'], key = lambda x:int(x[1][0]))
        assert int(match[1][0]) <= int(big[1][0]) #should definitely be true if we got here.
        #now split up the sequence based on the cigar.
        indeces = [0] + [int(i[0]) for i in cigar] + [-1]
        nseqs = [entry[9][indeces[i]:indeces[i] + indeces[i+1]]
            for i in range(len(indeces)-1)] 
        #feed the appropriate subsequences to find_align
        #remove empties from nseqs 
        nseqs = [n for n in nseqs if n != '']
        nbs = find_align(nseqs[big[0]], nseqs[match[0]])
        #build a new cigar entry over the big entry and pop the original.
        #split the cigar over the index and then knit it all back together
        #using S as the marker because it's processed the same way once it gets passed back to the primary parser
        if nbs != None:
            ncig = cigar[:big[0]] + [(k, 'S') for k in nbs] + cigar[big[0]+1:]
            #reknit the cigar, overwrite entry, and return.
            entry[5] = ''.join([''.join([str(v) for v in k]) for k in ncig])
        return entry
    else:
        return entry #it's good as is.

def find_align(bigseq, smallseq):
    '''
    Perform an alignment using Bio.pairwise2 of two sequence strings to find the best matching subsequence 
    returns the lengths of matching sequence and nonmatching sequences at either end.
    Currently does not filter based on alignment quality but doesn't allow gapping.
    '''
    assert len(bigseq) >= len(smallseq)
    aln = pairwise2.align.localms(bigseq, smallseq, 1, -2, -5,-1, one_alignment_only = True)
    if len(aln) == 0:
        return None
    else:
        aln = aln[0]
    #collect the start and stop of alignment of the alignment of the second string to the first
    #three new sections: trim up to where trim matches match, match, trim after matching match ends
    return [aln[-2], aln[-1]-aln[-2], len(bigseq)-aln[2]]

def create_bed(aln, d, cstops):
    """
    Create a bed entry based on alignment plus distance. Used downstream for retrieving areas of genome for local realignment.
    """
    #bed format is chromosome (or in this case region or scaffold or whatever.), name, start, stop
    #rind is the .fai file that tells me where my boundaries are for defining the bed
    name = aln[0]
    chro = aln[2]
    start = max(int(aln[3]) - d, 1)
    stop = min(int(aln[3]) + len(aln[9]) + d, cstops[chro])
    bstr = '\t'.join([chro, str(start), str(stop), name])
    return bstr

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference_index', help = 'Path to a .fai produced by samtools index for checking chromosome lengths when defining features.')
    parser.add_argument('-f', '--fasta_pref', help = 'Choose a name for the fasta file. Default is split.fa', default = 'split.fa')
    parser.add_argument('-b', '--bed_pref', help = 'Choose a name for the bed file. Default is regions.bed', default = 'regions.bed')
    parser.add_argument('-d', '--distance', type = int, help = 'Set to a value of distance around alignments to use for local area remapping of secondary alignments downstream.', default = 10)
    parser.add_argument('sam', help = 'Path to files containing the primary single alignments for all reads')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    seen = set()
    cstops = {}
    with open(args.reference_index) as inf:
        for entry in inf:
            spent = entry.strip().split()
            cstops[spent[0]] = int(spent[1])
    with open(args.sam) as samf:
        with open(args.fasta_pref, 'w+') as outf:
            with open(args.bed_pref, 'w+') as outb:
                for entry in samf:
                    if entry[0] != '@': #skip headers
                        spent = entry.strip().split()
                        #remove entries with flags > 1000, unaligned reads *, hard clipped, deletions and insertions.
                        badcig = any([k != -1 for k in [spent[5].find(let) for let in 'HID']])
                        if int(spent[1]) > 1000 or spent[5] == '*' or badcig or int(spent[4]) == 0:
                            continue
                        if spent[0] in seen:
                            name = spent[0] + '_1' #its the second read primary alignment
                        else:
                            name = spent[0] + '_0'
                            seen.add(spent[0])
                            #adding additional functionality: now creating a bed file based on the read index areas.
                            bent = create_bed(spent, args.distance, cstops)
                            print(bent, file = outb)
                        # for l in create_fasta(additional_split(spent), name):
                        for l in create_fasta(spent, name):
                            print(l, file = outf)

if __name__ == "__main__":
    main()
