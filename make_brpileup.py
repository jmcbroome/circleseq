#!/usr/bin/env python3

#build an artificial sam file of consensus mappings out of split.fa, the regions from the initial alignment, and the coordinates of those regions in the bed file.
#use N's at any point that's not very certain by whatever metric.
#IMPORTANT: SCRIPT ASSUMES PHRED + 33 QUALITY ENCODING. WILL BREAK WITH OTHER ENCODINGS.

#import
import argparse
import skbio.alignment as skaln
from multiprocessing import Pool
from functools import partial

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-c', '--consensus', help = "Set prefix of output 'sam' file for consensus alignments.")
    parser.add_argument('-t', '--threads', type = int, help = 'Set number of threads to use. Default 1', default = 1)
    parser.add_argument('bed', help = 'Path to a bed file containing coordinate information matching the primary fasta file.')
    parser.add_argument('regions', help = 'Path to a fasta file containing reference alignment regions from the original round of alignment')
    parser.add_argument('split', help = 'Path to file containing the split sequences (.fa)')
    args = parser.parse_args()
    return args

def rand_aln_p(bl, sl):
    #use a probability function that yields the chance that SL will randomly align perfectly to BL as a filter for too-short fragments. 
    # assert sl < bl
    if sl < bl:
        return 1-(1-.25**sl)**(bl-sl)
    else:        
        return 1-(1-.25**bl)**(sl-bl)

def phred_code_qual(sym):
    if isinstance(sym, str):
        return ord(sym) + 33
    elif isinstance(sym, int):
        return chr(sym-33)
    else:
        print("Error: quality score conversion not understood. Check input")
        raise ValueError

def bedread(bed_path):
    bd = {}
    with open(bed_path) as bedf:
        for entry in bedf:
            spent = entry.strip().split()
            read = spent[3]
            reg = tuple(spent[0:3])
            bd[read] = bd.get(read, [])
            if reg not in bd[read]:
                bd[read].append(reg)
    return bd

def regionread(region_path):
    rd = {}
    cur = None
    with open(region_path) as regf:
        for entry in regf:
            if entry[0] == '>':
                cur = entry.strip()[1:].partition('::')[0] #forwards compatibility.
            else:
                rd[cur] = entry.strip()
    return rd

def splitread(split_path):
    sd = {}
    cur = None
    with open(split_path) as splf:

        for entry in splf:
            if entry[0] == '+':
                continue #useless lines.
            if entry[0] == '@':
                cur = entry.partition('_')[0][1:]
                sd[cur] = sd.get(cur, []) #initialize entry.
                tmp = []
            elif len(tmp) == 0: #add the actual sequence.
                tmp.append(entry.strip())
            elif len(tmp) == 1:
                tmp.append(entry.strip())
                sd[cur].append(tmp)
    return sd

def generate_consensus(region, seqs):
    region_counts = basecount(region,seqs)
    #to create the consensus, simply iterate through region counts and indeces.
    consensus = []
    quality = [] #actually just count. makes stratifying easier.
    for i, bcs in sorted(region_counts.items(), key = lambda x:x[0]):
        #get the best base in bcs.
        b,q = best_base(bcs)
        consensus.append(b)
        quality.append(str(q))
    return ''.join(consensus), ''.join(quality)

def best_base(basecounts):
    #if all are 0, return an N.
    if all([c==0 for c in basecounts.values()]):
        return ('N',0)
    nonz = [(b,c) for b,c in basecounts.items() if c > 0]
    if len(nonz) == 1:
        return nonz[0]
    else:
        return ('N',0) #for now, ignoring ANY level of ambiguity- outside of bad bases, ofc. This section can be easily altered in the future.

def basecount(region, seqs):
    counts = {i:{b:0 for b in 'ACGT'} for i in range(len(region)+1)}
    #build the align object for the region
    ssw = skaln.StripedSmithWaterman(region, gap_open_penalty = 5, gap_extend_penalty = 1, match_score = 1, mismatch_score = -2, zero_index = False)
    for seq, qual in seqs: 
        if rand_aln_p(len(region), len(seq)) < .0001:
            aln = ssw(seq)
            #convert results to look like pairwise2 output.
            #if aln['optimal_alignment_score'] <= len(seq) - 4:
                #print('Skipping for low score')
             #   continue
            #if indels in the cigar- e.g. not a simple parse- skip it.
            try:
                matched = int(aln['cigar'][:-1])
            except:
                #print('indel in alignment {}, continuing'.format(aln['cigar']))
                continue
            if matched < len(seq) - 4:
                index = int(aln['query_begin'])
                for tar_ind in range(aln['target_begin'], aln['target_end_optimal']):
                    if index + tar_ind in counts:
                        base = aln['target_sequence'][tar_ind]
                        if base in 'ACGT':
                            counts[index + tar_ind][base] += 1
            # aln = pairwise2.align.localms(region, seq, 1, -2, -5,-1, one_alignment_only = True)
            #set a minimum allowable alignment score to count this. If no more than 1 mismatch, say, minimum score should be length of alignment -1 .
            # if len(aln) == 0 or aln[0][2] <= len(seq) - 4: #score for a sequence with a single snp mismatch is length-3
                # continue
            # else:
                # aln = aln[0]
            # if len(aln[1].strip('-')) >= len(seq) - 4: #maximum allowed to be trimmed off by the local, say. if it barely aligns any of it, don't trust it at all.
                # for i, base in enumerate(aln[1]): #the small sequence
                    # if base in 'ACGT'): #skip the nonmap '-', obvs.
                        # try:
                            # counts[i][base] += 1 #think on how to incorporate quality scores here.
                        # except KeyError:
                            #alignment probably just ran off the edge. ignore it.
                            # continue
                            # print(len(region), region)
                            # print(len(seq), seq)
                            # print(aln)
                            # print(i)
                            # raise AssertionError
    return counts

def sam_entry(seqqual, coord, name):
    #build a fake sam consensus entry.
    # print(name)
    # print(coord[0])
    coord = coord[0]
    # print(seqqual)
    return name + "\t0\t" + coord[0] + '\t' + str(coord[1]) + '\t60\t' + str(len(seqqual[0])) + "M\t*\t0\t0\t" + seqqual[0] + '\t' + seqqual[1]

def mapper(input_iter):
    # k, bed_d, reg_d, spl_d = input_iter

    # try:
        # coord = bed_d[k]
    # except KeyError:
        # print("Error: read {} does not have matching bed information. Continuing".format(k))
        # return '' #this read doesn't have matching info- possibly because it failed to map originally.
    # try:
        # reg = reg_d[k]
    # except KeyError:
        # print("Error: read {} does not have matching reference region information. Continuing".format(k))
        # return ''
    # splts = spl_d[k]
    k, coord, reg, splts = input_iter
    cons = generate_consensus(reg, splts) #version of the reference with either Ns for ambiguous/no information or the 
    samstr = sam_entry(cons, coord, k)
    return samstr

def main():
    #first step is to read all files into similar dictionary objects
    #each object key is the base read name
    #values vary by dictionary
    args = argparser()
    bed_d = bedread(args.bed) #values are [(chro, start, stop)] for bed regions assigned to that name. May be multiple regions.
    reg_d = regionread(args.regions) #values are straight up just a sequence of ACTG, to act as a base for the local alignment.
    spl_d = splitread(args.split) #values are all split up short sequences assigned to that read name and their quality scores
    with open(args.consensus, 'w+') as outf:
        with Pool(args.threads) as p:
            # mapper_d = partial(mapper, bed_d = bed_d, reg_d = reg_d, spl_d = spl_d)
            arguments = [(k, bed_d[k], reg_d[k], spl_d[k]) for k in spl_d.keys() if k in reg_d and k in bed_d]
            samstrs = p.imap_unordered(mapper, arguments)
            for s in samstrs:
                if s != '':
                    print(s, file = outf)

if __name__ == "__main__":
    main()
