#!/usr/bin/env python3

#build an artificial sam file of consensus mappings out of split.fa, the regions from the initial alignment, and the coordinates of those regions in the bed file.
#use N's at any point that's not very certain by whatever metric.
#IMPORTANT: SCRIPT ASSUMES PHRED + 33 QUALITY ENCODING. WILL BREAK WITH OTHER ENCODINGS.

import argparse

from pandas import qcut
import skbio.alignment as skaln
from multiprocessing import Pool
import re

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--consensus', help = "Set prefix of output 'sam' file for consensus alignments.")
    parser.add_argument('-t', '--threads', type = int, help = 'Set number of threads to use. Default 1', default = 1)
    parser.add_argument('bed', help = 'Path to a bed file containing coordinate information matching the primary fasta file.')
    parser.add_argument('regions', help = 'Path to a fasta file containing reference alignment regions from the original round of alignment')
    parser.add_argument('split', help = 'Path to file containing the split sequences (.fa)')
    args = parser.parse_args()
    return args

def rand_aln_p(bl, sl):
    """
    This simplified function yields the chance that a short sequence will randomly align to a longer sequence when both are composed of random letters. 
    Used to filter fragments that are too short to be reliably mappable over the reference region.
    """
    if sl < bl:
        return 1-(1-.25**sl)**(bl-sl)
    else:        
        return 1-(1-.25**bl)**(sl-bl)

def phred_code_qual(sym):
    """
    Convert between phred quality encodings (string character vs integer value).
    """
    if isinstance(sym, str):
        return ord(sym) - 33
    elif isinstance(sym, int):
        return chr(sym+33)
    else:
        print("Error: quality score conversion not understood. Check input")
        raise ValueError

def bedread(bed_path):
    """
    Read a named bed file (has a fourth "name" column) in as a simple dictionary of tuple locations corresponding to names. Assumes unique names.
    """
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
    """
    Read a fasta in as a simple dictionary. Used to process reference sequences for local realignment.
    """
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
    """
    Read a fastq in as a simple dictionary. Used to process fragments for local realignment.
    """
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
    """
    Compile a consensus sequence from local alignments.
    """
    region_counts, flag = basecount(region,seqs)
    #to create the consensus, simply iterate through region counts and indeces.
    consensus = []
    quality = [] #actually just count. makes stratifying easier.
    for i, bcs in sorted(region_counts.items(), key = lambda x:x[0]):
        #get the best base in bcs.
        b,q = best_base(bcs)
        consensus.append(b)
        quality.append(str(q))
    return [''.join(consensus), ''.join(quality)], flag

def best_base(basecounts):
    """
    Identify the base with the highest support for a given location. Used to compile consensus sequences.
    """
    #if all are 0, return an N.
    if all([c==0 for c in basecounts.values()]):
        return ('N',0)
    nonz = sorted([(b,c) for b,c in basecounts.items() if c > 0],key = lambda x:x[1], reverse = True)
    if len(nonz) == 1:
        return nonz[0]
    elif len(nonz) > 1:
        return ('N', nonz[0][1]) #retain depth but don't count as a real change if there's any controversy.
    else:
        return ('N',0) 

def basecount(region, seqs):
    """
    Perform a local smith-waterman realignment of sequence fragments to the targeted reference region and compile a consensus sequence from the resulting alignments. 4
    Removes fragments that can't be reliably mapped locally. 
    Ignores fragments which align with an indel and flags for quality control. Additionally flags poorly aligned consensi.
    Ignores bases which are below a quality score of 13 when compiling a consensus.
    """
    counts = {i:{b:0 for b in 'ACGT'} for i in range(len(region)+1)}
    flag = False #flag is a bool representing whether this particular sequence had an issue of any kind. Used to flag sam output.
    #print("ErrorChecker from basecount in make_brpileup: Number of sequences, lengths of sequences", len(seqs), [len(v[0]) for v in seqs])
    for seq, qual in seqs:
        if rand_aln_p(len(region), len(seq)) <= .0001: # effectively disabling filter for testing.
            seqaln = skaln.StripedSmithWaterman(seq, gap_open_penalty = 5, gap_extend_penalty = 2, match_score = 1, mismatch_score = -3, zero_index = False)
            aln = seqaln(region)
            #if indels in the cigar- e.g. not a simple parse- skip it.
            if aln['cigar'].count('I') == 0 and aln['cigar'].count('D') == 0:
                matched = 0
                for sec in aln['cigar'].split("M"):
                    reg = re.split('[A-Z]',sec)
                    if len(reg[-1]) > 0:
                        matched += int(reg[-1])
            else:
                flag = True
                #print('indel in query alignment {}, continuing'.format(aln['cigar']))
                continue
            #print("ErrorChecker from basecount in make_brpileup: Length of sequence, length aligned, length of target", len(seq),matched,len(region))
            index = aln['target_begin']
            if index != -1: #apparently this means no mapping?
                try:
                    rseq = aln.target_sequence[aln.target_begin:aln.target_end_optimal]
                except:
                    flag = True
                    continue #if I don't have a target sequence aligned value, I didn't get an alignment and I don't want to contiue.
                qseq = aln.query_sequence[aln.query_begin:aln.query_end]
                qqual = qual[aln.query_begin:aln.query_end]
                #print('\t'.join([str(v) for v in [len(aln.query_sequence),len(qseq),len(aln.target_sequence),len(rseq),aln.cigar,aln.optimal_alignment_score,aln.optimal_alignment_score/(len(rseq)+1)]]))
                if aln.optimal_alignment_score / (len(rseq)+1) < .5: #about .8% of the alignments are below this threshold in quality and shouldn't be used.
                    flag = True
                    continue
                for i, rb in enumerate(rseq):
                    assert len(qseq) == len(rseq)
                    try:
                        qb = qseq[i]
                        if qb in 'ACGT' and phred_code_qual(qqual[i]) > 12: #12 is the cutoff of <5% of being an error, e.g. qscore 13+ == <.05 of being error
                            counts[index + i + 1][qb] += 1
                    except:
                        print("Error trying to count alignment to region", region, aln.aligned_target_sequence, aln.query_sequence, aln['cigar'])#, index, i, rb, qb)
    return counts, flag

def trim_sam(seqqual, coord):
    """
    Remove no-information Ns and 0s from the alignment before saving output. Used as a final step for more parseable output.
    """
    #remove leading and trailing bases not actually mapped from the consensus alignment and update the starting location and total matched length appropriately.
    qual = seqqual[1]
    #some tricky code from stackoverflow... one line generator expressions!
    front = next((i for i, ch  in enumerate(qual) if ch != '0'),0)
    back = len(qual)-next((i for i, ch  in enumerate(qual[::-1]) if ch != '0'),0)
    new_seq = seqqual[0][front:back]
    new_qual = qual[front:back]
    assert len(new_seq) == len(new_qual)
    new_coords = [coord[0]]
    new_coords.append(str(int(coord[1]) + front))
    new_coords.append(str(int(coord[2]) - back))
    return tuple(new_coords), (new_seq,new_qual)

def sam_entry(seqqual, coord, name, flag = False, trim = True):
    """
    Construct a new sam entry from a compiled consensus sequence. Flag consensi which encountered low alignment quality problems in the local realignment step as 512 (failed QC).
    """
    coord = coord[0]
    if not flag:
        fv = '0'
    else:
        fv = '512' #ones with issues will be flagged as being below alignment quality
    #remove bases with quality 1 or less from the data if requested and update alignment coordinates appropriately.
    if trim:
        coord, seqqual = trim_sam(seqqual, coord)
    return name + "\t" + fv + "\t" + coord[0] + '\t' + str(coord[1]) + '\t60\t' + str(len(seqqual[0])) + "M\t*\t0\t0\t" + seqqual[0] + '\t' + seqqual[1]

def mapper(input_iter):
    """
    Coordinator function for consensus construction and output.
    """
    k, coord, reg, splts = input_iter
    cons, flag = generate_consensus(reg, splts) #version of the reference with either Ns for ambiguous/no information or the 
    samstr = sam_entry(cons, coord, k, flag)
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
            arguments = [(k, bed_d[k], reg_d[k], spl_d[k]) for k in spl_d.keys() if k in reg_d and k in bed_d]
            samstrs = p.imap_unordered(mapper, arguments)
            for s in samstrs:
                if s != '':
                    print(s, file = outf)

if __name__ == "__main__":
    main()
