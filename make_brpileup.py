#!/usr/bin/env python3

#build an artificial sam file of consensus mappings out of split.fa, the regions from the initial alignment, and the coordinates of those regions in the bed file.
#use N's at any point that's not very certain by whatever metric.
#IMPORTANT: SCRIPT ASSUMES PHRED + 33 QUALITY ENCODING. WILL BREAK WITH OTHER ENCODINGS.

#import
import argparse
import skbio.alignment as skaln
from multiprocessing import Pool
from functools import partial
import re
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
        return ord(sym) - 33
    elif isinstance(sym, int):
        return chr(sym+33)
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

def generate_consensus(region, seqs, id = None):
    region_counts = basecount(region,seqs, id)
    #to create the consensus, simply iterate through region counts and indeces.
    consensus = []
    quality = [] #actually just count. makes stratifying easier.
    for i, bcs in sorted(region_counts.items(), key = lambda x:x[0]):
        #get the best base in bcs.
        b,q = best_base(bcs)
        consensus.append(b)
        quality.append(str(q))
    #try an additional step removing useless ns from the end to make certain stuff easier. just the end because don't want to mess with trying to update the index of the start.
    #stripped = ''.join(consensus).rstrip("N")
    #stripped_quality = quality[:len(stripped)]
    #assert len(stripped) = len(stripped_quality)
    return ''.join(consensus), ''.join(quality)
    #return stripped, stripped_quality

def best_base(basecounts):
    #if all are 0, return an N.
    if all([c==0 for c in basecounts.values()]):
        return ('N',0)
    #nonz = [(b,c) for b,c in basecounts.items() if c > 0]
    nonz = sorted([(b,c) for b,c in basecounts.items() if c > 0],key = lambda x:x[1], reverse = True)
    if len(nonz) == 1:
        return nonz[0]
    elif len(nonz) > 1:
        return ('N', nonz[0][1]) #retain depth but don't count as a real change if there's any controversy.
    else:
        return ('N',0) 

def basecount(region, seqs, id = None):
    counts = {i:{b:0 for b in 'ACGT'} for i in range(len(region)+1)}
    #build the align object for the region
    #ssw = skaln.StripedSmithWaterman(region, gap_open_penalty = 5, gap_extend_penalty = 1, match_score = 1, mismatch_score = -1, zero_index = False)
    #print("ErrorChecker from basecount in make_brpileup: Number of sequences, lengths of sequences", len(seqs), [len(v[0]) for v in seqs])
    for seq, qual in seqs:
        #if True:
        if rand_aln_p(len(region), len(seq)) <= .0001: # effectively disabling filter for testing.
            #print('len of seq', len(seq))
            #aln = ssw(seq)
            seqaln = skaln.StripedSmithWaterman(seq, gap_open_penalty = 5, gap_extend_penalty = 2, match_score = 1, mismatch_score = -3, zero_index = False)
            aln = seqaln(region)
            #print('prop of seq aligned', (aln['target_end_optimal']-aln['target_begin'])/len(seq))
            #convert results to look like pairwise2 output.
            #if aln['optimal_alignment_score'] <= len(seq) - 4:
                #print('Skipping for low score')
             #   continue
            #if indels in the cigar- e.g. not a simple parse- skip it.
            if aln['cigar'].count('I') == 0 and aln['cigar'].count('D') == 0:
                #matched = int(aln['cigar'][:-1])
                matched = 0
                for sec in aln['cigar'].split("M"):
                    reg = re.split('[A-Z]',sec)
                    if len(reg[-1]) > 0:
                        matched += int(reg[-1])
            else:
                #print('indel in query alignment {}, continuing'.format(aln['cigar']))
                continue
            #print("ErrorChecker from basecount in make_brpileup: Length of sequence, length aligned, length of target", len(seq),matched,len(region))
            #if matched >= len(seq) - 4: # - 4: #should always be true atm.
            if True: 
                index = aln['target_begin']
                #for tar_ind in range(aln['target_begin'], aln['target_end_optimal']):
                if index != -1: #apparently this means no mapping?
                    try:
                        rseq = aln.target_sequence[aln.target_begin:aln.target_end_optimal]
                    except:
                        continue #if I don't have a target sequence aligned value, I didn't get an alignment and I don't want to contiue.

                        #rseq = aln['target_sequence'][aln['target_begin']:aln['target_end_optimal']]
                    #try:
                        #qseq = aln.aligned_query_sequence
                    #rseq = aln.aligned_target_sequence
                    #except:
                    #print("Failure to Align?", len(seq), aln['target_begin'], len(region), aln['query_begin'], aln['cigar']
                        #print("Cant use aligned query sequence object- skipping. len,cigar ", len(aln['query_sequence']), aln['cigar'])
                        #qseq = aln['query_sequence'][aln['query_begin']:aln['query_end']]
                        #continue
                    qseq = aln.query_sequence[aln.query_begin:aln.query_end]
                    qqual = qual[aln.query_begin:aln.query_end]
                    #print('\t'.join([str(v) for v in [len(aln.query_sequence),len(qseq),len(aln.target_sequence),len(rseq),aln.cigar,aln.optimal_alignment_score,aln.optimal_alignment_score/(len(rseq)+1)]]))
                    #if aln.optimal_alignment_score != len(rseq) + 1:
                        #print(rseq.strip(), aln['cigar'], qseq.strip(), qqual.strip(), sep = '\t', end = '\n')
                    if aln.optimal_alignment_score / (len(rseq)+1) < .5: #about .8% of the alignments are below this threshold in quality and shouldn't be used.
                        continue
                    for i, rb in enumerate(rseq):
                        assert len(qseq) == len(rseq)
                        try:
                            qb = qseq[i]
                            if qb in 'ACGT' and phred_code_qual(qqual[i]) > 20: #trying a more stringent cutoff, 20 is .01 # 12: #12 is the cutoff of <5% of being an error, e.g. qscore 13+ == <.05 of being error
                                counts[index + i + 1][qb] += 1
                                #if rb in 'CG' and qb in 'AT':
                                #    print(id[0], id[1], id[2], aln['cigar'], index, i+1, rb, qb, phred_code_qual(qqual[i]))
                        except:
                            print("Error trying to count alignment to region", region, aln.aligned_target_sequence, aln.query_sequence, aln['cigar'])#, index, i, rb, qb)
                        #except:
                            #continue
                            #print("Index issues", 'countslen', len(counts), 'len rseq', len(rseq), 'target start', index, 'current base', i, aln['cigar'], qseq, rseq)
                    #    except:
                     #       print("Index error with counts of length {}, index {}, position in query {}".format(len(counts),index,i))
 #                   try:
#                        qb = aln['query_sequence'][tar_ind - aln['target_begin'] + aln['query_begin']]
  #                  except:
                        #print("Index {} not in sequence of length {}, {} trimmed".format(tar_ind - aln['target_begin'] + aln['query_begin'], len(aln['query_sequence']), aln['query_begin']))
                        #continue
                    #if tar_ind in counts:
                     #   if qb in 'ACGT':
                     #       counts[tar_ind][qb] += 1
                #for tar_ind in range(aln['query_begin'], aln['query_end']):
                #    if index + tar_ind in counts:
                #        base = aln['query_sequence'][tar_ind]
                #        if base in 'ACGT':
                #            counts[index + tar_ind][base] += 1
            #else:
             #   print("ErrorChecker from basecount in make_brpileup: Sequence length overmatched? {}".format(len(seq)), matched, seq)
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
    cons = generate_consensus(reg, splts, id = (coord[0][0],coord[0][1],k)) #version of the reference with either Ns for ambiguous/no information or the 
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
