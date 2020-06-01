#!/usr/bin/env python
from Bio import SeqIO as sqio
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', help = 'path to the gtf')
    parser.add_argument('-g', '--genome', help = 'path to the genome fasta')
    parser.add_argument('-i', '--ids', help = 'file with gene ids to use. Optional', default = None)
    parser.add_argument('-o', '--output', help = 'name of the output file.', default = 'v2.lookupplus')
    parser.add_argument('-d', '--detail', type = bool, help = 'toggle to further break down nonsynonymous mutation into conservative and radical changes. Default False', default = False)
    args = parser.parse_args()
    return args

class Coding:
    '''
    This class contains the methods and information required to go from a bundle of split GTF entries to a translateable full CDS
    Used for synonymity analysis and annotation of mutations.
    Dependent on Biopython.
    '''
    translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    aagroups = {'aliphatic':['G','A','V','L','I'], 'hydroxyl':['S','C','U','T','M'], 'cyclic':['P'], 'aromatic':['F','Y','W'], 'basic':['H','K','R'], 'acidic':['D','E','N','Q'], 'N/A':['Stop']}

    def __init__(self, gtfe):
        '''
        Save the location, gene, and strand information for these entries and assert that they match appropriately.
        '''
        self.gid = gtfe[0][-1].split()[1].strip("';\"") #they should all have the same gene on being passed in.
        self.strand = gtfe[0][6] #they should also all have the same strand.
        self.cds = []
        for ge in gtfe:
            #double check their gene and strand attributes match.
            assert gtfe[0][-1].split()[1].strip("';\"") == self.gid
            assert gtfe[0][6] == self.strand
            #proceed to parse it
            if ge[2] == 'start_codon':
                self.start = (ge[0], int(ge[3]), int(ge[4]))
            elif ge[2] == 'stop_codon':
                self.stop = (ge[0], int(ge[3]), int(ge[4]))
            else:
                #need to retain both the location and the frame.
                self.cds.append([ge[0], int(ge[3]), int(ge[4]), int(ge[7])])
    
    def locate(self):
        '''
        Return a tuple representing the overall start and stop location of the sequence
        Uses the original coordinates from the GFF.
        '''
        if self.strand == '+':
            return (self.start[0], self.start[1], self.stop[2])
        elif self.strand == '-':
            return (self.stop[0], self.stop[1], self.start[2])
        
    def extract_sequence(self, fasta):
        '''
        Extract the stitched coding sequence for each attribute in cds, plus the stop codon, from the target fasta.
        Reverse complements as necessary.
        Asserts that it starts with a start and stops with a stop, otherwise will raise an error.
        Returns the stitched coding sequence as a Biopython.SeqRecord object.
        '''
        genome = sqio.to_dict(sqio.parse(fasta, format = 'fasta'))
        if self.strand == '+':
            seq = self._read_positive(genome)
        elif self.strand == '-':
            seq = self._read_negative(genome)
        #sanity check the sequence
        if seq.seq[:3] != 'ATG':
            print("Sequence {} has an incorrect start {}".format(self.gid, seq.seq[:3]))
        if seq.seq[-3:] not in ['TAA', 'TAG', 'TGA']:
            print('Sequence {} has an incorrect stop {}'.format(self.gid, seq.seq[-3:]))
        if len(seq.seq) %3 != 0:
            print("Sequence {} is not divisible by 3".format(self.gid))
        seq.id = self.gid
        seq.description = 'CDS extracted by Coding custom class'
        seq.name = fasta
        return seq
    
    def _read_positive(self, genome):
        #read the fasta on the positive strand.
        #remember the start codon is included in the cds, but the end codon is not.
        #seqs = [genome[chro][start-1+frame:stop] for chro, start, stop, frame in self.cds]
        seqs = [genome[chro][start-1:stop] for chro, start, stop, frame in self.cds]
        fseq = seqs[0].seq
        for sr in seqs[1:]:
            fseq += sr.seq
        fseq += genome[self.stop[0]][self.stop[1]-1:self.stop[2]]
        return fseq 
    
    def _read_negative(self, genome):
        rcds = self.cds[::-1]
        #seqs = [genome[chro][start-1:stop-frame] for chro, start, stop, frame in rcds]
        seqs = [genome[chro][start-1:stop] for chro, start, stop, frame in rcds]        
        fseq = genome[self.stop[0]][self.stop[1]-1:self.stop[2]]
        for sr in seqs: 
            fseq += sr.seq
        return fseq.reverse_complement()
    
    def get_index(self):
        '''
        This method returns the vector of all indeces for bases that would go into a stitched sequence, in order that they're read.
        Doesn't need a fasta file. 1-based coordinates.
        '''
        idx = []
        if self.strand == '+':
            #start is included in the first cds, so we can pass it.
            for chro, start, stop, frame in self.cds:
                #for i in range(start-1+frame, stop):
                for i in range(start-1, stop):
                    idx.append(i)
            idx.extend([i for i in range(self.stop[1]-1, self.stop[2])])
        elif self.strand == '-':
            idx.extend([i for i in range(self.stop[1]-1, self.stop[2])])
            rcds = self.cds[::-1]
            for chro, start, stop, frame in rcds:
                for i in range(start-1,stop): #-frame):
                    idx.append(i)
        return idx
    
    def annotate(self, genome, detail = False):
        '''
        This method uses the GFF to create an annotation structure representing the effects of all possible changes to this coding sequence.
        Returns a dictionary keyed with (index, ref, alt) that returns one of several possible effect codes.
        '''
        fullseq = self.extract_sequence(genome)
        index = self.get_index()
        #zip these vectors together, then divide them into groups of 3.
        combined = list(zip(index, str(fullseq.seq)))
        codons = [combined[i:i+3] for i in range(0, len(index), 3)]
        complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
        effects = {}
        for c in codons:
            intent = self.translate.get(''.join([cv[1] for cv in c]), None)
            if intent == None:
                print("Can't translate original", c)
                continue
            for i,ib in enumerate(c):
                #i represents the location within the codon
                #gi represents the location in the genome
                #b is the actual base.
                gi, b = ib
                for ab in 'ACGT':
                    if ab != b.upper():
                        nc = [cv[1] for cv in c]
                        nc[i] = ab
                        change = self.translate.get(''.join(nc), None)
                        if change == None:
                            print("Can't translate alternative", nc)
                            continue
                        #finally, compare the intent and the change.
                        #have to complement the bases if this gene is negative-strand before saving the data.
                        if self.strand == '+':
                            key = (int(gi)+1, b, ab) #+1 brings it into line with the 1-based indexing of the pileup
                        elif self.strand == '-':
                            key = (int(gi)+1, complement[b], complement[ab])
                        if intent == change:
                            effects[key] = 'synonymous' #no change, definitely okay
                        elif intent != 'Stop' and change == 'Stop': #stop gained mid-sequence, very bad
                            effects[key] = 'stop_gained'
                        elif intent == 'Stop' and change != 'Stop': #stop lost at the end, also bad
                            effects[key] = 'stop_lost'
                        else:
                            if detail:
                                #check whether the intent and the change are in the same group.
                                intg = [k for k,v in self.aagroups.items() if intent in v][0] #should be length 1 anyways.
                                if change in self.aagroups[intg]:
                                    #its conservative.
                                    effects[key] = 'conservative'
                                else:
                                    #its radical
                                    effects[key] = 'radical'
                            else:
                                effects[key] = 'nonsynonymous' #different amino acid, maybe bad
        return effects

def group_data(gtf_path, ids = None):
    curgene = None
    curgene_ents = []
    with open(gtf_path) as inf:
        for entry in inf:
            spent = entry.strip().split('\t')
            gid = spent[-1].split()[1].strip("';\"")
            if ids != None:
                if gid not in ids:
                    continue
            #print(gid)
            if curgene == None:
                curgene = gid
            if gid == curgene and spent[2] in ['start_codon', 'CDS', 'stop_codon']:
                curgene_ents.append(spent)
            elif gid != curgene:
                curgene = gid
                if spent[0] in ['Scf_2L','Scf_2R','Scf_3L','Scf_3R','Scf_X'] and len(curgene_ents) > 0:
                    yield curgene_ents #skip scaffolds/chr4, for now at least
                curgene_ents = []
                                
def process_annd(annd, detail = False):
    #convert an annd dictoinary into a lookup+ list structure for printing to a file.
    #mirrors the original lookup concept but with additional columns representing stop gain and stop loss
    data = []
    coord_data = {} #first, need to reshape the dictionary to be grouped by coordinate for easier accession
    for k,effect in annd.items():
        coord_data[(k[0],k[1])] = coord_data.get((k[0],k[1]), []) + [[k[2], effect]]
    
    #then collect the bases which have each effect at an index into groups
    for c, ae in coord_data.items():
        syn = ','.join([a for a, e in ae if e == 'synonymous'])
        sg = ','.join([a for a, e in ae if e == 'stop_gained'])
        sl = ','.join([a for a, e in ae if e == 'stop_lost'])
        if detail: 
            con = ','.join([a for a, e in ae if e == 'conservative'])
            rad = ','.join([a for a, e in ae if e == 'radical'])
            data.append([c[0], c[1], syn, con, rad, sg, sl])
        else:
            non = ','.join([a for a, e in ae if e == 'nonsynonymous'])
            data.append([c[0], c[1], syn, non, sg, sl])
    return data
        
def create_lookup_plus(gtf_path, genome, outf = 'dsim_all.lookupplus', ids = None, detail = False):
    with open(outf, 'w+') as out:
        for cd in group_data(gtf_path, ids):
            cobj = Coding(cd)
            if len(cobj.extract_sequence(genome).seq)%3 != 0:
                continue #skip ones that come out weird for whatever reason.
            #if ids != None:
            #    if cobj.gid not in ids:
            #        continue
            print("Examining gene {}".format(cobj.gid))
            annd = cobj.annotate(genome, detail)
            #turn the annd into a sorted list of entries and print.
            dvecs = process_annd(annd, detail)
            for d in dvecs:
                print('\t'.join([str(v) for v in [cobj.start[0]] + d + [cobj.gid, cobj.strand]]), file = out)

def read_geneid(path):
    gids = set()
    with open(path) as inf:
        for entry in inf:
            gids.add(entry.strip())
    return gids

def get_paralogs(gids, species = 'sim', path_to_fb_convert = 'dmel_orthologs_in_drosophila_species_fb_2020_02.tsv'):
    #flybase has a particular file they provide that lists all paralogs.
    #this function returns the set of IDs passed in through GIDs.
    pcv = {}
    with open(path_to_fb_convert) as inf:
        for entry in inf:
            if entry[0] != '#':
                spent = entry.strip().split()
                melid = spent[0]
                spec = spent[6].split('\\')[0]
                if spec == 'Dsim' and spent[0] in gids:
                    pcv[melid] = spent[5]
    return pcv

def main():
    args = argparser()
    if args.detail:
        print("using detailed annotation")
    if args.ids != None:
        gids = read_geneid(args.ids)
        #print(len(gids))
        pcv = get_paralogs(gids) #not going to be reusable as is in the future, lazy for now.
        print("Using {} genes".format(len(pcv)))
        create_lookup_plus(args.annotation, args.genome, args.output, pcv.values(), detail = args.detail)
    else:
        create_lookup_plus(args.annotation, args.genome, args.output, detail = args.detail)


if __name__ == '__main__':
    main()


