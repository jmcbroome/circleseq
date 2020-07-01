#!/usr/bin/env python3
#this script implements the code from "Exploration Round 2", which creates a mutation dataframe from likelihood-annotated vcfs.

import numpy as np
import pandas as pd
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help = 'set to True to print status updates', default = True)
    parser.add_argument('-s', '--snpeff', type = bool, help = 'Indicate whether the annotated vcf has a snpeff annotation that should be included in the frame. Default False', default = False)
    parser.add_argument('-l', '--lookup', help = 'Path to a lookup format file used for custom mutation effect annotation. These files are created by gtf_to_lookup.py from a GTF file for your species. Optional', default = None)
    parser.add_argument('files', nargs = '+', help = 'paths to any number of annotated vcf files to include in the frame.')
    parser.add_argument('-o', '--output', help = 'Name to save the dataframe object to. Default is mutations.tsv', default = 'mutations.tsv')
    parser.add_argument('-m', '--maxdepth', type = int, help = 'Maximum depth to include in the frame. Default 2000', default = 2000)
    parser.add_argument('-n', '--mindepth', type = int, help = 'Minimum depth to include in the frame. Default 20', default = 20)
    parser.add_argument('-t', '--genecount', type = int, help = 'Maximum number of mutations allowed in a single gene. Default 500', default = 500)
    parser.add_argument('-c', '--cluster', type = int, help = 'Minimum distance in basepairs required between neighboring somatic mutations. Does not affect germline mutations. Default 50', default = 50)
    parser.add_argument('-g', '--badgenes', help = 'Path to a file containing undesired gene IDs to exclude. Default is None', default = None)
    parser.add_argument('-a', '--shared', help = 'Filter somatic mutations which are shared at least this many times as possible leaky germline. Default 1', default = 1)
    parser.add_argument('-f', '--frequency', type = int, help = 'Set to an integer value for the minimum number of consensus circles supporting an alterantive allele. Default 1', default = 1)
    parser.add_argument('-k', '--mark_only', type = bool, help = 'Set to True to add columns indicating whether each filter would remove a given mutation. Default False', default = False)
    args = parser.parse_args()
    return args

def parse_snpf(snpf):
    '''takes a snpeff annotation data line, sans the "ANN=" at the start, and returns a series of strings representing its state'''
    #variants are separated by comma.
    altd = {}
    vard = snpf.split(',')
    for vd in vard:
        data = vd.split("|") #fields are delineated by whatever this symbol is called
        #the only ones I really care about are the first three fields
        #these are, respectively, the alternative, the type, and the impact.
        #lazy parsing is just to straight up return these. I can always one-hot encode later if I want.
        altd[data[0]] = (data[1],data[2],data[4]) #type, effect, GID, respectively.
    return altd

def construct_base(files, snpf = False):
    mutdf = {k:[] for k in ['Strain','Stage','SampleNum','Chro','Loc','Ref','Alt','Depth','SampleFreq','PredFreq','Prob']}
    if snpf:
        mutdf.update({k:[] for k in ['Type', 'Impact', "GID"]})
    ##THIS IS CURRENTLY WRITTEN TO WORK WITH SIMULANS
    ##IT WILL NEED TO BE EDITED TO WORK WITH MELANOGASTER, IF MEL ISN'T CONSISTENT WITH THE 2/3 R/L SCHEME
    chro_lookup = {'2L':'NT_479533.1', '2R':'NT_479534.1', '3L':'NT_479535.1', '3R':'NT_479536.1', 'X':'NC_029795.1'} #a dictionary for fixing inconsistent chromosome ID's across files

    chro_lookup.update({v:v for v in chro_lookup.values()}) #return the same id if it's already good
    for fn in files:
        info = fn.split('_')[1]#[:-4]
        #print(fn, info)
        strain = info[0]
        stage = info[1]
        if stage == 'f':
            stage = 'a' #handling some weird file names.
        snum = info[2]
        with open(fn) as inf:
            for entry in inf:
                if entry[0] != '#':
                    chro, loc, x, ref, alts, y, pv, info = entry.strip().split()
                    chro = chro.strip("Scf_").strip("_pilon")
                    if chro not in chro_lookup:
                        continue #not one of the main chromosomes.
                    else:
                        chro = chro_lookup[chro]
                    #further parse the info line.
                    if snpf:
                        if "ANN=" in info: #many lines do not receive annotation.
                            info, snpf = info.split("ANN=")
                            #parse the snpf information.
                            snpfd = parse_snpf(snpf)
                        else:
                            snpfd = {} #need an empty dict for default get later.

                    info = info.split(';')#[:-1] #dump an empty space at the end there.
                    depth = int(info[0][3:])
                    sfs = [float(v) for v in info[1][3:].split(',')]
                    try:
                        pfpd = {i[0]:[float(v[2:].strip("[]")) for v in i[2:].split(',')] for i in info[2:] if i != ''}
                    except:
                        print('QC: issue parsing info {}'.format(info))
                    if pfpd == {}:
                        print(entry)
                        print(info)
                        continue
                    for i, a in enumerate(alts.split(',')):
                        sf = sfs[i]
                        pf, p = pfpd[a]
                        mutdf['Strain'].append(strain)
                        mutdf['Stage'].append(stage)
                        mutdf['SampleNum'].append(snum)
                        mutdf['Chro'].append(chro)
                        mutdf['Loc'].append(loc)
                        mutdf['Ref'].append(ref)
                        mutdf['Alt'].append(a)
                        mutdf['Depth'].append(depth)
                        mutdf['SampleFreq'].append(sf)
                        mutdf['PredFreq'].append(pf)
                        mutdf['Prob'].append(p)
                        if snpf:
                            mutdf['Type'].append(snpfd.get(a, ['noncoding_variant'])[0]) #a simplified use of SNPEff to make annotation easier, i don't care about various splice sites and so forth atm.
                            mutdf['Impact'].append(snpfd.get(a, ['x','LOW'])[1])
                            mutdf['GID'].append(snpfd.get(a, [None, None, 'none'])[2])
    mutdf = pd.DataFrame(mutdf)
    mutdf['SSN'] = (mutdf['Strain'] + mutdf["Stage"] + mutdf['SampleNum'])
    mutdf['Somatic'] = mutdf.SampleFreq < .25 #kind of an arbitrary threshold.
    mutdf['LogPredFreq'] = np.log(mutdf.PredFreq)
    return mutdf

def parse_lookupplus(path):
    lookd = {}
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split('\t')
            #store these results as a dictionary of dictionaries, outer key chro-loc, inner keys each of the bases, gene id, strand
            #when accessing with mutdf, grab the mutation's spot, check the reference, then check the results for the alternative
            subd = {}
            try:
                subd['ref'] = spent[2]
                subd['syn'] = spent[3].split(',')
                subd['non'] = spent[4].split(',')
                subd['sgain'] = spent[5].split(',')
                subd['sloss'] = spent[6].split(',')
                subd['gid'] = spent[7]
                subd['strand'] = spent[8]
                lookd[(spent[0], int(spent[1]))] = subd
            except:
                print("Can't parse entry")
                print(entry)
                continue
    return lookd

def annotate_mutdf(mutdf, lookd, prefix = ''):
    '''
    Add columns to a dataframe representing the custom annotation. 
    One column for discordant or the effect, one column for gene id, and one column for strand of that gene.
    '''
    ##ANOTHER LINE THAT WILL NEED TO BE EDITED TO FIX CHROMOSOME NAMES POTENTIALLY
    chro_lookup = {'NT_479533.1':'Scf_2L', 'NT_479534.1':'Scf_2R', 'NT_479535.1':'Scf_3L', 'NT_479536.1':'Scf_3R', 'NC_029795.1':'Scf_X'}

    effects = []
    gids = []
    strands = []
    for i,d in mutdf.iterrows():
        key = (chro_lookup[d.Chro], int(d.Loc))
        try:
            lmd = lookd[key]
        except KeyError:
            #probably intergenic
            effects.append('Intergenic') #was originally NoneTypes, but that causes problems with saving the annotation.
            gids.append('None')
            strands.append('None')
            continue
        #check for discordancy.
        if lmd['ref'] != d.Ref: #can't trust any of these sites.
            effects.append("discord")
            gids.append(lmd['gid'])
            strands.append(lmd['strand'])
            continue
        #if no discordancy, record the effects for this alternative.
        safeent = []
        for et in ['syn', 'non', 'sgain', 'sloss']:
            if d.Alt in lmd[et]:
                safeent.append(et)
                #should only ever record exactly once, if it ends up short I know I have bad entries somewhere.
        if len(safeent) == 1: #if it didn't record exactly once, its screwy. continue.
            effects.append(safeent[0]) 
            gids.append(lmd['gid'])
            strands.append(lmd['strand'])
    #update mutdf.
    mutdf[prefix+'MyAnn'] = effects
    mutdf[prefix+'GID'] = gids
    mutdf[prefix+'Strand'] = strands
    return mutdf

def get_sites(lookup):
    pairs = []
    for a in 'ACGT':
        for b in 'ACGT':
            if a != b:
                pairs.append((a,b))

    def count_sites(lookup):
        sd = {b:0 for b in pairs}
        nd = {b:0 for b in pairs}
        with open(lookup) as inf:
            for entry in inf:
                spent = entry.strip().split('\t')
                if len(spent) != 9:
                    continue
                if spent[0] != 'Chro':
                    if spent[3] != '':
                        for v in spent[3].split(','):
                            if v != '':
                                key = (spent[2],v)
                                if key in sd:
                                    sd[key] += 1
                        #if spent[2] in sd:
                            #sd[spent[2]] += len([v for v in spent[-2].split(',') if v != ''])
                    if spent[4] != '\n':
                        for v in spent[4].split(','):
                            if v != '':
                                key = (spent[2],v)
                                if key in nd:
                                    nd[key] += 1
                        #if spent[2] in sd:
                            #nd[spent[2]] += len([v for v in spent[-1].strip().split(',') if v != ''])
        return sd, nd
    sd,nd = count_sites('dsim_all.lookupplus')
    sync = sum(sd.values())
    nonc = sum(nd.values())
    return sync, nonc, sd, nd

def create_frame(args):
    #this function is the first part of main, and creates and returns an unfiltered dataframe from any number of input files.
    mutdf = construct_base(args.files, args.snpeff)
    if args.lookup != None:
        lookd = parse_lookupplus(args.lookup)
        mutdf = annotate_mutdf(mutdf, lookd) #this uses the whole proteome with no prefix, though the annotation function can be imported into an interactive environment and the prefix argument used to look at subsets of various genes
    return mutdf

def filter_frame(mutdf, args):
    #this function is the second part of main, and filters down the mutation dataframe using a few different quality filters.

    #filter out densely-clustered mutations
    #this is the most stringent filter by far and goes first.
    keepvec = []
    for chro in mutdf.Chro.value_counts().index:
        subdf = mutdf[mutdf.Chro == chro].reset_index()
        locv = sorted(subdf.Loc.astype(int))
        for i in range(len(locv)):
            if i != 0 and subdf.loc[i,'Somatic']:
                if locv[i] - locv[i-1] >= args.cluster:
                    keepvec.append(True)
                else:
                    keepvec.append(False)
            else:
                keepvec.append(True)
    if args.mark_only:
        mutdf['ClusterF:' + str(args.cluster)] = keepvec #false = removed.
    else:
        mutdf = mutdf[keepvec]
    #filter out genes with many mutations
    somg = mutdf[mutdf.Somatic].GID.value_counts()
    targets = [i for i in somg.index if somg[i] > args.genecount and i != "None"] #intergenic does not count!
    if args.mark_only:
        mutdf['MutGeneF:' + str(args.genecount)] =  (~mutdf.GID.isin(targets)) #false = removed.
    else:
        mutdf = mutdf[~mutdf.GID.isin(targets)]
    #filter out sites that are excessively high depth
    if args.mark_only:
        mutdf['HDepF:' + str(args.maxdepth)] = (mutdf.Depth <= args.maxdepth)
    else:
        mutdf = mutdf[mutdf.Depth <= args.maxdepth]
    #and sites with low depth as well.
    if args.mark_only:
        mutdf['MDepF:' + str(args.mindepth)] = (mutdf.Depth > args.mindepth)
    else:
        mutdf = mutdf[mutdf.Depth >= args.mindepth]
    #filter out mutations belonging to genes you want to exclude, listed in a column text file of gene IDs.
    badgenes = []
    if args.badgenes != None:
        with open(args.badgenes) as inf:
            for entry in inf:
                badgenes.append(entry.strip())
    if args.mark_only and args.badgenes != None:
        mutdf['SGeneF'] = (~mutdf.GID.isin(badgenes)) #false = removed.
    else:
        mutdf = mutdf[~mutdf.GID.isin(badgenes)]
    #and remove somatic sites that are shared with any other individual
    locvc = (mutdf[mutdf.Somatic].Chro + mutdf[mutdf.Somatic].Loc.astype(str)).value_counts() 
    targets_som = set([l for l in locvc.index if locvc[l] > args.shared])
    if args.mark_only:
        mutdf['ShareF:' + str(args.shared)] = (~mutdf.Loc.isin(targets_som))
    else:
        mutdf = mutdf[~mutdf.Loc.isin(targets_som)] #this will also remove germline mutations which are in the same location as any shared somatic mutation, but germline doesn't count towards sharing itself.
    #and filter on the raw number of consensus circles supporting.
    if args.mark_only:
        mutdf['SampleC:'+str(args.frequency)] = (round(mutdf.SampleFreq.astype(float) * mutdf.Depth.astype(int)) >= args.frequency)
    else:
        mutdf = mutdf[(round(mutdf.SampleFreq.astype(float) * mutdf.Depth.astype(int)) >= args.frequency)]
    return mutdf

def main():
    args = argparser()
    mutdf = create_frame(args)
    mutdf = filter_frame(mutdf, args)
    mutdf.to_csv(args.output, sep = '\t')

if __name__ == '__main__':
    main()