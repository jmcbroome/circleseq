#!/usr/bin/env python3

#script for visualization of the text format file created by the circleseq gene annotation secondary analysis pipeline.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

def get_ratios(translate):
    #for each position and for all positions, count the number of synonymous (silent) and nonsynonymous (other) mutations
    counts = {pos:{'syn':0,'non':0} for pos in [0,1,2,12,'all']}
    for codon, aa in translate.items():
        #permute this codon in all ways possible, record results.
        for i, base in enumerate(codon):
            for abase in 'ACGT':
                if base != abase:
                    ncodon = list(codon)
                    ncodon[i] = abase
                    if translate[''.join(ncodon)] != aa:
                        counts[i]['non'] += 1
                        counts['all']['non'] += 1
                        if i == 0 or i == 1:
                            counts[12]['non'] += 1
                    else:
                        counts[i]['syn'] += 1
                        counts['all']['syn'] += 1
                        if i == 0 or i == 1:
                            counts[12]['syn'] += 1
    for pos, sn in counts.items():
        s = sn['syn']
        n = sn['non']
        sumv = s + n
        for t in ['syn','non']:
            counts[pos][t] = sn[t] / sumv
    return counts

def parse_custom(path):
    genes = []
    with open(path) as inf:
        #initialize values for collection.
        gened = {}
        for entry in inf:
            if entry[0] != '#':
                spent = entry.strip().split()
                if entry[0:4] == "FBgn": #its a gene line
                    if gened != {}: #save the dictionary and reset
                        genes.append(gened)
                    gened = {} 

                    gened['GeneID'] = spent[0]
                    gened['Terms'] = spent[1].split(',')
                    gened['Exonl'] = int(spent[2])
                    gened['Intronl'] = int(spent[3])
                    gened['Start'] = int(spent[4])
                    gened['End'] = int(spent[5])
                    gened['Mutations'] = []
                else: #its a mutation line
                    md = {}
                    md['coord'] = int(spent[0]) #do not currently trust this value.
                    md['pos'] = int(spent[1])
                    md['type'] = spent[2]
                    md['cats'] = spent[3].split("&")
                    md['imp'] = spent[4]
                    md['dep'] = int(spent[5])
                    md['count'] = int(spent[6]) #probably always 1 at this point.
                    gened['Mutations'].append(md)
    return genes
                    
def select_ontology(genevec, ontkey):
    ontd = {"Term":ontkey, 'Genes':[], "Exonl":0, "Intronl":0, "Mutations":[]}
    for gene in genevec:
        if ontkey in gene['Terms']:
            ontd['Genes'].append(gene)
            ontd['Exonl'] += gene['Exonl']
            ontd['Intronl'] += gene['Intronl']
            ontd["Mutations"].extend(gene['Mutations'])
    return ontd
                    
def calculate_rates(data, ratios):
    #calculate a series of different rates for a generic set of mutations and context values.
    #assumption of attributes is only that it has to have exonl, intronl, and mutations. depth is calculated in line.
    adepth = np.mean([d['dep'] for d in data['Mutations']])
    #intronic = # of entries with pos -1 / (average depth * length of introns)
    intronic = len([d for d in data['Mutations'] if d['pos'] == -1]) / (adepth * data['Intronl'])
    #all synonymous = # of entries with synonymous / (length of exons * ratio of all mutations which are silent * average depth)
    allsyn = len([d for d in data['Mutations'] if 'synonymous_variant' in d['cats']]) / (data['Exonl'] * ratios['all']['syn'] * adepth)
    allnon = len([d for d in data['Mutations'] if 'missense_variant' in d['cats'] or 'stop_gained' in d['cats']]) / (data['Exonl'] * ratios['all']['non'] * adepth)
    syn3 = len([d for d in data['Mutations'] if 'synonymous_variant' in d['cats'] and d['pos'] == 2]) / (data['Exonl'] * 1/3 * ratios[2]['syn'] * adepth)
    non12 = len([d for d in [sd for sd in data['Mutations'] if sd['pos'] != 2] if 'missense_variant' in d['cats'] or 'stop_gained' in d['cats']]) / (data['Exonl'] * 2/3 * ratios[12]['non'] * adepth)
    if syn3 != 0 and non12 != 0:
        dnds = non12/syn3 #hopefully, approximately.
    else:
        dnds = None
    overall = len(data['Mutations']) / ((data['Exonl'] + data['Intronl'])*adepth)
    return adepth, intronic, allsyn, allnon, syn3, non12, dnds, overall

def construct_ont_table(alldata, ratios):
    #collect ontologies.
    onts = []
    for gened in alldata:
        onts.extend(gened['Terms'])
    onts = list(set(onts))
    #collect data per ontology and place in rows.
    #row structure will be ontology, # of genes in it, total # of mutations, exon len, intron len, adepth, intronic rate, synonymous at 3 rate, nonsynonymous at 12 rate, dnds (the previous 2), overall.
    structure = {k:[] for k in 
                 ['Term', 'GeneCount', 'MutCount', 'Avgdepth', 'ExonLen', 'IntronLen',
                  'IntronRate', 'SynRate', 'Syn3Rate', 'NonRate', 'Non12Rate', 'DnDs', 'Overall']}
    for o in onts:
        ontdata = select_ontology(alldata, o)
        if len(ontdata['Genes']) > 5:
            adepth, intronic, allsyn, allnon, syn3, non12, dnds, overall = calculate_rates(ontdata, ratios)
            if dnds != None:
                structure['Term'].append(o)
                structure["GeneCount"].append(len(ontdata['Genes']))
                structure['MutCount'].append(len(ontdata['Mutations']))
                structure['Avgdepth'].append(adepth)
                structure['ExonLen'].append(ontdata['Exonl'])
                structure['IntronLen'].append(ontdata['Intronl'])
                structure['IntronRate'].append(intronic)
                structure['SynRate'].append(allsyn)
                structure['Syn3Rate'].append(syn3)
                structure['NonRate'].append(allnon)
                structure['Non12Rate'].append(non12)
                structure['DnDs'].append(dnds)
                structure['Overall'].append(overall)
    return pd.DataFrame(structure)

def construct_gene_table(alldata, ratios):
    structure = {k:[] for k in 
                 ['GeneID', "Terms", 'Start', 'End','MutCount', 'Avgdepth', 'ExonLen', 'IntronLen',
                  'IntronRate', 'SynRate', 'Syn3Rate', 'NonRate', 'Non12Rate', 'DnDs', 'Overall']}
    for gene in alldata:
        adepth, intronic, allsyn, allnon, syn3, non12, dnds, overall = calculate_rates(gene, ratios)
        if dnds != None:
            structure['GeneID'].append(gene['GeneID'])
            structure['Terms'].append(','.join(gene['Terms']))
            structure['Start'].append(gene['Start'])
            structure['End'].append(gene['End'])
            structure['MutCount'].append(len(gene['Mutations']))
            structure['Avgdepth'].append(adepth)
            structure['ExonLen'].append(gene['Exonl'])
            structure['IntronLen'].append(gene['Intronl'])
            structure['IntronRate'].append(intronic)
            structure['SynRate'].append(allsyn)
            structure['Syn3Rate'].append(syn3)
            structure['NonRate'].append(allnon)
            structure['Non12Rate'].append(non12)
            structure['DnDs'].append(dnds)
            structure['Overall'].append(overall)
    return pd.DataFrame(structure, index = structure['GeneID'])

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-i', '--input', help = 'Path to custom data file for visualization')
    parser.add_argument('-p', '--prefix', help = 'String to prefix all saved files.', default = 'circleseq')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    #use the translation table to establish neutral baselines for random selection of mutations from among synonymous and nonsynonymous muation
    alldata = parse_custom(args.input)
    translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    ratios = get_ratios(translate)
    #collect ontology data into a frame.
    ontdf = construct_ont_table(alldata, ratios)
    ontdf = ontdf.dropna()
    for col in ontdf.columns:
        if col != 'Term':
            ontdf['Log'+col] = np.log10(ontdf[col].astype(float))
    genedf = construct_gene_table(alldata, ratios)
    genedf = genedf.dropna()

    for col in genedf.columns:
        if col != 'GeneID' and col != 'Terms':
            genedf['Log'+col] = np.log10(genedf[col].astype(float))

    #insert seaborn calls on the dataframes below.
    sns.scatterplot(x='LogSyn3Rate',y='LogNon12Rate',hue='DnDs',data=ontdf)
    plt.savefig(args.prefix + 'OntologyMutationRates.png')
    plt.clf()
    sns.scatterplot(x='LogSyn3Rate',y='LogNon12Rate',hue='DnDs',data=genedf)
    plt.savefig(args.prefix + "GeneMutationRates.png")
    plt.clf()
    #the above two graphs indicate that synonymous and nonsynonymous rates for each gene, and each ontology, are correlated. They further what proportion of genes fall below or above the DnDs threshold for functional conservation.
    sns.scatterplot(x='LogExonLen',y='LogDnDs',data=ontdf)
    plt.savefig(args.prefix + 'OntologyDnDsVariability.png')
    plt.clf()
    sns.scatterplot(x='LogExonLen',y='LogDnDs',data=genedf)
    plt.savefig(args.prefix + 'GeneDnDsVariability.png')
    plt.clf()
    #the above two graphs indicate whether DnDs for individual genes or ontologies is sufficiently informative and low variance based on size.


if __name__ == "__main__":
    main()
