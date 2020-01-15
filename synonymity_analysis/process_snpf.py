#!/usr/bin/env python3

#import
import argparse
import gffutils
import numpy as np
from sqlite3 import OperationalError
import pandas as pd
from multiprocessing import Pool
import sys
#define functions/classes
translate = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
            'TAT':'Y','TAC':'Y','TAA':'Stop','TAG':'Stop','TGT':'C','TGC':'C','TGA':'Stop','TGG':'W',
            'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G', 
            None:'0'} #if a returned AA has a 0, means the AA at that index is unknown

def parse_annvars(path):
    '''
    Parse a vcf with annotations added by SnpEff to extract all coding sequence changes
    Recorded fields are the gene involved, depth of coverage, type X->Y, the category -sense, the "impact", and whether its at a 0, 1, or 2 site.
    '''
    variants = {}
    novariant_depth = 0
    if path == None:
        inf = sys.stdin
    else:
        inf = open(path)
    for entry in inf:
        if entry[0] != '#':
            var_pres = False
            vcfi = entry.strip().split() #first layer includes location and other basic information
            ref = vcfi[3]
            coord = int(vcfi[1])
            info = vcfi[7].split(";") #second layer includes depth and annotation fields
            depth = int(info[0].strip("DP="))
            if info[-1][:4] == 'ANN=': #if there is annotation data included in info.
                varlist = info[-1].split(',')
                for var in varlist:
                    anndata = var.split('|')
                    try:
                        allele, annotation, putative_impact, gene_name, gene_id, feature_type, feature_id, transcript_biotype, rank, HGVSc, HGVSp, cDNA_poslen, CDS_poslen, Prot_poslen, dist, warnings = anndata
                    except:
                        print(len(anndata))
                        print(info[-1])
                        print(anndata)
                        continue
                    if annotation not in ['upstream_gene_variant','downstream_gene_variant','intergenic_region']: #none of these annotations make sense for the rest of this.
                        # gene_id = gene_id.strip("CHR_START-")
                        alt = allele[-1]
                        if len(alt) != 1:
                            print("Alternate invalid length")
                            print(anndata)
                            continue
                        # assert len(alt) == 1 #should be true and if its not I need to look into those.
                        var_pres = True
                        if gene_id not in variants:
                            variants[gene_id] = []
                        #collect the data for this particular mutation.
                        #process the CDS_poslen field into a 0,1,2 position number.
                        #keep in mind that with a 1 index as with the original field, 0 is the last of the 3 fields (1,2,0)
                        mtype = ref + '->' + alt
                        try:
                            loc, tlen = [int(s) for s in CDS_poslen.split('/')]
                            assert tlen%3 == 0 #This should always be true...
                            if loc%3 == 0: #convert the 1,2,0 to a 0,1,2 order for python.
                                pos = '2' 
                            else:
                                pos = str(loc%3 - 1) 
                        except:
                            pos = '-1' #basically means its not in a coding sequence- its intronic.
                        mutdata = (str(depth), mtype, annotation, putative_impact, pos, str(coord))
                        variants[gene_id].append(mutdata)
            if not var_pres:
                #there were no variants here, count the bases seen and move on.
                novariant_depth += depth
    if path != None:
        inf.close()
    return variants, novariant_depth

def get_mtypes():#helper function.
    for b in 'ACGT':
        for a in 'ACGT':
            if b != a:
                yield b + '->' + a

def combine_annotations(variants, db, novariant_depth, outf, verbose = False):
    '''
    Combine the ontology annotations from the gff database for each gene ID with the set of mutations associated with all genes with that ontology term.
    '''
    #get a list of all the types so I can track them.
    # mtypes = list(get_mtypes())
    # ontdata = {'NonVar':novariant_depth}
    # depthdata = {'NonVar':novariant_depth}
    with open(outf, 'w+') as out:
        colvec = ['GeneID','Ontologies','ExonLen','IntronLen','Start','End']
        colvec2 = ['Coordinate', "Position", 'Type', 'Category', 'Impact', 'Depth', 'Count']
        print('#'+"\t".join(colvec), file = out) #a header line to remind me what goes where.
        print('#'+"\t".join(colvec2), file = out)
        for gene, muts in variants.items():
            genedata = {}
            try:
                feature = db[gene]
                print("Collecting data for {}".format(gene))
            except:
                print("GeneID {} not available in database, continuing".format(gene))
                continue
            onts = feature['Ontology_term']
            #collect basic gene information
            genedata['GeneID'] = gene
            genedata['Ontologies'] = ','.join(onts)
            codelen = 0
            intronlen = 0
            for c in db.children(gene, featuretype = 'exon'):
                codelen += abs(c.end - c.start) #length of exons in this gene.
            for c in db.children(gene, featuretype = 'intron'):
                intronlen += abs(c.end - c.start) #length of introns in this gene.
            genedata['ExonLen'] = codelen
            genedata['IntronLen'] = intronlen
            genedata['Start'] = db[gene].start
            genedata['End'] = db[gene].end
            #collect mutation specific data and depths for average depth calculations.
            genedata['Mutations'] = {}
            for m in muts:
                #format is depth, mtype, category, impact, position
                depth, mtype, cat, imp, pos, coord = m #depth at any bases and overall number of mutations seen is correlated, ofc.
                dat = (coord, pos, mtype, cat, imp, depth)
                if dat not in genedata['Mutations']:
                    genedata['Mutations'][dat] = 1
                else:
                    genedata['Mutations'][dat] += 1 #the number of times this particular mutation was seen in the dataset. Germline will be at all depths, or about half of depths rather for heterozygotes.
            #save the data.
            print('\t'.join([str(genedata[k]) for k in colvec]),file = out) #print gene parameter data
            for m,c in genedata['Mutations'].items():
                print('\t'.join(list(m) + [str(c)]), file = out)


# def save_file(geneds, outf):
#     #write the output of the above to a text file.
#     with open(outf) as out:
#         colvec = ['GeneID','Ontologies','ExonLen','IntronLen','Start','End']
#         colvec2 = ['Coordinate', "Position", 'Type', 'Category', 'Impact', 'Depth', 'Count']
#         print('#'+"\t".join(colvec), file = out) #a header line to remind me what goes where.
#         print('#'+"\t".join(colvec2), file = out)
#         for gened in geneds:
#             print('\t'.join([gened[k] for k in colvec]),file = out) #print gene parameter data
#             for m,c in gened['Mutations'].items():
#                 print('\t'.join(list(m) + [c]), file = out)

        #old code below.
        # for o in onts:
                # ontdata[o] = {k:0 for k in ['Silent','Missense','Nonsense','0','1','2'].extend(mtypes)}
            #now recording number of bases included in this ontology for denominator calculations.
            #create a couple more values for key from the above
            # for m in muts:
                #format is depth, mtype, category, impact, position
                # depth, mtype, cat, imp, pos = m #depth at any bases and overall number of mutations seens is correlated, ofc.
                #use the rest to construct a tuple key and add depth.
                # key = (o, mtype, cat, imp, pos)
                # ontdata[key] = ontdata.get(key,[0,0])
                # ontdata[key][0] += 1
                # ontdata[key][1] += depth
    # print(ontdata)
    # print("Total bases seen: {}".format(novariant_depth))
    # save_df(ontdata, outf)

# def save_df(ontd, outf):
#     #flatten, make a dataframe, and save to csv
#     #outer layer is an ontology term, inner layer is everything else in counts.
#     reform = {k:[] for k in ['Term','Mutation','Category','Impact','Position','Count','Depth']}
#     for key, counts in ontd.items():
#         if key == 'NonVar':
#             #give it a unique set of values for feature space.
#             for k in reform:
#                 if k == 'Depth':
#                     reform[k].append(counts)
#                 if k == 'Count':
#                     # reform[k].append(ontd['Length'])
#                     reform[k].append(0)
#                 else:
#                     reform[k].append("OntParams") #OntParams means this entry contains the total bases covered for all depths (in depths) and the length of bases in this ontology (in counts). 
#         # elif key == 'Length':
#             # continue #covered in the above special entry
#         else:
#             term, mtype, cat, imp, pos = key
#             count, depth = counts
#             reform['Term'].append(term)
#             reform['Mutation'].append(mtype)
#             reform['Category'].append(cat)
#             reform['Impact'].append(imp)
#             reform['Position'].append(pos)
#             reform['Count'].append(count)
#             reform['Depth'].append(depth)
#     df = pd.DataFrame(reform)
#     df.to_csv(outf)

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = True)
    parser.add_argument('-a', '--annotation', help = 'Path to a genome annotation file in gff3 format.')
    parser.add_argument('-m', '--mutation', help = 'Path to a SnpEff annotated vcf. Default is stdin', default = None)
    parser.add_argument('-d', '--database', help = 'Name of the database file to create out of the annotation. If it already exists, will use it instead of overwriting. Default is annotations.db', default = 'annotations.db')
    parser.add_argument('-o', '--output', help = "Name of the output data table file for visualization and downstream processing. Default is rates.txt", default = 'rates.txt')
    # parser.add_argument('-t', '--threads', type = int, help = 'Number of threads allowed. Default 1', default = 1)
    # parser.add_argument('-n', '--ontologies', nargs = '+', help = 'Name specific ontologies to evaluate. Default evaluates all ontologies', default = None)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    try:
        db = gffutils.create_db(args.annotation, dbfn=args.database, force=False, keep_order=False, merge_strategy='merge', sort_attribute_values=False)
    except OperationalError:
        db = gffutils.FeatureDB(args.database)
    if args.verbose:
        print("Parsing variant file")
    variants, novariant_depth = parse_annvars(args.mutation)
    if args.verbose:
        print("Combining variants with ontology database")
    combine_annotations(variants, db, novariant_depth, outf = args.output, verbose = args.verbose)
    # save_file(geneds, outf = args.output)
if __name__ == "__main__":
    main()
