import argparse

#this script is intended to parse in a pileup and an annotated vcf and/or a bed
#it will subset to mutations annotated in the vcf and/or sites included in a bed file span

def parse_bed(bedpath):
    bedd = {}
    with open(bedpath) as inf:
        for entry in inf:
            spent = entry.strip().split()
            if spent[0] not in bedd:
                bedd[spent[0]] = []
            bedd[spent[0]].append((int(spent[1]),int(spent[2])))
    ##sort the bed entries. (They should be pre-sorted atm)
    #beddsort = {}
    #for k,v in bedd.items():
    #    beddsort[k] = sorted(v)
    return bedd

def parse_vcf(vcfpath):
    vcfd = {}
    with open(vcfpath) as inf:
        for entry in inf:
            if entry[0] != '#':
                spent = entry.strip().split()
                #parsing that end thing is going to be a pain.
                #use a simple string search for the phrase I want and mark it that way. 
                key = (spent[0], int(spent[1]))
                if 'synonymous_variant' in spent[-1]:
                    vcfd[key] = vcfd.get(key, []) + ['synonymous']
                if 'missense_variant' in spent[-1]:
                    vcfd[key] = vcfd.get(key, []) + ['nonsynonymous']
    return vcfd

def parse_pileup(pilepath, bedd, vcf = None):
    synents = []
    nonents = []
    basekeep = False #state-based, switch in a pairwise iteration method with bed intervals, faster than checking each
    curgene_ind = 0
    curchrom = None
    with open(pilepath) as pilin:
        for entry in pilin:
            spent = entry.strip().split()
            chro = spent[0]
            #if we're onto a new chromosome, we need to reset the curgene_ind
            if chro != curchrom:
                curgene_ind = 0
                curchrom = chro
                print("Subsetting chromosome {}".format(chro))
            loc = int(spent[1])
            #below if statement is no good, way way way way too slow- use a state based pairwise iteration instead
            ###if any([e[0] <= int(loc) <= e[1] for e in bedd[chro]]): #this is not going to be fast, unfortunately.
            #instead, we check if this base has reached the start of the 'current gene span' as noted in bedd[chro][curgene_ind]
            #if it is, we set the basekeep state to true and record until we are no longer in the current span then we set it to false and incremement the curgene ind
            #preprocessing sorted the bed and merged overlapping entries so hopefully there will be minimal weirdness.
            try:
                cbednt = bedd[chro][curgene_ind]
            except:
                continue #there are no more genes in this chromosome, presumably. Just iterate through entries until you get out.
                #alternatively, the chromosome the pileup is looking at may not be one of the main ones, which we also want to skip.
            if int(loc) >= cbednt[1]: #the end statement goes first in a case where I somehow get out ahead of a span without setting basekeep appropriately
                basekeep = False
                curgene_ind += 1

            elif int(loc) >= cbednt[0]:
                basekeep = True

            if basekeep:
                if vcf != None:
                    if (chro,loc) in vcf:
                        if 'synonymous' in vcf[(chro,loc)]:
                            synents.append(entry.strip())
                        if 'nonsynonymous' in vcf[(chro,loc)]:
                            nonents.append(entry.strip())
                    else: #its not a variant site but it is in a gene, so it should go in both since it can probably make both kinds theoretically. 
                        #this isn't necessarily true, but its true often enough that I'm going to shortcut this for right now.
                        synents.append(entry.strip())
                        nonents.append(entry.strip())
                else:
                    synents.append(entry.strip())
    if vcf != None:
        return synents, nonents
    else:
        return nonents

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--annotation', help = 'Path to a SNPEff annotated vcf. Optional', default = None)
    parser.add_argument('-b','--bed', help = 'Path to a bed file containing regions to be retained')
    parser.add_argument('-p', '--pileup', help = 'Path to a pileup for subsetting')
    parser.add_argument('-o', '--output', help = 'Prefix for output files', default = 'filtered')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    bedd = parse_bed(args.bed)
    if args.annotation != None:
        vcfd = parse_vcf(args.annotation)
        print("VCF parsed")
        synents, nonents = parse_pileup(args.pileup, bedd, vcfd)
        with open(args.output + '_synonymous.txt', 'w+') as outf:
            for ge in synents:
                print(ge, file = outf) 
        with open(args.output + '_nonsynonymous.txt', 'w+') as outftwo:
            for ge in nonents:
                print(ge, file = outftwo) 
    else:
        synents = parse_pileup(args.pileup, bedd)
        with open(args.output + '.txt', 'w+') as outf:
            for ge in synents:
                print(ge, file = outf)
    
if __name__ == '__main__':
    main()