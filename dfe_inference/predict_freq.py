#!/usr/bin/env python3
#this script uses the MLE process with an input pileup to produce a predicted alpha and beta parameter for the real distribution of mutation frequencies in the sample. Includes method and tolerance parameters.
#this iteration does NOT incorporate per-site error into the model. 
#however, it does include a basic count filter to ignore sites only seen exactly once (and thus remove error from the model, hypothetically). 
import numpy as np
from scipy.optimize import minimize
from scipy.stats import beta, binom, betabinom
from scipy.special import gammaln, comb
import scipy.integrate as integrate
import argparse
import sys
import math

def argparser():
    parser = argparse.ArgumentParser()
    #the script will try all combinations of the entries below and save the results to an output.
    parser.add_argument('-v', '--verbose', type = bool, help = "Set to True to print status updates. Default True", default = False)
    parser.add_argument('-a', help = 'initial alpha guess to try', type = float, default = 0.1)
    parser.add_argument('-b', help = 'initial beta guess to try', type = float, default = 10)
    parser.add_argument('-c', '--cutoff', type = int, help = 'Minimum non-zero count value to include in the dataset. Default -1 (no filter)', default = -1)
    parser.add_argument('-d', '--depth', type = int, help = 'minimum depth value, default no filtering', default = 0)
    parser.add_argument('-m', '--method', help = 'method to use', default = 'Nelder-Mead')
    parser.add_argument('-f', '--frequency', type = float, help = 'set a maximum frequency cutoff value. Default is .25', default = .25)
    parser.add_argument('-l', '--tolerance', type = float, help = 'value to use for tolerance', default = .00001)
    args = parser.parse_args()
    return args
#define a global variable which will contain all the n choose k values I'll need
#update it as necessary.
dcm = {}

def logg(value):
    if value < 0:
        return gammaln(0) #if it tries to set alpha or beta to negative numbers it can fly off the handle, lock it to 0 instead. Hopefully won't be a problem?
    return gammaln(value)

def calc_first_term(a,b):
    #uses/updates nck to calculate the first term value for the equation, result is passed to 'get site integral'.
    #first term is d c m * gamma(a+b)/(gamma(a)*gamma(b))
    #dcm will be accessed per site, but the other part of the term can be precalculated and used as a constant for the run
    return logg(a+b)-(logg(a)+logg(b))

def get_site_integral(seen, depth, a, b, first):
    #get the likelihood of producing this site from this distibution
    #a b are beta parameters.        
    #calculate it out- equation in whole is d c m * gamma(a+b)/(gamma(a)*gamma(b)) * (gfunc(a+seen)*gfunc(b+depth-seen)) / gfunc(a+b+depth)
    #dcm key should already be in there since its constant for any sample regardless of parameters.
    
    dcv = math.log(dcm[(seen,depth)])
    t1 = logg(a+seen)
    t2 = logg(b+depth-seen)
    t3 = logg(a+b+depth)
    value = dcv + first + t1 + t2 - t3 #more complicated than it needs to be because I was printing some values out for QC.
    #if value > 10000:
    #    print(dcv, first, t1, t2, t3)
    return value
        
def objective(optparams, sitedata):
    #return the negative log likelihood of the a and b parameters for the given sitedata
    #use a dynamic programming solution to minimize redundant calculation
    #two levels to this; first, sites that have been calculated have their total results saved
    #second, the first term of the simplified equation is precalculated and used for the entire run except for the DcM part, which is globally dynamic and won't be recalculated after the first run.
    integrals = {}
    a,b = optparams
    fterm = calc_first_term(a,b)
    print(a,b)
    #print(a,b,fterm)
    for seen, depth in sitedata:
        key = (seen,depth)
        if key not in integrals:
            integ = get_site_integral(seen, depth, a, b, fterm)
            integrals[key] = [-integ,1] #initialize entry for new combo.
        else:
            integrals[key][1] += 1 #increment the count for this value by one instead of recalculating the integral, which is slow.
    #calculate tv by summing all pairs of entries multiplied together
    tv = np.sum([ent[0] * ent[1] for ent in integrals.values()])
    #tv = np.sum([-np.log(get_site_integral(seen, depth, a, b, e)) for seen, depth in sitedata]) #non-dynamic version
    return tv

def optimize(sites, depths, a, b, method, tol):
    #convert e to a vector representing per-site error rate. for now, its flat, could be an RVS or a specific genome-based vector later.
    result = minimize(fun = objective, x0 = np.array([a,b]), args = np.array(list(zip(sites,depths))), method = method, tol = tol) #minimize the negative sum of log likelihoods for all sites
    return result

def parse_pileup(mind = 0, cutoff = -1, freqcap = .25): 
    #pulls from sys.stdin
    sites = []
    depths = []
    for entry in sys.stdin:
        try:
            chro, loc, ref, depth, alts, quals = entry.strip().split()
        except:
            print("Bad Entry, Skipping:", entry)
            continue
        if int(depth) > mind:
            for base in 'ACGT':
                if alts.count(base) > cutoff and alts.count(base)/int(depth) < freqcap:
                    dcm.update({(int(alts.count(base)),int(depth)):comb(N=int(depth),k=int(alts.count(base)), exact=True)})
                    sites.append(int(alts.count(base)))
                    depths.append(int(depth))
    return sites, depths

def main():
    args = argparser()
    sites, depths = parse_pileup(args.depth, args.cutoff, args.frequency)
    if args.verbose:
        print("{} sites parsed, optimizing model".format(len(sites)))
    result = optimize(sites, depths, args.a, args.b, args.method, args.tolerance)
    print(result)

if __name__ == "__main__":
    main()
