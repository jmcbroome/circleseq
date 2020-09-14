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
import time

def argparser():
    parser = argparse.ArgumentParser()
    #the script will try all combinations of the entries below and save the results to an output.
<<<<<<< HEAD
    parser.add_argument('-v', '--verbose', action = 'store_true', help = "Use to print status updates.")
=======
    parser.add_argument('-v', '--verbose', action = 'store_true', help = "Print status updates.")
>>>>>>> d44592d81925e24dc02ed313370132ff6645c98e
    parser.add_argument('-a', help = 'initial alpha guess to try', type = float, default = 0.1)
    parser.add_argument('-b', help = 'initial beta guess to try', type = float, default = 10)
    parser.add_argument('-c', '--cutoff', type = int, help = 'Minimum non-zero count value to include in the dataset. Default -1 (no filter)', default = -1)
    parser.add_argument('-d', '--depth', type = int, help = 'minimum depth value, default no filtering', default = 0)
    parser.add_argument('-m', '--method', help = 'method to use', default = 'Nelder-Mead')
    parser.add_argument('-f', '--frequency', type = float, help = 'set a maximum frequency cutoff value. Default is .25', default = .25)
    parser.add_argument('-l', '--tolerance', type = float, help = 'value to use for tolerance', default = .01)
    parser.add_argument('-e', '--error', type = float, help = 'initial guess at error rate', default = .00001)
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

def get_site_integral_noerr(seen, depth, a, b, first):
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
        
def get_site_integral(seen, depth, a, b, first, err):
    #get the likelihood of producing this site from this distibution
    #for now treating e as an add-on to the likelihood of seeing a mutation at a site
    #a b are beta parameters.        
    #return integrate.quad(lambda f: beta.pdf(f, a, b) * binom.pmf(seen, depth, f+err), 0, 1)[0] * first #this version is VERY VERY slow
    return integrate.quad(lambda f: f**(a-1) * (1-f)**(b-1) * (f + err)**seen * (1 - f - err)**(depth - seen), 0, 1)[0] * first

def objective(optparams, sitedata):
    #return the negative log likelihood of the a and b parameters for the given sitedata
    #use a dynamic programming solution to minimize redundant calculation
    #two levels to this; first, sites that have been calculated have their total results saved
    #second, the first term of the simplified equation is precalculated and used for the entire run except for the DcM part, which is globally dynamic and won't be recalculated after the first run.
    integrals = {}
    a,b,err = optparams
    fterm = calc_first_term(a,b)
    a = abs(a)
    b = abs(b)
    err = abs(err)
    #print(a,b,err)
    #print(max(a,0),max(b,0),min(.001, max(0.0000001, err)))
    integral_time = []
    for seen, depth, count in sitedata:
        key = (seen,depth) #theoretically in the count structure there should be no repetition of seen/depth combinations... but I don't feel like rewriting the dynamic programmer right now and this should still work.
        if key not in integrals:
            t = time.time()
            #otherwise, need to perform a more complex integral.
            #integ = get_site_integral(seen, depth, a, b, fterm, err)
            #try a different way to incorporate error.
            #integ = get_site_integral_noerr(seen, depth, max(a,0), max(b,0), fterm)                
            integ = get_site_integral_noerr(seen, depth, a, b, fterm)                
            #add a binomial with an optimizable error parameter
            #stored dynamically, of course.
            #constrain the error term.
            
            #integ += binom.pmf(k = seen, n = depth, p = min(.001, max(0.0000001, err))) #attempting to maximize the sum (minimize the negative sum) of probabilities between these distributions for all sites.
            #integ += binom.pmf(k = seen, n = depth, p = err)
            #integ *= binom.pmf(k = seen, n = depth, p = err)
            #integ = max(integ, binom.pmf(k = seen, n = depth, p = err))
            integral_time.append(time.time() - t)
            integrals[key] = [-integ,count] #initialize entry for new combo.
        else:
            integrals[key][1] += count #increment this count by the number of times this particular depth/count combination appears in the base data
    #print the integral time for quality control.
    #print('QC: Integral Time with error {} is {} mean, {} stdev'.format(err, np.mean(integral_time), np.std(integral_time)))
    #calculate tv by summing all pairs of entries multiplied together
    tv = np.sum([ent[0] * ent[1] for ent in integrals.values()])
    #tv = np.sum([-np.log(get_site_integral(seen, depth, a, b, e)) for seen, depth in sitedata]) #non-dynamic version
    return tv

def optimize(sites, depths, counts, a, b, err, method, tol):
    #convert e to a vector representing per-site error rate. for now, its flat, could be an RVS or a specific genome-based vector later.
    if method == 'Nelder-Mead':
        result = minimize(fun = objective, x0 = np.array([a,b,err]), args = np.array(list(zip(sites,depths,counts))), method = method, tol = tol) #minimize the negative sum of log likelihoods for all sites
    else:   
        result = minimize(fun = objective, x0 = np.array([a,b,err]), args = np.array(list(zip(sites,depths,counts))), bounds = ((0.00000000001, np.inf), (0.000000000001, np.inf), (0.000000000001, 1)), method = method, tol = tol) #minimize the negative sum of log likelihoods for all sites
    return result

def parse_pileup(mind = 0, cutoff = -1, freqcap = .25): 
    #pulls from sys.stdin
    #sites = []
    #depths = []
    sd = {} #updating to use a more memory efficient structure
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
                    #sites.append(int(alts.count(base)))
                    #depths.append(int(depth))
                    kt = (alts.count(base), int(depth))
                    sd[kt] = sd.get(kt,0) + 1
    #return sites, depths
    return sd

def main():
    args = argparser()
    #sites, depths = parse_pileup(args.depth, args.cutoff, args.frequency)
    sdict = parse_pileup(args.depth, args.cutoff, args.frequency)
    if args.verbose:
        print("QC: Unique Entries in site tracker " + str(len(sdict)), file = sys.stderr) 
    #reform the dictionary structure into a trio of lists that can be used to create a numpy array
    #this is a compact represenentation of the count distribution of the original data.
    sites = []
    depths = []
    counts = []
    for sd, c in sdict.items():
        counts.append(c)
        sites.append(sd[0])
        depths.append(sd[1])
    result = optimize(sites, depths, counts, args.a, args.b, args.error, args.method, args.tolerance)
    print(result)

if __name__ == "__main__":
    main()
