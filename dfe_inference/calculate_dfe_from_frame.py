#!/usr/bin/env python3
#this script implements the full pipeline from "Round 2 of DFE Estimation"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import dadi
from dadi import Godambe
import math
import Optimize_Functions
import argparse
import Selection
from scipy.stats import binom, gamma
import gffutils

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help = 'set to True to print status updates', default = True)
    parser.add_argument('-s', '--samples', type = int, help = 'Number of samples to create a virtual SFS for. Default 1000', default = 1000)
    parser.add_argument('-d', '--maxdepth', type = int, help = 'Maximum depth to include in the analysis. Default 2000', default = 2000)
    parser.add_argument('-m', '--mindepth', type = int, help = 'Minimum depth to include in the analysis. Default 0', default = 0)
    parser.add_argument('-c', '--cutoff', type = float, help = 'Exclude sites with this sample frequency and higher. Default 1', default = 1)
    parser.add_argument('-r', '--ratio', type = float, help = 'Ratio of S to NS thetas for the target species. Default 3.2', default = 3.2)
    parser.add_argument('-b', '--bootstraps', type = int, help = 'Number of bootstraps to use for Godambe uncertainty analysis. Default 25', default = 25)
    parser.add_argument('-bs', '--bootstrap_size', type = float, help = 'Size of bootstraps to use. Default .1', default = .1)
    parser.add_argument('-n', '--model_name', help = 'Name to save the optimized demographic information under. Default "demog_opt"', default = "demog_opt")
    parser.add_argument('-f', '--frame', help = 'Path to the mutation dataframe produced by the earlier steps in the pipeline. Required.')
    parser.add_argument('-p', '--prefix', help = 'Prefix to use when saving residual and SFS comparison plots. Default "plot"', default = 'plot')
    parser.add_argument('-t', '--type', help = 'Choose a feature to split the mutation dataframe into for multiple runs. Optional', default = None)
    parser.add_argument('-g', '--go', help = 'Set a path to a gffdb file created by gffutils to divide mutations by GO term for analysis. Optional', default = None)
    parser.add_argument('-o', '--mode', help = 'Prefix for a target annotation column to use in the frame. Default no prefix', default = '')
    parser.add_argument('-e', '--germline', type = bool, help = 'Set to True to infer DFE for germline mutations instead of somatic mutations (cutoff between types set by argument). Default False', default = False)
    args = parser.parse_args()
    return args

def get_binom(prob, samples = 10000):
    return binom.rvs(n = samples, p = prob) #random variate instead of mean to add some noise to the plot.

def sfs_from_binomial(mutdf, sub, cutoff = 1, samples = 10000, maxd = 2000, mind = 0, mode = 'MyAnn', germ = False):
    if not germ:
        sfvc = mutdf[(mutdf.SampleFreq < cutoff) & (mutdf.Depth > mind) & (mutdf.Depth < maxd) & (mutdf[mode] == sub) & (mutdf.PredFreq > 1e-6)].PredFreq.apply(get_binom, samples = samples).apply(np.around).value_counts()
    else:
        sfvc = mutdf[(mutdf.SampleFreq >= cutoff) & (mutdf.Depth > mind) & (mutdf.Depth < maxd) & (mutdf[mode] == sub) & (mutdf.PredFreq > 1e-6)].PredFreq.apply(get_binom, samples = samples).apply(np.around).value_counts()
    afs = [sfvc[i] if i in sfvc.index else 0 for i in range(0,samples+1)]
    return dadi.Spectrum(afs)

def get_best_demog(path, tv = -1):
    #the file may contain a series of optimizations. we only want the latest set.
    subopts = []
    curopt = []
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split()
            if spent[0] == 'Model':
                subopts.append(curopt)
                curopt = []
            if spent[0] != 'Model':
                curopt.append([float(spent[2]), float(spent[-2]), [float(v) for v in spent[-1].split(',')]])
    #get the best entry
    subopts.append(curopt) #get the last one.
    best = max(subopts[tv], key = lambda x:x[0])
    return best

def exponential_development(params, ns, pts):
    nu, T = params
    xx = dadi.Numerics.default_grid(pts) 
    phi = dadi.PhiManip.phi_1D(xx)
    zygote = lambda t: 2 * (12000/2) ** (t/T) #2 is the starting number of chromosomes in development, 6000 is where doubling ends
    phi = dadi.Integration.one_pop(phi, xx, nu = zygote, T = 13) #13 generations of doubling initially.
    #this is followed by a period of drift with gradual population expansion to some potentially much larger number.
    #this second stage is what I'm not certain about, that should be optimized, including selection potentially.
    phi = dadi.Integration.one_pop(phi, xx, nu = nu, T = T)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
#the above model should be fit on the synonymous site spectra.

#there's also a version of this model that includes a period of selection *after* the expansion.
def develop_then_select(params, ns, pts):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts) 
    phi = dadi.PhiManip.phi_1D(xx)
    zygote_func = lambda t: 2 * (12000/2) ** (t/T) #2 is the starting number of chromosomes in development, 6000 is where doubling ends
    phi = dadi.Integration.one_pop(phi, xx, nu = zygote_func, T = 13) #13 generations of doubling initially.
    #this is followed by a period of drift with gradual population expansion to some potentially much larger number.
    #this second stage is what I'm not certain about, that should be optimized, including selection potentially.
    phi = dadi.Integration.one_pop(phi, xx, nu = nu, T = T, gamma = gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def make_bootstrap_sfs_binom(mutdf, sub, cutoff = 1, samples = 10000, maxd = 2000, number = 25, size = .1):
    spectrums = [] #not correct plural but "specta" is something else in this notebook
    for i in range(number):
        subdf = mutdf[(mutdf.MyAnn == sub) & (mutdf.SampleFreq < cutoff)].sample(frac = size, replace = True)
        spectrums.append(sfs_from_binomial(subdf, sub, cutoff = cutoff, samples = samples, maxd = maxd))
    return spectrums

def neugamma(mgamma, pneu, alpha, beta):
    mgamma = -mgamma
    if (0 <= mgamma) and (mgamma < 1e-4):
        return pneu/(1e-4) + (1-pneu) * Selection.gamma_dist(-mgamma, alpha, beta)
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * (1-pneu)

def compute_sfs(mutdf, args):
    '''
    The first step in the procedure. Uses binomial draws with the predicted sample frequency to estimate the number of inviduals in a virtual population that would have each mutation, then computes the SFS from that distribution.
    '''
    coln = args.mode + 'MyAnn'
    non_sfs = sfs_from_binomial(mutdf, 'non', cutoff = args.cutoff, maxd = args.maxdepth, mind = args.mindepth, samples = args.samples, mode = coln, germ = args.germline)
    syn_sfs = sfs_from_binomial(mutdf, 'syn', cutoff = args.cutoff, maxd = args.maxdepth, mind = args.mindepth, samples = args.samples, mode = coln, germ = args.germline)
    return non_sfs, syn_sfs

def fit_demography(syn_sfs, args):
    '''
    The second step in the procedure. Fits a dadi demographic model, with hard-coded assumptions about development and growth, to the synonymous mutation set.
    Uses the dadi_pipeline wrapper package.
    '''
    pts = [args.samples + 10, args.samples + 20, args.samples + 30]
    #constrain optimization to Nu > 1 only. My initial round of 6k cells is definitely not shrinking.
    p_labels = "nu, T"
    upper = [1000, 1000]
    lower = [1, .001]
    Optimize_Functions.Optimize_Routine(syn_sfs, pts, str(args.samples), args.model_name, exponential_development, 3, 2, fs_folded = False, param_labels = p_labels, in_upper = upper, in_lower = lower)

def create_spectra(demog_params, theta, args):
    ns = np.array([args.samples])
    pts_l = [args.samples + 10, args.samples + 20, args.samples + 30]
    spectra = Selection.spectra(demog_params, ns, develop_then_select, pts_l = pts_l, int_bounds = (1e-5, 50), Npts = 500, echo = True, mp = True)
    theta_ns = theta * args.ratio
    return spectra, theta_ns

def fit_simple_dfe(spectra, theta_ns, non_sfs):
    sel_params = [.2, 1000]
    lower_bound = [1e-3, 1e-2]
    upper_bound = [10, 50000]
    p0 = dadi.Misc.perturb_params(sel_params, lower_bound = lower_bound, upper_bound = upper_bound)
    popt = Selection.optimize_log(p0, non_sfs, spectra.integrate, Selection.gamma_dist, theta_ns, lower_bound = lower_bound, upper_bound = upper_bound, verbose = len(sel_params), maxiter = 250)
    return popt

def fit_neugamma_dfe(spectra, theta, non_sfs):
    neugamma_vec = np.frompyfunc(neugamma, 4, 1)
    sel_params = [.2, .2, 1000]
    lower_bound = [1e-3, 1e-3, 1e-2]
    upper_bound = [1,1,50000]
    p0 = dadi.Misc.perturb_params(sel_params, lower_bound = lower_bound, upper_bound = upper_bound)
    popt_b = Selection.optimize_log(p0, non_sfs, spectra.integrate, neugamma_vec, theta, lower_bound = lower_bound, upper_bound = upper_bound, verbose = len(sel_params), maxiter = 30)
    return popt_b

def save_results(non_sfs, syn_sfs, spectra, simple_popt, complex_popt, like, theta, demog_params, theta_ns, uncert, prefix, args):
    #collect predicted sfs for each model.
    simple_predicted_non_sfs = spectra.integrate(simple_popt[1], Selection.gamma_dist, theta_ns)
    n, a, s = complex_popt[1]
    complex_predicted_non_sfs = spectra.integrate([a,s], Selection.gamma_dist, theta_ns)
    
    #dfe sfs success graph
    sns.lineplot(x = range(len(simple_predicted_non_sfs)), y= simple_predicted_non_sfs, color = 'blue')
    sns.lineplot(x = range(len(complex_predicted_non_sfs)), y= complex_predicted_non_sfs, color = 'orange')
    sns.lineplot(x = range(len(non_sfs)), y = non_sfs, color = 'black')
    plt.savefig(prefix + "_sfs_comparison.png") #these lines should be similar.
    plt.clf()
    #now calculate residuals for each model.
    dadi.Plotting.plot_1d_comp_Poisson(simple_predicted_non_sfs, non_sfs)
    plt.savefig(prefix + "_simple_residuals.png")
    plt.clf()
    dadi.Plotting.plot_1d_comp_Poisson(complex_predicted_non_sfs, non_sfs)
    plt.savefig(prefix + "_complex_residuals.png")
    plt.clf()
    #print some text summarizing what we have here.
    columns = ['Samples', 'Like', 'Theta', 'DemogParams', 'DemogUncert', 'SimpleDFE', 'ComplexDFE', 'CompNeuProp']
    output = [args.samples, like, theta, demog_params, uncert, list(simple_popt[-1]), [a,s], n]
    print('\t'.join(columns))
    print('\t'.join([str(s) for s in output]))

def perform_analysis(mutdf, prefix, args):
    non_sfs, syn_sfs = compute_sfs(mutdf, args)
    
    if args.verbose:
        print("Fitting demographic parameters...")
    fit_demography(syn_sfs, args)
    like, theta, demog_params = get_best_demog('.'.join([str(args.samples), args.model_name, 'optimized', 'txt']))
    if args.verbose:
        print("Bootstrapping for Godambe uncertainty")
    bootstraps = make_bootstrap_sfs_binom(mutdf, 'syn', samples = args.samples) #for Godambe uncertainty of demographic parameters.
    uncert = Godambe.GIM_uncert(dadi.Numerics.make_extrap_func(exponential_development), [args.samples + 10, args.samples + 20, args.samples + 30], bootstraps, demog_params, syn_sfs)

    if args.verbose:
        print("Calculating spectra")
    spectra, theta_ns = create_spectra(demog_params, theta, args)
    if args.verbose:
        print("Fitting DFE models")
    simple_popt = fit_simple_dfe(spectra, theta_ns, non_sfs)
    complex_popt = fit_neugamma_dfe(spectra, theta_ns, non_sfs)
    #not currently using Godambe uncertainty for dfe parameters, though hypothetically possible, haven't worked out issues with implementation
    if args.verbose:
        print("Collecting final results")
    #plot and save final results.
    save_results(non_sfs, syn_sfs, spectra, simple_popt, complex_popt, like, theta, demog_params, theta_ns, uncert, prefix, args)

def master(args):
    if args.verbose:
        print("Reading and calculating SFS")
    mutdf = pd.read_csv(args.frame, sep = '\t', low_memory = False)
    if args.type != None:
        #subset the mutdf and run everything
        subtypes = mutdf[args.type].value_counts().index
        for st in subtypes:
            if args.verbose:
                print("Evaluating subtype {} of column {}".format(st, args.type))
            subdf = mutdf[mutdf[args.type] == st]
            prefix = args.prefix + '_' + st
            perform_analysis(subdf, prefix, args)
    else:
        perform_analysis(mutdf, args.prefix, args)
    if args.verbose:
        print("Analysis complete.")

def split_go_term(mutdf, args):
    #split the mutdf into subframes stored in a dictionary.
    gdb = gffutils.FeatureDB(args.go)
    god = {}
    for gid in mutdf.GID.value_counts().index:
        terms = gdb[gid].attributes['Ontology_term']
        for t in terms:
            if t not in god:
                god[t] = set()
            god[t].add(gid)
    god = {k:v for k,v in god.items() if 50 < len(v) < 1000} #hard-coded limits to remove outlier terms or terms without enough information.
    tdfs = {}
    for term, gv in god.items():
        subdf = mutdf[mutdf.GID.isin(gv)]
        tdfs[term] = subdf
    return tdfs

def master_go(args):
    if args.verbose:
        print("GO mode: Reading and splitting by terms.")
    mutdf = pd.read_csv(args.frame, sep = '\t', low_memory = False)
    mutdf = mutdf[mutdf.GID != "None"] #no reason to retain noncoding/unannotated mutations for this analysis.
    gd = split_go_term(mutdf, args)
    for term, sdf in gd.items():
        print("Evaluating {} mutations assigned to term {}".format(sdf.shape[0], term))
        if args.type != None: #can still divide by arbitrary categories.
            #subset the sdf and run everything
            subtypes = sdf[args.type].value_counts().index
            for st in subtypes:
                if args.verbose:
                    print("Evaluating subtype {} of column {}".format(st, args.type))
                subdf = sdf[sdf[args.type] == st]
                prefix = args.prefix + '_' + st
                perform_analysis(subdf, prefix, args)
        else:
            perform_analysis(sdf, args.prefix, args)
    if args.verbose:
        print("Analysis complete.")

def main():
    args = argparser()
    if args.go == None:
        master(args)
    else:
        master_go(args)

if __name__ == '__main__':
    main()


