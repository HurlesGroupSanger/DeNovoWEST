'''
@queenjobo @ksamocha
16/08/2019

Functions to calculate p-value in enrichment DeNovoWEST test

'''


import logging as lg
import pandas as pd
import numpy as np
import itertools
from scipy import stats


def calc_p0(gl,s_obs):
    '''exact calculate P(S >= sobs | N = 0)P(N = 0)'''
    if s_obs == 0:
        p0 = 1*stats.poisson.pmf(0,gl)
    else:
        p0 = 0
    return(p0)

def calc_p1(gl,s_obs,rates):
    '''exact calculate P(S >= s_obs | N = 1)P(N = 1)'''
    p1c = rates['prob'][rates['weight']>=s_obs].sum()/rates['prob'].sum()
    p1 = p1c * stats.poisson.pmf(1,gl)
    return(p1)
    
    
def calc_pn(gl,s_obs,rates,n,nsim,weights_sorted):
    '''simulation to approximate  P(S >= s_obs | N = n)P(N = n)'''
    pndnm = stats.poisson.pmf(n,gl)
    nsim = max([int(round(nsim*pndnm)),100])
    if np.sum(weights_sorted[-n:])<s_obs:
        pscore = 0.0
    elif np.sum(weights_sorted[0:n] >= s_obs):
        pscore = 1.0
    else:
        s = sim_score(gl,s_obs,rates,n,nsim)
        pscore = float(s)/nsim
    pn = pndnm * pscore
    return(pn)

def sim_score(gl,s_obs,rates,n,nsim):
    '''simulate drawing n random mutations from genes based on mutation probabilites 
        and calculating severity score
    '''
    #sample position - sum rates across position and selection position
    posprob = rates.groupby(['pos'])['prob'].sum()
    #randomly select mutations and corresponding weights based on mutation rate probability and sum across genes nsim times
    scores = np.sum(np.random.choice(rates['weight'],(n,nsim),p = rates['prob']/gl),axis = 0)
    #number of times that the gene score is greater than what we observe
    s = np.sum(scores>=s_obs)
    return(s)

def get_pvalue(rates,s_obs,nsim):
    '''P(S>=s_obs) under null mutation rate model'''
    #poisson parameter for gene
    gl = rates['prob'].sum()
    exp = np.sum(rates['prob']*rates['weight'])
    #sort weights for use later
    weights_sorted = np.sort(rates['weight'])
    #calculate prob of seeing as or more extreme score if there are 0 mutations or 1 mutations
    p0 = calc_p0(gl,s_obs)
    p1 = calc_p1(gl,s_obs,rates)
    #keep a running p value
    ptot = p0 + p1
    info = "finished all sims"
    #simulate for 2 to 250 mutations in the gene
    for i in range(2,250):
        # set exception if the observed is less than expected (calculating a 1-way p-value)
        if s_obs < exp:
            ptot = 1
            info = "observed < expected, pvalue set at 1"
            break
        #calculate probability of observing gene score as or more extreme as what we observe | i mutations    
        pi = calc_pn(gl,s_obs,rates,i,nsim,weights_sorted)
        #add to running p value
        ptot = ptot + pi
        #
        #to speed up sims further can uncomment the following lines. This ends the simulation if cross nominal p-value threshold
        #if p value is over threshold then stop going further
        #if ptot > 0.01:
        #    info = "pvalue above reasonable threshold"
        #    break
        #
        #if probability of seeing as or more extreme is small enough (but not 0 which may mean not enough simulations) and
        # the number of mutations is larger than the poisson rate then we can break simulations 
        picdf = 1-stats.poisson.cdf(i,gl)
        if picdf< 10**(-12) and i>gl:
            info = "probability of observing this or more number of mutations is too small"
            break
        #if seeing our observed score is not possible with this number of mutations then break
        if np.sum(weights_sorted[0:i])>s_obs:
            ptot = ptot + picdf
            info ="observed is lower than whole distribution - probability is 1"
            break
    if ptot == 0:
        ptot = 10**(-14)
        info = "p value was 0, set at 10^-14"
    return(ptot,info,exp)
