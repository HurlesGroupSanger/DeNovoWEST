'''
@queenjobo @ksamocha

16/08/2019

Functions to assign weights to different variant classes

'''

import pandas as pd
import logging as lg
import numpy as np
import itertools
from scipy import stats
import sys

def get_loess_weight_dic(wdicpath,cqs, print_flag):
    '''
    create constrained/unconstrained dictionaries mapping CADD score to LOESS enrichment weight for missense variants
    '''
    wdic = pd.read_table(wdicpath)
    wdic = wdic[wdic.cq.isin(cqs)]
    con = dict(zip("True"+wdic.score.round(3).astype(str),wdic.con))
    uncon = dict(zip("False"+wdic.score.round(3).astype(str),wdic.uncon))
    con.update(uncon)
    conmax = wdic.con[wdic.score == 45]
    unconmax = wdic.uncon[wdic.score == 45]
    return(con,conmax,unconmax)

def get_missense_weights(rates,wdicpath,print_flag):
    '''
    get weights for missense variants
    '''
    cqname = ["missense","missense_variant"]
    wdic,conmax,unconmax = get_loess_weight_dic(wdicpath,cqname, print_flag)

    if 'weight' in rates.columns:
        # weight column has already been defined, so only replace for consequences under study
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.constrained.astype(str) + rates.score.round(3).astype(str), rates['weight'])
    else:
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.constrained.astype(str) + rates.score.round(3).astype(str), 'NaN')

    rates['weight'] = np.where(rates.cq.isin(cqname), rates['weight'].map(wdic), rates['weight'])
    rates['weight'] = np.where((rates.constrained == True)  & (rates.score > 45) & (rates.cq.isin(cqname)), conmax, rates['weight'])
    rates['weight'] = np.where((rates.constrained == False)  & (rates.score > 45) & (rates.cq.isin(cqname)), unconmax, rates['weight'])
    return(rates)
    
def get_nonsense_weights(rates,wdicpath, print_flag):
    '''
    get weights for nonsense variants
    '''
    cqname = ["nonsense","stop_gained"]
    wdic,conmax,unconmax = get_loess_weight_dic(wdicpath,cqname,print_flag)
    if 'weight' in rates.columns:
        # weight column has already been defined, so only replace for consequences under study
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.constrained.astype(str) + rates.score.round(3).astype(str), rates['weight'])
    else:
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.constrained.astype(str) + rates.score.round(3).astype(str), 'NaN')
    rates['weight'] = np.where(rates.cq.isin(cqname), rates['weight'].map(wdic), rates['weight'])
    if print_flag:
        sys.stderr.write('after line two\n{0}\n\n'.format(rates.head()))
    rates['weight'] = np.where((rates.constrained == True)  & (rates.score > 45) & (rates.cq.isin(cqname)), conmax, rates['weight'])
    rates['weight'] = np.where((rates.constrained == False)  & (rates.score > 45) & (rates.cq.isin(cqname)), unconmax, rates['weight'])
    return(rates)

def get_other_weight_dic(wdicpath):
    '''
    get weight dic for non missense variants
    '''
    wdic = pd.read_table(wdicpath)
    wdic = wdic[wdic.cq.isin(["synonymous","splice_lof","inframe","frameshift"])]
    con = dict(zip("True"+wdic.cq,wdic.con))
    uncon = dict(zip("False"+wdic.cq,wdic.uncon))
    con.update(uncon)
    return(con)

def get_other_weights(rates,wdicpath):
    '''
    get weights for non missense variants - just enrichment of entire class as opposed to stratified by CADD
    '''
    othercq = ["synonymous","splice_lof","inframe","frameshift"]
    con= get_other_weight_dic(wdicpath)
    rates.weight.loc[rates.cq.isin(othercq)] = rates.constrained.astype(str) + rates.cq
    rates.weight.replace(con,inplace=True)
    return(rates)

def get_indel_weights(wdicpath):
    wdic = pd.read_table(wdicpath)
    indel_weights = wdic[wdic.cq.str.contains("frame")]
    return(indel_weights)

def get_weights(rates,wdicpath,print_flag):
    '''
    given rates df get weights columns using con, uncon dictionaries mapping to enrichment weights
    '''
    rates = get_missense_weights(rates,wdicpath, print_flag)
    rates = get_nonsense_weights(rates,wdicpath, print_flag)
    rates = get_other_weights(rates,wdicpath)
    indel_weights = get_indel_weights(wdicpath)
    return(rates,indel_weights)
