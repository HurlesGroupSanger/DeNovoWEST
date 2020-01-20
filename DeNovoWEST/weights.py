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
    #separate by missense constraint
    con = wdic[wdic.constrained == True]
    uncon = wdic[wdic.constrained == False]
    #here key is based on shethigh (T/F) + missenseconstrained (T/F) + CADD score rounded to 3 figures
    shethigh_con = dict(zip("True"+ "True" + con.score.round(3).astype(str),con.shethigh))
    shetlow_con = dict(zip("False"+ "True" + con.score.round(3).astype(str),con.shetlow))
    shethigh_uncon = dict(zip("True" + "False" + uncon.score.round(3).astype(str),uncon.shethigh))
    shetlow_uncon = dict(zip("False" + "False" + uncon.score.round(3).astype(str),uncon.shetlow))
    #merge dictionaries into one
    alldic = {**shethigh_con, **shetlow_con, **shethigh_uncon, **shetlow_uncon}
    shethigh_con_max = con.shethigh.max()
    shetlow_con_max = con.shetlow.max()
    shethigh_uncon_max = uncon.shethigh.max()
    shetlow_uncon_max = uncon.shetlow.max()
    caddmax = wdic.score.max()
    return(alldic,shethigh_con_max,shetlow_con_max,shethigh_uncon_max,shetlow_uncon_max,caddmax)

def get_missense_weights(rates,wdicpath,print_flag):
    '''
    get weights for missense variants
    '''
    cqname = ["missense","missense_variant"]
    wdic,shcm,slcm,shum,slum,caddmax = get_loess_weight_dic(wdicpath,cqname, print_flag)
    if 'weight' in rates.columns:
        # weight column has already been defined, so only replace for consequences under study
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.shethigh.astype(str) + rates.constrained.astype(str) + rates.score.round(3).astype(str), rates['weight'])
    else:
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.shethigh.astype(str) + rates.constrained.astype(str) + rates.score.round(3).astype(str), 'NaN')
    #when cadd score exceeds bin replace with max for that subset
    rates['weight'] = np.where(rates.cq.isin(cqname), rates['weight'].map(wdic), rates['weight'])
    rates['weight'] = np.where((rates.constrained == True)  & (rates.shethigh == True) & (rates.score > caddmax) & (rates.cq.isin(cqname)), shcm, rates['weight'])
    rates['weight'] = np.where((rates.constrained == True)  & (rates.shethigh == False) & (rates.score > caddmax) & (rates.cq.isin(cqname)), slcm, rates['weight'])
    rates['weight'] = np.where((rates.constrained == False)  & (rates.shethigh == True) & (rates.score > caddmax) & (rates.cq.isin(cqname)), shum, rates['weight'])
    rates['weight'] = np.where((rates.constrained == False)  & (rates.shethigh == False) & (rates.score > caddmax) & (rates.cq.isin(cqname)), slum, rates['weight'])
    print(rates.head(10))
    return(rates)
    
def get_nonsense_weights(rates,wdicpath, print_flag):
    '''
    get weights for nonsense variants
    '''
    cqname = ["nonsense","stop_gained"]
    wdic,shcm,slcm,shum,slum,caddmax = get_loess_weight_dic(wdicpath,cqname, print_flag)
    if 'weight' in rates.columns:
        # weight column has already been defined, so only replace for consequences under study
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.shethigh.astype(str) + rates.constrained.astype(str) + rates.score.round(3).astype(str), rates['weight'])
    else:
        rates['weight'] = np.where(rates.cq.isin(cqname), rates.shethigh.astype(str) + rates.constrained.astype(str) + rates.score.round(3).astype(str), 'NaN')
    rates['weight'] = np.where(rates.cq.isin(cqname), rates['weight'].map(wdic), rates['weight'])
    if print_flag:
        sys.stderr.write('after line two\n{0}\n\n'.format(rates.head()))
    rates['weight'] = np.where((rates.constrained == False)  & (rates.shethigh == True) & (rates.score > caddmax) & (rates.cq.isin(cqname)), shum, rates['weight'])
    rates['weight'] = np.where((rates.constrained == False)  & (rates.shethigh == False) & (rates.score > caddmax) & (rates.cq.isin(cqname)), slum, rates['weight'])
    return(rates)

def get_other_weight_dic(wdicpath):
    '''
    get weight dic for non missense variants
    '''
    wdic = pd.read_table(wdicpath)
    wdic = wdic[wdic.cq.isin(["synonymous","splice_lof","inframe","frameshift"])]
    shethigh = dict(zip("True"+ "False" + wdic.cq,wdic.shethigh))
    shetlow = dict(zip("False"+ "False" + wdic.cq,wdic.shetlow))
    alldic = {**shethigh, **shetlow}
    return(alldic)

def get_other_weights(rates,wdicpath):
    '''
    get weights for non missense variants - just enrichment of entire class as opposed to stratified by CADD
    '''
    othercq = ["synonymous","splice_lof","inframe","frameshift"]
    con= get_other_weight_dic(wdicpath)
    rates.weight.loc[rates.cq.isin(othercq)] = rates.shethigh.astype(str) + rates.constrained.astype(str) + rates.cq
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
