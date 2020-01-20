'''
16/08/2019 @queenjobo @ksamocha

Script to run enrichment test part of the DeNovoWEST testing framework

For usage see 'python DNE_test.py -h'

make a change
'''

__version__ = 1.2
__date__ = '2019-01-09'


#IMPORTS--------------

import logging
import random
import argparse
import pandas as pd
import numpy as np
import itertools
import datetime
from scipy import stats
from probabilities import *
from weights import *

#set seed
random.seed(31)
np.random.seed(31)

#get todays date
now = datetime.datetime.now()
today = str(now).split(' ')[0]

#default paths
wdicpath = "input/loess_weight_dic_v4.tab"
denovopath = "input/de_novos.txt.gz"
outpath = "DeNovoWEST_enrichment_results_" + today + ".tab"


#FUNCTIONS---------------

def get_options():
    '''
    parse options
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--weightdic', default=wdicpath,
        help='path to dictionary of weights to assign variants')
    parser.add_argument('--rates',
        help='path to site specific rates per gene annotated with CADD and constraint info')
    parser.add_argument('--denovos', default=denovopath,
        help='path to denovos annotated with constraint and CADD')
    parser.add_argument('--output',default = outpath,help = "output file to save to")
    parser.add_argument('--nsim',default = 10000000,type = int,help = 'number of simulations for each mutation score')
    parser.add_argument('--nmales',type = int, help = 'number of males in study')
    parser.add_argument('--nfemales', type = int, help = 'number of females in study')
    return parser.parse_args()

    
def load_denovos(path):
    ''' 
    get a DataFrame of DDD de novos
    '''
    logging.info("Loading de novos")
    denovos = pd.read_table(path, sep='\t') # added sep
    denovos['chrom'] = denovos['chrom'].astype(str)
    denovos['pos'] = denovos['pos'].astype(int)
    denovos = fix_consequences(denovos)
    return(denovos)


def fix_consequences(denovos):
    '''
    function to consolidate consequences across cohorts
    '''
    cdic = {'cq':{"frameshift_variant":"frameshift","inframe_insertion":"inframe","inframe_deletion":"inframe","missense_variant":"missense","stop_gained":"nonsense","synonymous_variant":"synonymous","splice_acceptor_variant":"splice_lof","splice_donor_variant":"splice_lof","splice_region_variant":"splice_region","conserved_exon_terminus_variant":"splice_lof"}}
    denovos = denovos[denovos.consequence.isin(cdic['cq'].keys())]
    denovos['cq'] = denovos['consequence']
    denovos = denovos.replace(cdic)
    return(denovos)


def load_rates(ratespath):
    '''
    load rates file
    '''
    #rates = pd.read_table(ratespath,names = ['symbol','chrom','pos','ref','alt','cq','prob','raw','score','maf','constrained'])
    logging.info("Loading rates")
    rates = pd.read_table(ratespath)
    rates['chrom'] = rates['chrom'].astype(str)
    rates['pos'] = rates['pos'].astype(int)
    return(rates)

    
def correct_for_x_chrom(rates, male_n, female_n):
    ''' 
    correct for X chromosome
    '''
    autosomal = 2 * (male_n + female_n)
    female_transmit = male_n + female_n
    male_transmit = female_n
    
    # get scaling factors using the alpha from the most recent SFHS (Scottish
    # Family Health Study) phased de novo data.
    alpha = 3.4
    male_k = 2 / (1 + (1 / alpha))
    female_k = 2 / (1 + alpha)
    
    # correct the non-PAR chrX genes for fewer transmissions and lower rate
    # (dependent on alpha)
    chrX = rates['chrom'].isin(['X', 'chrX'])
    x_factor = ((male_transmit * male_k) + (female_transmit * female_k)) / autosomal
    x_factor = pd.Series([x_factor] * len(chrX), index=rates.index)
    x_factor[~chrX] = 1
    
    rates['prob'] *= x_factor
    
    return rates

    
def get_expected_rates(rates, male, female):
    '''
    get expected mutation rates 
    '''
    logging.info("get expected")
    autosomal = 2 * (male + female)
    rates['prob'] = rates['prob'] * autosomal
    
    return correct_for_x_chrom(rates, male, female)

    
def test_gene(rates,denovos,gene,nsim,indel_weights):
    '''
    test single gene given rates object and de novos
    '''
    generates = rates[rates.symbol == gene]
    generates = get_indel_rates(generates,indel_weights)
    # here add inframe/frameshift rates
    s_obs = denovos[denovos.symbol == gene].weight.sum()
    pval,info,s_exp = get_pvalue(generates,s_obs,nsim)
    return(s_exp,s_obs,pval,info)

    
def get_indel_rates(generates,indel_weights):
    '''
    get gene specific indel rates
    '''
    missense_rate_con = generates[(generates.cq == "missense") & generates.constrained].prob.sum()
    missense_rate_uncon = generates[(generates.cq == "missense") & ~generates.constrained].prob.sum()
    nonsense_rate = generates[generates.cq == "nonsense"].prob.sum()
    inframe_rate_con = missense_rate_con * 0.03
    inframe_rate_uncon = missense_rate_uncon * 0.03
    frameshift_rate = nonsense_rate * 1.3
    inframe_weight_con = float(indel_weights[indel_weights.cq == "inframe"].con)
    inframe_weight_uncon = float(indel_weights[indel_weights.cq == "inframe"].uncon)
    constraint = generates.constrained.any()
    if constraint:
        frameshift_weight =float(indel_weights[indel_weights.cq == "frameshift"].con)
    else:
        frameshift_weight =float(indel_weights[indel_weights.cq == "frameshift"].uncon)
    indelrates = pd.DataFrame([["inframe",inframe_rate_con,inframe_weight_con,True],["inframe",inframe_rate_uncon,inframe_weight_uncon,False],["frameshift",frameshift_rate,frameshift_weight,constraint]],columns = ["cq","prob","weight","constrained"])
    generates = generates.append(indelrates,ignore_index = True)
    return(generates)

    
def test_all_genes(rates,denovos,nsim,indel_weights):
    '''
    test all genes for enrichment
    '''
    logging.info("Starting tests:  "+ str(nsim) + " simulations")
    genes = rates.symbol.unique()
    gnes = denovos.symbol.unique()
    results = []
    for gene in genes:
        logging.info(gene)
        #skip gene if no observed de novos in our dataset
        if denovos[denovos.symbol == gene].weight.sum() == 0:
            continue
        s_exp,s_obs,pval,info = test_gene(rates,denovos,gene,nsim,indel_weights)

        # extract hgnc_id
        hgnc = denovos.loc[denovos['symbol'] == gene, 'hgnc_id'].iloc[0]
        results.append((gene,hgnc,s_exp,s_obs,pval,info))
    return(results)

    
def main():
    args = get_options()

    #initialise log file
    logging.basicConfig(filename=args.output.split(".")[0]+".log",level=logging.DEBUG)

    #number of females/males in DDD
    male = args.nmales
    female = args.nfemales

    #load rates and get expected
    rates = load_rates(args.rates)
    rates = get_expected_rates(rates, male, female)

    #get weights
    logging.info("get weights")
    rates,indel_weights = get_weights(rates,args.weightdic,False)
    rates = rates[rates.weight.apply(np.isreal)]
    
    #load denovos
    denovos = load_denovos(args.denovos)

    #use de novo mutation weights
    denovos,_ = get_weights(denovos,args.weightdic,False)
    denovos = denovos[denovos.weight.apply(np.isreal)]

    # run gene tests
    results = test_all_genes(rates,denovos,args.nsim,indel_weights)

    #save results to file
    df = pd.DataFrame.from_records(results,columns = ["symbol","hgnc_id","expected","observed","p-value","info"])
    df.to_csv(args.output,sep = "\t",index = False)

#------------------SCRIPT--------------------------------------

if __name__=='__main__':
    main()
    
    
    
