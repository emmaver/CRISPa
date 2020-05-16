# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:01:48 2019

@author: emmav
"""

import pandas as pd
import csv
import re


def parse_benchmark(benchmarkfile):
    timestat={}
    time=pd.read_table(benchmarkfile, usecols=['s']) 
    timestat['tavg']=round(time['s'].mean(), ndigits=2)
    timestat['tsd']=round(time['s'].std(), ndigits=2)
    return timestat


def parse_CCF(sumfile, crisprfile):
   
    tsvsum = pd.read_table(sumfile)
    tsvcris = pd.read_table(crisprfile) #min evidence level only applied here 
    
    nbarrays = len(tsvcris.index) #0 if nothing found
    nbcas = 0
    nbie = 0
    nbif = 0
    nbic = 0
    nbo = 0 #other: mostly forbidden or type III or so, these are not counted in nbcas
        
    for i in range(0, len(tsvsum.index)):
        nbcas += tsvsum.loc[i, 'Nb Cas']
        cctypes = tsvsum.loc[i, 'Cas Types/Subtypes'] #if no hits: nan, otherwise: str
        if isinstance(cctypes, str):
            cctypes = cctypes.split(',')
            for typehit in cctypes:
                parse = re.search('Type([IV]*[A-Z]?).*n=(\d)', typehit)
                parse2 = re.search('CAS.\(n=(\d)', typehit)
                if parse != None:
                    cctype = parse.group(1)
                    nb = parse.group(2)
                    if cctype == 'IE':
                        nbie += int(nb)
                    elif cctype == 'IF':
                        nbif += int(nb)
                    elif cctype == 'IC':
                        nbic += int(nb)
                    elif isinstance(cctype, str):
                        nbo += int(nb)
                elif parse2 != None:
                    nbo += int(parse2.group(1))
                        
    nbcas -= nbo
                    
    crispr = {'nbarrays': nbarrays, 'nbcas': nbcas, 'typeIE': nbie, 'typeIF': nbif, 'typeIC': nbic, 'other': nbo}
    return crispr

def write_csv(sample, outfile, timestat, crispr):
    result={'sample': sample}
    result.update(timestat)
    result.update(crispr)
    with open(outfile, "w") as out:
        writer=csv.DictWriter(out, fieldnames=result.keys())
        writer.writerow(result)
        

sample = snakemake.wildcards.sample
path = '/scratch/leuven/325/vsc32567/MT/crispr_comparison/'
benchmarkfile = path+snakemake.input[1] #path/crispr/lr-or-sr-or-ncbi/CRISPRCasFinder/{sample}/benchmark.txt
lsr = re.search('/([ls]r|ncbi)/', snakemake.input[1]).group(1)

sumfile = path+'crispr/'+lsr+'/CRISPRCasFinder/'+sample+'/TSV/CRISPR-Cas_summary.tsv'
crisprfile = path+'crispr/'+lsr+'/CRISPRCasFinder/'+sample+'/TSV/Crisprs_REPORT.tsv'
outputfile = path+snakemake.output[0]

timestat = parse_benchmark(benchmarkfile)
crispr = parse_CCF(sumfile, crisprfile)
write_csv(sample, outputfile, timestat, crispr)
