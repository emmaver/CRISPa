# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:01:48 2019

@author: emmav
"""

import pandas as pd
import csv
import re
import sys


def parse_CCF(sumfile, crisprfile, casfile):
        
    tsvsum = pd.read_table(sumfile)
    tsvcris = pd.read_table(crisprfile) #min evidence level only applied here 
    
    #amount of arrays
    nbarrays = len(tsvcris.index) #0 if nothing found, but maybe Cas is found anyway? 
    
    #total amount of spacers
    nbspacers = 0
    for i in range(0, len(tsvcris.index)):
        nbspacers += tsvcris.loc[i, 'Spacers_Nb']
    
    #cas: total amount of systems with min. 3 genes, type I-E/I-F/I-C, other types, false positives
    nbcas = 0
    nbie = 0
    nbif = 0
    nbic = 0
    nbo = 0 
    ostr=[]
    
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
                        ostr.append(cctype)
                elif parse2 != None:
                    nbo += int(parse2.group(1))
                        
    nbcas -= nbo
             
    #parse cas report (not a tsv)
    allgenes=[]
    with open(casfile) as cas:
        c=cas.read().split('#'*44)
        for seq in c:
            hits=re.findall('Summary system.+\[(.+)\]',seq)
            if hits != None:
                for hit in hits:
                    genes=list(set(hit.split('; ')))
                    allgenes.append(genes)
                    
    
    crispr = {'nbarrays': nbarrays, 'nbspacers': nbspacers, 'nbcas': nbcas, 'typeIE': nbie, 'typeIF': nbif, 'typeIC': nbic, 'nbother': nbo, 'otherstr': list(set(ostr)), 'allgenes': allgenes}
    return crispr

def write_csv(sample, outfile, crispr):
    result={'sample': sample}
    result.update(crispr)
    with open(outfile, "w") as out:
        writer=csv.DictWriter(out, fieldnames=result.keys())
        writer.writerow(result)
        

sample = snakemake.wildcards.sample
path = '/scratch/leuven/325/vsc32567/MT/crispr_comparison/'
lsr = re.search('/([ls]r|ncbi)/', snakemake.input[1]).group(1)

sumfile = path+'crispr/'+lsr+'/CRISPRCasFinder/'+sample+'/TSV/CRISPR-Cas_summary.tsv'
crisprfile = path+'crispr/'+lsr+'/CRISPRCasFinder/'+sample+'/TSV/Crisprs_REPORT.tsv'
casfile = path+'crispr/'+lsr+'/CRISPRCasFinder/'+sample+'/TSV/Cas_REPORT.tsv'
outputfile = path+snakemake.output[0]

#sumfile = sys.argv[1]
#crisprfile = sys.argv[2]
#casfile = sys.argv[3]
#sample = sys.argv[4]
#outputfile = "{x}_ccf.csv".format(x=sample)

crispr = parse_CCF(sumfile, crisprfile, casfile)
write_csv(sample, outputfile, crispr)



