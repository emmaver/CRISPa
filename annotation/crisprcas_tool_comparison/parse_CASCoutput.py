# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:45:13 2019

@author: emmav
"""

import pandas as pd
import csv
import os
import re


def parse_benchmark(benchmarkfile):
    timestat = {}
    time = pd.read_table(benchmarkfile, usecols = ['s']) 
    timestat['tavg'] = round(time['s'].mean(), ndigits = 2)
    timestat['tsd'] = round(time['s'].std(), ndigits = 2)
    return timestat


def parse_CASC(report):
    crispr = {}
        
    if not os.path.exists(report):
        crispr = {'nbarrays': 0, 'nbcas': 0}
    else:
        #res = pd.read_table(results)
        #crispr['nbarrays'] = res.shape[0]
        #cas = res['code'].tolist()
       #nbcas = 0
        #for i in range(0,len(cas)):
            #remove hit if no cas hit, no repeat hit & no proper statistics
          #  if cas[i]  == 0:
          #      res.drop(i, axis = 0, inplace = True)
          #  if cas[i] in range(4,8):
           #     nbcas += 1
        #crispr['nbcas'] = nbcas 
        
        rep = open(report).read()
        arrays = re.search('Bona fide CRISPR arrays = (\d+)', rep).group(1)
        crispr['nbarrays'] = arrays
        
        cas = re.search('Arrays with Cas protein upstream = (\d+)', rep).group(1)
        crispr['nbcas'] = cas
    return crispr


def write_csv(sample, outfile, timestat, crispr):
    result = {'sample': sample}
    result.update(timestat)
    result.update(crispr)
    subtype = {'typeIE': None, 'typeIF': None, 'typeIC':None}
    result.update(subtype)
    with open(outfile, "w") as out:
        writer = csv.DictWriter(out, fieldnames = result.keys())
    
        writer.writerow(result)
    
    
    #columns: sample, runtime, runtime-sd, #arrays, #Cas, type I-E, type I-F, type I-C


sample = snakemake.wildcards.sample
path = '/mnt/blastdb/emma/crispr_comparison/'
benchmarkfile = path+snakemake.input[1] #path/crispr/lr-or-sr/CASC/{sample}/benchmark.txt
lsr = re.search('/([ls]r)/', snakemake.input[1]).group(1)
if lsr == 'sr':
    asmfile = 'contigs'
else:
    asmfile = 'assembly'

report = path+'crispr/'+lsr+'/CASC/'+sample+'/'+asmfile+'.report.md'
#result = path+'crispr/'+lsr+'/CASC/'+sample+'/'+asmfile+'.result.txt'
outputfile = path+snakemake.output[0]

timestat = parse_benchmark(benchmarkfile)
crispr = parse_CASC(report)
write_csv(sample, outputfile, timestat, crispr)

