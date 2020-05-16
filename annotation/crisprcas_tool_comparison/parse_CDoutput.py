# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:05:58 2019

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

def parse_CD(detectfile):
    
    with open(detectfile) as detect:
        #empty file if nothing is found
        lines = detect.read() 
        hits = lines.split('\n//\n\n\n')
        if hits[-1] == '': #last one is empty string, also if file was empty
                del hits[-1]
        if len(hits) == 0:
            crispr = {'nbarrays': 0, 'nbcas': 0, 'typeIE': None, 'typeIF': None, 'typeIC': None}
        else:
            nbcas = 0
            nbie = 0
            nbif = 0
            nbic = 0
            for i in range(0, len(hits)):
                if re.search('Questionable array : ([A-Z]+)', hits[i]).group(1) == 'YES':
                    del hits[i]
                nbcas += int(re.search('Score Detail : 1:(\d)', hits[i]).group(1))
                cctype = re.search('Array family : ([A-Z\-]{2,3})', hits[i]).group(1)
                if cctype == 'I-E':
                    nbie += 1
                elif cctype == 'I-F':
                    nbif += 1
                elif cctype == 'I-C':
                    nbic += 1
            if nbie == 0 & nbif == 0 & nbic == 0: #means that type is not detected -> no output
                crispr = {'nbarrays': len(hits), 'nbcas': nbcas, 'typeIE': None, 'typeIF': None, 'typeIC': None}
            crispr = {'nbarrays': len(hits), 'nbcas': nbcas, 'typeIE': nbie, 'typeIF': nbif, 'typeIC': nbic}
        return crispr
    

def write_csv(sample, outfile, timestat, crispr):
    result={'sample': sample}
    result.update(timestat)
    result.update(crispr)
    with open(outfile, "w") as out:
        writer=csv.DictWriter(out, fieldnames=result.keys())
        writer.writerow(result)
   

sample = snakemake.wildcards.sample
path = '/mnt/blastdb/emma/crispr_comparison/'
benchmarkfile = path+snakemake.input[1] #path/crispr/lr-or-sr/CRISPRDetect/{sample}/benchmark.txt
lsr = re.search('/([ls]r)/', snakemake.input[1]).group(1)

detectfile = path+'crispr/'+lsr+'/CRISPRDetect/'+sample+'/'+sample
outputfile = path+snakemake.output[0]

timestat = parse_benchmark(benchmarkfile)
crispr = parse_CD(detectfile)
write_csv(sample, outputfile, timestat, crispr)


