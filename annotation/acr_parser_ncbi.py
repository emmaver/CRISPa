# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:17:18 2020

@author: emmav
"""

import pandas as pd
import csv
import re
import sys

def parse_acrfinder(filepath):
    #first argument is homology based, second is GBA (just strings, maybe the file doesn't exist)
    HOM=file_path+'_homology_based.out'
    GBA=file_path+'_guilt_by_association.out'
    
    #from homology
    hom_acr_genes=[]
    IE=0
    IF=0
    #IC=0
    
    nb_hits=0
    other_genes=0
        
    #from GBA
    low_acr=0
    
    GBA_IE_med=0
    GBA_IF_med=0
    GBA_IC_med=0
    
    GBA_IE_high=0
    GBA_IF_high=0
    GBA_IC_high=0
    
    phage_prot=[]
    
    
    try:
        hom = pd.read_table(HOM)
        genes_hom = hom['Acr_Hit|pident'].values.tolist()
        
        
        for i in range(0,len(genes_hom)):
            if genes_hom[i] != '---':
                hom_acr_genes.append(genes_hom[i].split('|')[0])  
                
        nb_hits=len(hom_acr_genes)
        other_genes=len(genes_hom)-nb_hits

        
        for gene in hom_acr_genes:
            types=re.findall('I[EFC]',gene)
            for type in types:
                if type == 'IE':
                    IE += 1
                elif type == 'IF':
                    IF += 1
                #elif type == 'IC':
                   # IC += 1
            
        hom_acr = {'nb_IE': IE, 'nb_IF': IF, 'Acr-genes': hom_acr_genes, 'nb_hits': nb_hits, 'nb_other_genes': other_genes}
    
    except FileNotFoundError:
        print('No homology-based results')
        hom_acr = {'nb_IE': IE, 'nb_IF': IF, 'Acr-genes': hom_acr_genes, 'nb_hits': nb_hits, 'nb_other_genes': other_genes}
    
    
    try:
        #GBA file: some lines are problematic due to missing tab
        GBA_new = GBA.replace('.out','_adapted.out')
        with open(GBA) as gba, open(GBA_new,'w') as out:
            text = gba.read()
            new = re.sub('Confidence\t(\d)','Confidence\t\t\\1', text)
            out.write(new)
        
        gba = pd.read_table(GBA_new)
                
        low = gba.loc[gba['#Classification'] == 'Low Confidence']
        medium = gba.loc[gba['#Classification'] == 'Medium Confidence']
        high = gba.loc[gba['#Classification'] == 'High Confidence']
    
        #Acr/Aca -> Acr candidates per confidence level and per type or no type (low)
        low_acr = len(low['Acr/Aca'])
        
        for i in medium.index:
            if medium.loc[i, 'Acr/Aca'] == 'Acr':
                target=medium.loc[i, 'Self Target Outside 5000 BP']
                subtypes=set(list(re.findall('Type([IV]*[A-Z]?)', target)))
                for type in subtypes:
                    if type == 'IE':
                        GBA_IE_med += 1
                    elif type == 'IF':
                        GBA_IF_med += 1
                    elif type == 'IC':
                        GBA_IC_med += 1
        
        for i in high.index:
            if high.loc[i, 'Acr/Aca'] == 'Acr':
                target=high.loc[i, 'Self Target w/in 5000 BP']
                subtypes=set(list(re.findall('Type([IV]*[A-Z]?)', target)))
                for type in subtypes:
                    if type == 'IE':
                        GBA_IE_high += 1
                    elif type == 'IF':
                        GBA_IF_high += 1
                    elif type == 'IC':
                        GBA_IC_high += 1
            
        #MGE/Prophage MetaData
        for hit in gba['MGE/Prophage MetaData'].values.tolist():
            m=re.match('PHAGE_Pseudo_(.+)(\|\d+\|\w+)(\|\w+_[\d.e|-]+)',hit)
            if m != None:
                prot = m.group(1)+m.group(3)
                phage_prot.append(prot)
        
        gba_acr = {'LC': low_acr, 'MC_IE': GBA_IE_med, 'MC_IF': GBA_IF_med, 'MC_IC': GBA_IC_med, 'HC_IE': GBA_IE_high, 'HC_IF': GBA_IF_high, 'HC_IC': GBA_IC_high, 'phage': phage_prot}
    
    except FileNotFoundError:
        print('No CRISPR-Cas systems, thus no GBA results')
        gba_acr = {'LC': low_acr, 'MC_IE': GBA_IE_med, 'MC_IF': GBA_IF_med, 'MC_IC': GBA_IC_med, 'HC_IE': GBA_IE_high, 'HC_IF': GBA_IF_high, 'HC_IC': GBA_IC_high, 'phage': phage_prot}
         
    return hom_acr, gba_acr


def write_csv(sample, outfile, hom_acr, gba_acr):
    result={'sample': sample}
    result.update(hom_acr)
    result.update(gba_acr)
    with open(outfile, "w") as out:
        writer=csv.DictWriter(out, fieldnames=result.keys())
        writer.writerow(result)
        

sample = sys.argv[1]
file_p = sys.argv[2]
s='_'.join(sample.split('_')[:2])
file_path=file_p+s
outputfile = sys.argv[3]

hom_acr, gba_acr = parse_acrfinder(file_path)
write_csv(sample, outputfile, hom_acr, gba_acr)
