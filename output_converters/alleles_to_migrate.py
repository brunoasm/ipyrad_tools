#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script produces a file for migrate-n in SNP format.
It takes as input a alleles.loci file from ipyrad with phased alleles
A file with population assignments is also be provided 
(following the same format as ipyrad's population assignment file, or a table, see help) 
written by B. de Medeiros, starting in 2018
parts of the code were based on pyrad source code
'''
import argparse
import pandas
import os
from collections import defaultdict
import sys
import numpy as np


#first, parse arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', 
                    help = 'path to *.alleles.loci file')
parser.add_argument('-p', '--popfile', 
                    help = 'path to population assignment file in pyrad format')
parser.add_argument('-t', '--table', 
                    help = 'path to csv table with population information')
parser.add_argument('-s', '--sample-field', 
                    help = 'name of column with sample names in csv table', 
                    default = 'sample')
parser.add_argument('-g', '--population-field', 
                    help = 'name of column with population names in csv table', 
                    default = 'population')
parser.add_argument('-f', '--filter', 
                    help = 'Optional. filter loci to only those found across all populations', 
                    action = 'store_true')
parser.add_argument('-o', '--pop-order', 
                    nargs = '+', 
                    help = 'Optional, list all populations in the order wanted for the output. You may include populations with no sequences in the alleles file.')
parser.add_argument('-q', '--quiet', 
                    help = 'Do not print messages about filtered loci',
                    action = 'store_true')



args = parser.parse_args()


if args.popfile is not None != args.table is not None:
    raise Exception('Provide either popfile or table.')
    

outpath = os.path.basename(args.input).split(".")[0] + ".migrate"



snps_list = [{}]
## read data
with open(args.input,'r') as infile:
    for line in infile:
        if '//' in line:
            snps_list.append({})
        elif line.rstrip(): #skip empty lines, sometimes ipyrad adds them
            sample, seq = line.rstrip().split()
            snps_list[-1][sample] = seq
    #the last element will be empty, remove it:
    #snps_list.pop()
    
    snps_table = pandas.DataFrame(snps_list).transpose()
    

###parse population assignment
if args.popfile:
    with open(args.popfile, 'r') as pop_file:
        taxa = defaultdict(list)
        for line in pop_file:
            pop, sample = [x.rstrip('\s\n') for x in line.split()]
            if sample + + '_0' in snps_table.index:
                taxa[pop].append(sample + '_0')
                taxa[pop].append(sample + '_1')

elif args.table:
    info = pandas.read_csv(args.table)
    taxa = defaultdict(list)
    for i, row in info.iterrows():
        if row.loc[args.sample_field] + '_0' in snps_table.index:
            taxa[row.loc[args.population_field]].append(row.loc[args.sample_field] + '_0')
            taxa[row.loc[args.population_field]].append(row.loc[args.sample_field] + '_1')
            
#print taxa

## filter out loci not found in every population
if args.filter:
    for locus_name,locus in snps_table.iteritems():
        all_pops = True
        for pop in taxa.keys():
            temp = locus.loc[locus.index.isin(taxa[pop]) &
                             -locus.isin([np.nan])]
            if not len(temp): #if any pop has NaNs for all samples, mark as false and break
                all_pops = False
                break
        if not all_pops:
            snps_table.drop(locus_name, axis = 1, inplace = True)
            if not args.quiet:
                sys.stderr.write('locus ' + str(locus_name) + ' removed: not in all populations\n')
            
            

#now let's reduce pop dictionary only to samples with data:
for pop in taxa.keys():
    samples = [samp for samp in taxa[pop] if samp in snps_table.index]
    if samples:
        taxa[pop] = samples
    else:
        del taxa[pop]


#this is for testing only. Keeping only first locus:
#snps_table = snps_table.iloc[:,:20]
#snps_table = snps_table.iloc[:,19:20]
#print snps_table

## write output
with open(outpath, 'w') as outfile:
    #first line is number N pop, N loci and project title
    if args.pop_order:
        npops = str(len(args.pop_order))
    else:
        npops = str(len(taxa))
        
    outfile.write(npops + ' '  + str(snps_table.shape[1]) + \
                  ' ' + os.path.basename(args.input).split(".")[0] + '\n')
    
    #second line is number of sites for each locus. Here, all 1
    outfile.write(' '.join(snps_table.apply(lambda x: max([str(len(a)) for a in x if a is not np.nan]), axis=0)) + '\n')
    
    #now we will loop through populations
    #if an order was provided, we will go with that order
    
    if not args.pop_order:
        pops = taxa.iterkeys()
    else:
        pops = args.pop_order
        
    for pop in pops:
        if pop in taxa.keys(): #if population in sequence file, add sequences
            #for each population, first line:
            #n_individuals pop_name
            smaller_table = snps_table.loc[taxa[pop]]
            #print smaller_table
            n_ind = smaller_table.apply(lambda col: str(smaller_table.shape[0] - sum([x is np.nan for x in col])), axis = 0)
            outfile.write(' '.join(n_ind) + ' ' + pop + '\n')
            
            #remaining lines, samples and SNPS, skipping Ns:
            for LOCUS in snps_table.columns:
                for sample in taxa[pop]:
                    if snps_table.loc[sample,LOCUS] is not np.nan:
                        outfile.write("{:<10}".format(sample[:10]))
                        outfile.write(snps_table.loc[sample,LOCUS])
                        outfile.write('\n')
        else: #if not in sequence file, only add population name and no sequences
            outfile.write('0 ' + pop + '\n')