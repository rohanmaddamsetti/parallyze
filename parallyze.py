#!/usr/bin/env python  
#usage: python parallyze.py
'''
Things to do:
1. change genefreqs and genecoords to use the locus tag because many gene names are 'none'
2. add asserts for list sizes to make sure results are sane
3. change locations of mutations to be dropped to align with included/specified gd mutations (ie, synon, non-synon, non-coding, etc)
(done) 4. change from #muts/gene to #lines/gene
5. analytic solution for SNPs
6. dN/dS ratio
7. get genes that positions fall into. have # of times hit, not just hit/nothit (best method for this? locus tag? gene? (do intergenic now - intragenic later)
8. sum, across all lines, # of mutations per gene, for experimental data (done) and simulated data - divide by reps for avg. 
9. print bar graphs?
10. change lists/dictionaries to dicts/Classes 
'''     

import argparse
import os
import sys

import config
import config_default

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from numpy import random
import numpy as np
import operator
import random
import pandas as pd

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    '''returns object that holds procedure name, ref file name, and diff file names'''
    conf={}
    if config.PROCEDURE in ['1','2','3','4', '5', '6']:  
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb') or ref.endswith('.gbk')
        diffs = file_list(config.GENOME_DIFFS)
        assert len(diffs)!=0
        for diff_file in diffs:   
	        assert diff_file.endswith('.gd')
        conf['ref'] = ref
        conf['procedure'] = config.PROCEDURE
        conf['diffs'] = diffs
        return conf
    else:
        print >>sys.stderr, 'Invalid procedure {p}! Using default data'.format(p=config.PROCEDURE)
        conf['ref'] = config_default.REF_GENOME
        conf['procedure'] = config_default.PROCEDURE
        conf['diffs'] = config_default.GENOME_DIFFS
        return conf

def gene_rank_and_mutate_parameters():
    params = {}
    params['synonymous'] = input("Include synonymous mutations in the analysis? 0 for N, 1 for Y: ")
    params['noncoding'] = input("Include noncoding SNPs in the analysis? 0 or 1: ")
    params['pseudogene'] = input("Include SNPs in psuedogenes in the analysis? 0 or 1: ")
    params['intergenic'] = input("Include intergenic SNPs in the analysis? 0 or 1: ")
    params['number_of_top_genes'] = int(input("How many of the most frequently mutated genes would you like displayed? "))
#   params['generations'] = input("For how many generations did your experiment run? ")
    params['replicates'] = input("How many replicates would you like? (If doing procedure 6, enter any number here) ")
    rank_mut_types = ['nonsynonymous']
    if params['synonymous'] == 1:
        rank_mut_types.append('synonymous')
    if params['noncoding'] == 1:
        rank_mut_types.append('noncoding')
    if params['pseudogene'] == 1:
        rank_mut_types.append('pseudogene')
    if params['intergenic'] == 1:
        rank_mut_types.append('intergenic')
    print
    return rank_mut_types, params['number_of_top_genes'], params['replicates']

def snpcount(diff_dict):
    '''input: 'mutations' i.e., parsed gd files
    returns dictionary of gdfiles, each containing a matrix
    (i.e., dict of dicts) of SNP mutations - to and from base'''
    matrixdict = {}
    for diff_name, mutlist in diff_dict.iteritems(): 
        snpmatrix={'A':{'G':0, 'C':0, 'T':0}, 'G':{'A':0, 'C':0, 'T':0},\
            'C':{'A':0, 'G':0, 'T':0}, 'T':{'A':0, 'C':0, 'G':0}}
        for mutation in mutlist:
            if mutation['mut_type']=='SNP':
                old_base = mutation['old_base']
                new_base = mutation['new_base']
                snpmatrix[old_base][new_base] = snpmatrix[old_base].get(new_base,0) + 1
        matrixdict[diff_name] = snpmatrix
    for filename in matrixdict:
        print 'file:', filename
        print 'to/from:', '\n', str_keyvalue(matrixdict[filename]), '\n'
    return matrixdict

def lines_gene_rank(filenames, params):
    '''input: mutations, from parse_gdfiles
    given user input #x, list # mutated lines/gene'''
    mut_genes = {}
    for fname in filenames: #if a file/line has any number of a particular mut, add 1
                            #so max number for a gene is the number of files
        for mut in filenames[fname]:
            if mut['mut_type'] == 'SNP' and mut['snp_type'] in params:
                for tag in mut['locus_tag']:
                    if tag not in mut_genes:
                        mut_genes[tag] = set()
                    mut_genes[tag].add(fname)
    srt_mut_genes = zip(mut_genes.keys(), mut_genes.values())
    srt_mut_genes = sorted(srt_mut_genes, key=lambda x: len(x[1]), reverse = True)
    print 'Some of the top mutating loci'
    for row in srt_mut_genes[:10]:
        print row[0], len(row[1]), row[1]
    return srt_mut_genes

def snpmutate(matrix, num_replicates, refseq_arr): 
    '''input: matrixdict and refseq as numpy array from snpcount and parse_ref, respectively
    called by get_mut_sites

    returns a dict with original base as key, numpy array as value:
    { 'A': np.array(num_replicates by num_mutatations A to others),
      'T': np.array(num_replicates by num_mutations T to others),
      ...
    }
    Where each array has num_replicates rows and columns corresponding to the number of mutations
    for the original base
    '''
    mut_sites = {}
    #print matrix
    for origbase in matrix:  
        num_muts = sum(matrix[origbase][newbase] for newbase in matrix[origbase])
        #print origbase, num_muts
        sites = np.zeros((num_replicates, num_muts), dtype=int)
        for rep in xrange(num_replicates):
            if rep % 50 == 0 and rep > 0: #prints progress report
                print '\t... rep', rep, 'for base', origbase
            try:
                rep_sites = np.random.choice(np.where(refseq_arr==origbase)[0], size=num_muts, replace=False) #[0] because returns a tuple, and we just want 1st element (list of indices). 
            except ValueError as e:
                print >>sys.stderr, e
                print >>sys.stderr, 'Error getting sites for replicate', rep
            else:
                sites[rep,:] = rep_sites
        mut_sites[origbase] = sites
    #print mut_sites
    return mut_sites
                
def get_mut_sites(matrices, refseq, num_replicates):
    mut_sites = {}
    refseq_arr = np.array([c for c in refseq]) 
    for filename in matrices:
        print '\n** generating mutation sites for', filename
        mut_sites[filename] = snpmutate(matrices[filename], num_replicates, refseq_arr)
    print
    '''
    mut_sites = { 'filename1': { 'A': array(row for reps, columns for mut position indices), 
                                 'T': array(...) ... },
                  'filename2': {...},
                  ...
                }
    '''
    return mut_sites

def write_gene_mut_counts(genecoords, mut_sites):
    header = 'gene, ' + ', '.join([filename for filename in mut_sites])
    with open('sim_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for locus_tag in genecoords:
            cds = genecoords[locus_tag]
            for filename in mut_sites:
                line_muts = 0
                for origbase in mut_sites[filename]:
                    line_muts += ((mut_sites[filename][origbase] >= start) & (mut_sites[filename][origbase] < end)).sum()
                row.append(line_muts)
            outfp.write(', '.join([str(c) for c in row]) + '\n')

'''

def write_gene_mut_counts(genecoords, mut_sites):
    header = 'gene, ' + ', '.join([filename for filename in mut_sites])
    with open('sim_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for start, end, name, tag in genecoords:
            row = [tag]
            for filename in mut_sites:
                line_muts = 0
                for origbase in mut_sites[filename]:
                    line_muts += ((mut_sites[filename][origbase] >= start) & (mut_sites[filename][origbase] < end)).sum()
                row.append(line_muts)
            outfp.write(', '.join([str(c) for c in row]) + '\n')

'''

def write_gd_gene_mut_counts(genecoords, gd_genes):
    header = 'gene, count'
    with open('exp_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for _, _, gene, tag in genecoords:
            muts = [tag]
            if tag in gd_genes:
                muts.append(gd_genes[tag])
            else:
                muts.append(0)
            outfp.write(', '.join([str(c) for c in muts]) + '\n')

def write_proc6_locus_mut_counts(linesmut):
    header = 'locus_tag; genomes'
    with open('locus_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for row in linesmut:
            locus = row[0]
            genomes = row[1]
            outfp.write('{}; '.format(locus))
            outfp.write(', '.join([str(g) for g in genomes]) + '\n')

def dnds_calculate(diff_dict):
    '''input: the output from parsed gd files
    goal: count the SNP mutation types for dN/dS, intergenic, etc
    output: dict of dict of a count'''
#    for fname in mutations:
#        for mut in mutations[fname]:
#            if mut['mut_type'] == 'SNP': 
#                    muttypes[gene_name] = muttypes.get(gene_name, 0) + 1
    '''
        for diff_name, mutlist in diff_dict.iteritems():
        muttype_dict = {'synonymous':0, 'nonsynonymous':0, 'intergenic':0, 'noncoding':0, 'pseudogene':0}
        for mutation in mutlist:
            if mutation['mut_type'] == "SNP":
                muttypes[mut['mut_type'][0]] = muttypes.get(mut['mut_type'],0) + 1
        muttypes[diff_name] = muttype_dict
    for fname in muttypes:
        print 'file:', fname
        print 'SNP mutation type', str.keyvalue(muttypes[fname]), '\n'
    return muttypes
    '''
    muttypes = {}
    for diff_name, mutlist in diff_dict.iteritems():
        muttype_dict = {'synonymous':0, 'nonsynonymous':0, 'intergenic':0, 'noncoding':0, 'pseudogene':0}
        for mutation in mutlist:
            if mutation['mut_type'] == "SNP":
                muttype_dict['mut_type'][0] = muttype_dict.get('mut_type',0) + 1
        muttype_dict['diff_name'] = muttype_dict #change this
    for fname in muttypes:
        print 'file:', fname
        print 'SNP mutation type', str.keyvalue(muttypes[fname]), '\n'
    return muttypes

#def introsteps(): #easily done?

def proc1(conf):
    '''simulated solution'''
    params, topgenes, reps = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords, total_bases = get_genecoordinates(record)
    with open('genecoords.txt', 'wb') as fp:
        fp.write(str(genecoords))
    mut_sites = get_mut_sites(matrices, refseq, reps)
    write_gene_mut_counts(genecoords, mut_sites)
    write_gd_gene_mut_counts(genecoords, genefreqs)

'''old proc3 assumptions: Synonymous mutations are neutral, Mutations are independent of one another, No defects to DNA repair, Mutation rate is constant across the genome, There is only one chromosome'''

def proc4(conf):
    '''dnds calculation'''
    params, topgenes, reps = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    dnds_calculate(mutations)

def proc5(conf):
    '''analytical solution'''
    params, topgenes, reps = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords, total_bases = get_genecoordinates(record)
    with open('genecoords.txt', 'wb') as fp:
        fp.write(str(genecoords))
    mut_sites = get_mut_sites(matrices, refseq, params['replicates']) #last should be reps
    write_gene_mut_counts(genecoords, mut_sites)
    write_gd_gene_mut_counts(genecoords, genefreqs)

def proc6(conf):
    '''Mike's: number of lines mutating in this particular gene'''
    params, topgenes, reps = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords, total_bases = get_genecoordinates(record)   #why is this commented out? cuz I'm not here yet?
    #with open('genecoords.txt', 'wb') as fp:
    #    fp.write(str(genecoords))
    write_proc6_locus_mut_counts(genefreqs)

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(dest='config', help="Config.py file should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref'], '\n','Genome diffs: ', conf['diffs'], '\n'
    if conf['procedure'] == '1':
	    proc1(conf)
    elif conf['procedure'] == '2':
	    proc2(conf)
    elif conf['procedure'] == '3':
	    proc3(conf)
    elif conf['procedure'] == '4':
	    proc4(conf)
    elif conf['procedure'] == '5':
        proc5(conf)
    elif conf['procedure'] == '6':
        proc6(conf)

main()

