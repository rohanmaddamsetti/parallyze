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

from config import SimpleConfig
from genomediff import parse_genomediff, GenomeDiff
from gene import Gene, get_genecoordinates
from parallyze_routines import *
import utils

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from numpy import random
import numpy as np
import operator

# TODO: Update for refactor
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

# TODO: Update for refactor
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

# TODO: Update for refactor
def write_proc6_locus_mut_counts(linesmut):
    header = 'locus_tag; genomes'
    with open('locus_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for row in linesmut:
            locus = row[0]
            genomes = row[1]
            outfp.write('{}; '.format(locus))
            outfp.write(', '.join([str(g) for g in genomes]) + '\n')

# TODO: Update for refactor
def proc1(conf):
    '''simulated solution'''
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

# NOTE: Updated to standards for refactor
def proc4(conf):
    '''
    Calculates dN/dS for all genes using the mutations in the
    provided genomediff files.
    '''
    record = utils.parse_genbank(conf.REF_GENOME)
    
    '''
    genomediffs will be the master dictionary of mutations.
    Each mutation stores the line it came from, and is uniquely
    id'd.
    '''
    genomediffs = {}
    for gd_file in conf.GENOMEDIFF_FILES:
        parse_genomediff(gd_file, record, genomediffs)

    dNdS_counts = calculate_dNdS(genomediffs)

# TODO: Update for refactor
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

# NOTE: Updated for refactor
def proc6(conf):
    '''Mike's: number of lines mutating in this particular gene'''
    record = utils.parse_genbank(conf.REF_GENOME)
    utils.print_genbank_summary(record)

    genomediffs = {}
    for gd_file in conf.GENOMEDIFF_FILES:
        parse_genomediff(gd_file, record, genomediffs)    
   
    counts = mutated_lines_per_gene(genomediffs, conf.snp_types)

    print '\n'.join([str(row) for row in counts])

    #genecoords, total_bases = get_genecoordinates(record)   #why is this commented out? cuz I'm not here yet?
    #write_proc6_locus_mut_counts(genefreqs)

config_keys = { 'REF_GENOME': str,
                'GENOMEDIFF_FILES': str,
                'SYNONYMOUS': bool,
                'NONCODING': bool,
                'PSEUDOGENE': bool,
                'INTERGENIC': bool,
                'REPLICATES': int,
                'GENES_TO_DISPLAY': int }

def main():
    parser = argparse.ArgumentParser(prog='parallyze')
    parser.add_argument('--config', dest='config', 
        help="parallyze config file", default='parallyze.conf')
    parser.add_argument('--procedure', dest='procedure', type=int, 
        default=6, choices=range(1,6))
    parser.add_argument('--synonymous', action='store_true')
    parser.add_argument('--noncoding', action='store_true')
    parser.add_argument('--pseudogene', action='store_true')
    parser.add_argument('--intergenic', action='store_true')
    parser.add_argument('--replicates', type=int)
    parser.add_argument('--genes_to_display', type=int)
    args = parser.parse_args()

    conf = SimpleConfig(config_keys, args.config)
    for fn in conf.GENOMEDIFF_FILES:
        assert os.path.isfile(fn)

    '''
    The user can override most of the options in config file
    with (optional) flags at runtime.
    '''

    if args.synonymous:
        conf.SYNONYMOUS = True
    if args.noncoding:
        conf.NONCODING = True
    if args.pseudogene:
        conf.PSEUDOGENE = True
    if args.intergenic:
        conf.INTERGENIC = True
    if args.replicates and args.replicates > 0:
        conf.REPLICATES = args.replicates
    if args.genes_to_display and args.genes_to_display >= 0:
        conf.GENES_TO_DISPLAY = args.genes_to_display

    '''
    Store the set of mutation types to use, for later reference
    in the program.
    '''

    conf.snp_types = set()
    conf.snp_types.add('nonsynonymous')
    if conf.SYNONYMOUS:
        conf.snp_types.add('synonymous')
    if conf.NONCODING:
        conf.snp_types.add('noncoding')
    if conf.PSEUDOGENE:
        conf.snp_types.add('pseudogene')
    if conf.INTERGENIC:
        conf.snp_types.add('intergenic')

    if args.procedure == 1:
	    proc1(conf)
    elif args.procedure == 2:
	    proc2(conf)
    elif args.procedure == 3:
	    proc3(conf)
    elif args.procedure == 4:
	    proc4(conf)
    elif args.procedure == 5:
        proc5(conf)
    elif args.procedure == 6:
        proc6(conf)

main()

