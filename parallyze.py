#!/usr/bin/env python  
#usage: python parallyze.py
'''
Things to do:
1. change genefreqs and genecoords to use the locus tag because many gene names are 'none'
2. add asserts for list sizes to make sure results are sane
3. change locations of mutations to be dropped to align with included/specified gd mutations (ie, synon, non-synon, non-coding, etc)
(done) 4. change from #muts/gene to #lines/gene
5. analytic solution for SNPs
6. graph for dN/dS ratio
7. get genes that positions fall into. have # of times hit, not just hit/nothit (best method for this? locus tag? gene? (do intergenic now - intragenic later)
8. sum, across all lines, # of mutations per gene, for experimental data (done) and simulated data - divide by reps for avg.? 
9. print bar graphs? 
10. condense/reformat/renumber/rearrage procedures
'''     

import argparse
import os
import sys
from pprint import pprint
import inspect
import time
from numpy import random
import numpy as np
import operator
import itertools
from copy import deepcopy

from config import SimpleConfig
from genomediff import parse_genomediff, GenomeDiff
from gene import Gene, get_genecoordinates
import utils
from parallyze_routines import *

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO
from Bio import Phylo

# TODO: Update for refactor
def simulationRM(conf): #old proc1
    #TODO: 
    '''simulated solution'''

    ## Rohan's simple solution:
    snp_total = BasicSNPCount(conf) ## just return the number of dN in genes.
    ## calculate statistics of parallel evolution.
    Statisticulate(conf, snp_total, star=True,reps=100)
    
    ''' 
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
    '''

'''Old Procedure 3: Genomes from multiple isolates from the same experimental evolution population.

1) A phylogeny must be provided as input.
2) Infer genotypes of all internal nodes by "using parsimonious assumptions" -- or better.
3) Count the number x of inferred mutations, and generate a 4x4 matrix of mutation probabilities.

        for 1 to N replicates:
        drop x mutations onto reference genome, and count number of independent mutations per gene.
        
        average the results to calculate the null distribution.

    ## PROBLEM: A true star phylogeny for replicate lineages is NOT recapitulated.
    ## This means that identical parallel mutations at the nucleotide level will
    ## be missed. Phylogeny software that assumes all samples come from the same
    ## time point will allow make errors through assumption.
    ## So, using outside packages is not appropriate for analyzing experimental
    ## evolution phylogenies at this point in my thinking.

    New Procedure 3: replicate the analysis in Lieberman 2011.

'''

def simulationEJB(conf): #old proc2
    '''EJB code to create histogram comparing actual and
         simulated genes mutated/lineage'''

    '''number of lines mutating in this particular gene'''


    snp_total = BasicSNPCount(conf) ## just return the number of dN in genes.
    print snp_total

def proc3(conf):

    print "IN PROC 3"
    #aln = SNPsToAlignment(conf)
    #phy = AlignmentToPhylogeny(aln) # phylogeny is a tree object.

    ##This line is just for testing purposes.
    phy = Phylo.read("temp/RAxML_bestTree.test", "newick")
    ## 3) from the Phylogeny,
    gtree = GenotypeTree(phy,conf)
    ## count independent mutations in gtree.
    #snp_total = BasicSNPCount(conf) ## just return the number of dN in genes.
    ## calculate statistics of parallel evolution.
    ## This function need to be rewritten to use GenotypeTree.
    #Statisticulate(conf, snp_total,star=False,reps=1000)

    #cleanup() # delete temp folder, doesn't work right now.

# NOTE: Updated to standards for refactor
def dNdS(conf): #old proc4
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
    print '\n'

    dNdS_counts,dNtotal,dStotal,dNdS1,dNdS2,dNdS3plus = calculate_dNdS(genomediffs)
    #print dNdS_counts
    print "dN:",dNtotal, "  dS:", dStotal, "  dN/dS:", float(dNtotal)/float(dStotal)
    print "dN/dS 1:", dNdS1, '\n', "dN/dS 2:", dNdS2, '\n', "dN/dS 3+:", dNdS3plus

# TODO: Update for refactor
def analyticalEJB(conf): #old proc5
    '''analytical solution'''

    record = utils.parse_genbank(conf.REF_GENOME)
    utils.print_genbank_summary(record)

    genomediffs = {}
    for gd_file in conf.GENOMEDIFF_FILES:
        parse_genomediff(gd_file, record, genomediffs = genomediffs)
    print '\n'

    snpcounting = snpcount(genomediffs, conf.GENOMEDIFF_FILES, conf.snp_types)

    '''
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
    '''

# NOTE: Updated for refactor
def mutationTally(conf, args, out_fn=None): #old proc6
    '''number of lines mutating in this particular gene'''
    record = utils.parse_genbank(conf.REF_GENOME)
    utils.print_genbank_summary(record)

    genomediffs = {}
    for gd_file in conf.GENOMEDIFF_FILES:
        parse_genomediff(gd_file, record, genomediffs=genomediffs)    
    print '\n'
  
    counts = mutated_lines_per_gene(genomediffs, conf.snp_types)

    if out_fn is None:
        #timestr=time.strftime('%Y_%m_%d')
        #out_fn = args.fname+timestr+'.tsv'
        out_fn = "Procedure6Output"+'.tsv'
    with open(out_fn, 'wb') as fp:
        for tag, data in counts:
            line = str(tag)
            line += '\t' + str(tuple(data['genes']))
            line += '\t' + str(len(data['lines']))
            line += '\t' + str(tuple(data['lines']))
            fp.write(line + '\n')
        #TODO: also include freq of mut types per mutating locus tag
            #e.g., 3 non-synon, 1 synon.

    #genecoords, total_bases = get_genecoordinates(record)   
    #why is this commented out? cuz I'm not here yet?
    #write_proc6_locus_mut_counts(genefreqs)

def infoRegion(conf, window_len=200): #old proc7
    '''Procedure 7: find most informative regions of the genome.
    take the union of all genome diffs; need position and originating line info.
    find windows that are most dense with SNPs for freq-seq.
    windows must: contain haplotypes that distinguish all (or many) LTEE pops.
    '''    
    ref_record = utils.parse_genbank(conf.REF_GENOME)
    mut_list = []
    conf.GENOMEDIFF_FILES.sort() # sort to ensure order is always the same
    for gd_file in conf.GENOMEDIFF_FILES:
        gd_dict = parse_genomediff(gd_file, ref_record)
        mut_list = mut_list + gd_dict.values()
    mut_list.sort(key=lambda x: x.position)
    dist_list = [y.position-x.position for x,y in window_iterator(mut_list)]
    ## add the last elt of dist_list (because chromosome is circular)
    dist_list.append(mut_list[0].position + len(ref_record.seq) - mut_list[-1].position)
    assert(sum(dist_list) == len(ref_record.seq))
    windows = []
    for i,mut in enumerate(mut_list):
        windows.append(makeWindow(i, mut_list, dist_list,window_len))
    ## filter out structural mutations.
    windows2 = []
    for w in windows:
        kinds = set([x.mut_type for x in w])
        if kinds == set(['SNP']):
            windows2.append(w)
    windows2.sort(key=lambda x: len(x),reverse=True)

    ## pick non-overlapping windows that cover all genomes as markers.
    unmarked_genomes = {k:1 for k in conf.GENOMEDIFF_FILES} # 1 if unmarked.
    markers = []
    while sum(unmarked_genomes.values()):
        best_window = []
        most_new_marks = 0
        for w in windows2:
            ## make sure mutations UNIQUELY IDENTIFY genome.
            ## CANNOT be the same new_base and position.
            unique_positions = [mut.position for mut in w]
            id_check = [(mut.position, mut.new_base) for mut in w]
            unique_muts = [mut for mut in w if id_check.count((mut.position, mut.new_base)) == 1]
            marks = [mut.fname for mut in unique_muts]
            new_marks = len([x for x in marks if unmarked_genomes[x] == 1])
            if new_marks > most_new_marks:
                most_new_marks = new_marks
                best_window = w
        markers.append(best_window)
        for mut in best_window:
            unmarked_genomes[mut.fname] = 0 # the genome is not unmarked anymore!
    ## print the ranges that are most informative.
    for i,m in enumerate(markers):
        print "MARKER", i+1, ":" 
        for mut in m:
            print mut.fname, mut.position, mut.old_base, mut.new_base, mut.gene_name

config_keys = { 'REF_GENOME': str,
                'GENOMEDIFF_FILES': str,
                'NONSYNONYMOUS': bool,
                'SYNONYMOUS': bool,
                'NONCODING': bool,
                'PSEUDOGENE': bool,
                'INTERGENIC': bool,
                'REPLICATES': int,
                'GENES_TO_DISPLAY': int }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', dest='config', 
        help="parallyze config file", default='parallyze.conf')
    parser.add_argument('--procedure', dest='procedure', 
                        default='mutationTally', choices=['simulationRM',
                        'simulationEJB','dNdS','mutationTally','infoRegion',
                        'analyticalEJB'])
    parser.add_argument('--nonsynonymous', action='store_true')
    parser.add_argument('--synonymous', action='store_true')
    parser.add_argument('--noncoding', action='store_true')
    parser.add_argument('--pseudogene', action='store_true')
    parser.add_argument('--intergenic', action='store_true')
    parser.add_argument('--replicates', type=int)
    parser.add_argument('--genes_to_display', type=int)
    parser.add_argument('--fname')
    args = parser.parse_args()
    #print args.fname

    conf = SimpleConfig(config_keys, args.config)
    for fn in conf.GENOMEDIFF_FILES:
        assert os.path.isfile(fn)

    '''
    The user can override most of the options in config file
    with (optional) flags at runtime.
    '''

    if args.nonsynonymous:
        conf.NONSYNONYMOUS = True
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
    if conf.NONSYNONYMOUS:
        conf.snp_types.add('nonsynonymous')
    if conf.SYNONYMOUS:
        conf.snp_types.add('synonymous')
    if conf.NONCODING:
        conf.snp_types.add('noncoding')
    if conf.PSEUDOGENE:
        conf.snp_types.add('pseudogene')
    if conf.INTERGENIC:
        conf.snp_types.add('intergenic')

    if args.procedure == 'simulationRM':
	    simulationRM(conf)
    elif args.procedure == 'simulationEJB':
	    simulationEJB(conf)
    elif args.procedure == 3:
	    proc3(conf)
    elif args.procedure == 'dNdS':
	    dNdS(conf)
    elif args.procedure == 'analyticalEJB':
        analyticalEJB(conf)
    elif args.procedure == 'mutationTally':
        mutationTally(conf,args)
    elif args.procedure == 'infoRegion':
        infoRegion(conf)

main()

