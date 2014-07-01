#!/usr/bin/env python  
#usage: python parallyze.py
'''
Things to do:
1. change genefreqs and genecoords to use the locus tag because many gene names are 'none'
2. add asserts for list sizes to make sure results are sane
3. change locations of mutations to be dropped to align with included/specified gd mutations (ie, synon, non-synon, non-coding, etc)
4. change from #muts/gene to #lines/gene
5. analytic solution for SNPs
6. dN/dS ratio
7. get genes that positions fall into. have # of times hit, not just hit/nothit (best method for this? locus tag? gene? (do intergenic now - intragenic later)
8. sum, across all lines, # of mutations per gene, for experimental data (done) and simulated data - divide by reps for avg. 
9. print bar graphs? 
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
#from Bio.Seq import MutableSeq

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

'''
def base_to_int(i):
    if i == 'A':
        return 0
    elif i == 'G':
        return 1
    elif i == 'C':
        return 2
    else:
        return 3

def int_to_base(i):
    if i == 0:
        return 'A'
    elif i == 1:
        return 'G'
    elif i == 2:
        return 'C'
    else:
        return 'T'

def seq_to_int(seq):
    converted_to_int=[base_to_int(b) for b in seq]
    return converted_to_int

def int_to_seq(seq):
    converted_to_seq=[int_to_base(b) for b in seq]
    return converted_to_seq
'''

def complementary_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

def make_record(ref_file):
    '''input:conf['ref'] 
    returns ref genome as a string 'refseq' '''

    '''
    SeqIO.parse is an iterator and so has a method named "next"
    which is called when you use it in a for loop, ie:
    for record in SeqIO.parse(stuff):
        do stuff
    We can call it explicitly, once, because we're assuming there is only
    one record in the SeqIO iterator. We will be *very* explicit and store
    the iterator itself as it, then  call next() on it like so:
    '''
    conf = get_config()
    it = SeqIO.parse(ref_file, "genbank")
    record = it.next()
    print '\n', record
    return record  

def get_refseq(record):
    # convert the biopython Seq object to a python string
    refseq = list(record.seq)
    length=len(refseq)    
    countA=refseq.count('A'); countG=refseq.count('G')
    countC=refseq.count('C'); countT=refseq.count('T')
    print '\n', 'Base distribution in reference: '
    print 'A:', countA, '   G:', countG, '   C:', countC, '   T:', countT
    # or seq = ''.join(seq) to save as string & overwrite the list
    ##print 'Seq as string [truncated]:', ''.join(seq)[:1000], '...'
    print "Number of bases: ", length
    print 'Seq as condensed string:', ''.join(refseq)[0:100], '...', ''.join(refseq)[length-100:length] 
    print
    return refseq

def get_genecoordinates(record):
    '''input: record file from reference (i.e., ref as record file)
    get all IDed genes and their beginning and end position;
    store in list of tuples'''
    geneinfo = []
    for i in record.features:
        if i.type == 'CDS': #or i[type] if dictionary - tab everything below
        #add non-CDS records later? (tRNA, lonely 'gene', repeat_region, misc_feature
            my_start = int(i.location.start)
            my_end = int(i.location.end)
            locus_tag = i.qualifiers['locus_tag'][0]
            try:
                 gene_name = i.qualifiers['gene'][0]
            except KeyError: 
                gene_name = 'none'
            mytuple = (my_start, my_end, gene_name, locus_tag)
            geneinfo.append(mytuple)
    print '\n', 'Reference genome gene list (1st 10): start position, end pos, gene, locus_tag', '\n', geneinfo[:10]
    return geneinfo

def str_keyvalue(data):
    s = '\n'.join([str(key)+': '+str(data[key]) for key in data])
    return s

#is this appropriate refgenome? there are multiple - maybe? REL606.6 
#need updated breseq?

def parse_gdfiles(filenames, refseq):
    '''input: conf['diffs'] and conf['ref']
    :returns dict of gdfile names, each of which contains a list of sequential #s;
    each # corresponds to a dictionary for each mutation (contains key and value)''' 
    alldiffs={}
    #disregarded_evidence = ['#', 'JC', 'RA', 'UN', 'MC', 'NOTE']
    for fname in filenames:
        alldiffs[fname] = []
        with open(fname) as fp:
            for line in fp: 
                #if line.startswith(i) for i in disregard_evidence:
                if line.startswith('#') or \
                    line.startswith('JC') or \
                    line.startswith('RA') or \
                    line.startswith('UN') or \
                    line.startswith('MC') or \
                    line.startswith('NOTE'):
                    continue
                line = line.split('\t')
                data = {}
                mut_type = line[0]
                mut_id = line[1]
                data['mut_type'] = line[0]
                data['mut_id'] = line[1]
                data['parent_ids'] = line[2].split(',')
                data['seq_id'] = line[3]
                data['position'] = int(line[4])-1
                if mut_type == 'SNP':
                    data['new_base'] = line[5]
                    for pair in line[6:]:
                        key,_,value = pair.partition('=')
                        data[key.strip()] = value.strip()
                    data['gene_name'] = [gene.strip() for gene in data['gene_name'].split('/')]
                    data['locus_tag'] = [tag.strip() for tag in data['locus_tag'].split('/')]

                #        gene_names = data['gene_name']
                #        data['gene_name'] = []
                #        for name in gene_names.split('/'):
                #            data['gene_name'].append(name.strip())
                #        data['gene_name'] =[value.strip()]
                #    if data['snp_type'] == 'intergenic':
                #        data[value] = value.split('/')
                    if data['snp_type'] in ['nonsynonymous', 'synonymous'] and data['gene_strand'] == '<':
                        data['new_base'] = complementary_base(line[5])
                    if data['snp_type'] in ['nonsynonymous', 'synonymous']:
                        data['codon_position'] = int(data['codon_position']) - 1
                        data['old_base'] = data['codon_ref_seq'][data['codon_position']]
#                        if data['old_base'] == data['new_base']:
#                            print 'WARNING: new base same as old base', ' ', data['snp_type'], ' inconsistencey in gdfile ', fname
#                            print str_keyvalue(data)
#                            print 'refseq old_base:', refseq[data['position']]
#                            start = data['position']-2
#                            end = start+4
#                            print 'refseq sequence:', '({}:{})'.format(start,end-1), ''.join(refseq[start:end]), '\n'      
                    elif data['snp_type'] in ['intergenic', 'pseudogene', 'noncoding']:
                        #print data['position']
                        data['old_base'] = refseq[data['position']]
#                        if data['snp_type'] == 'noncoding':
#                            print '^^'*30, fname, str_keyvalue(data), '^^'*30
                        if data['new_base'] == data['old_base']:
                            print 'WARNING: new base same as old base', ' ', data['snp_type'], '  problem with refgenome or indexing'
                    else:
                        print 'ERROR parsing SNP ', data['snp_type']
                elif mut_type == 'SUB':
                    data['size'] = int(line[5])
                    data['new_seq'] = line[6]
                elif mut_type == 'DEL':
                    data['size'] = int(line[5])
                elif mut_type=='INS':
                    data['new_seq'] = line[5]
                elif mut_type == 'MOB':
                    data['repeat_name'] = line[5]
                    data['strand'] = line[6]
                    data['duplication_size'] = int(line[7])
                elif mut_type == 'AMP':
                    data['size'] = int(line[5])
                    data['new_copy_number'] = int(line[6])
                elif mut_type == 'CON':
                    data['size'] = int(line[5])
                    data['region'] = line[6]
                elif mut_type == 'INV':
                    data['size'] = int(line[5])
                alldiffs[fname].append(data)
    #for gdname, mutation_list in alldiffs.iteritems():
        #print "File:", gdname
        #for list_position, mutation in enumerate(mutation_list):
            #if list_position < 1:
                #print "Mutation", list_position, "\n", str_keyvalue(mutation), "\n", "*" * 40
        ###print dbug_dict[f_mutid], '\n'
        ###dbug_dict = alldiffs[alldiffs.keys()[0]] 
    return alldiffs

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

def gene_rank_and_mutate_parameters():
    rank_and_mut_params = {}
    rank_and_mut_params['synonymous'] = input("Include synonymous mutations in the analysis? 0 for N, 1 for Y: ")
    rank_and_mut_params['noncoding'] = input("Include noncoding SNPs in the analysis? 0 or 1: ")
    rank_and_mut_params['pseudogene'] = input("Include SNPs in psuedogenes in the analysis? 0 or 1: ")
    rank_and_mut_params['intergenic'] = input("Include intergenic SNPs in the analysis? 0 or 1: ")
    rank_and_mut_params['number_of_top_genes'] = int(input("How many of the most frequently mutated genes would you like displayed? "))
#    rank_and_mut_params['generations'] = input("For how many generations did your experiment run? ")
    rank_and_mut_params['replicates'] = input("How many replicates would you like? ")
    print
    return rank_and_mut_params

def lines_gene_rank(filenames, params):
    '''input: mutations, from parse_gdfiles
    given user input #x, list # mutated lines/gene'''
    mut_genes = {}
    rank_mut_types = ['nonsynonymous']
    if params['synonymous'] == 1:
        rank_mut_types.append('synonymous')
    if params['noncoding'] == 1:
        rank_mut_types.append('noncoding')
    if params['pseudogene'] == 1:
        rank_mut_types.append('pseudogene')
    if params['intergenic'] == 1:
        rank_mut_types.append('intergenic')
    for fname in filenames: #if a file/line has any number of a particular mut, add 1
                            #so max number for a gene is the number of files
        for mut in filenames[fname]:
            if mut['mut_type'] == 'SNP' and mut['snp_type'] in rank_mut_types:
                for tag in mut['locus_tag']:
                    if tag not in mut_genes:
                        mut_genes[tag] = set()
                    mut_genes[tag].add(fname)
    srt_mut_genes = zip(mut_genes.keys(), mut_genes.values())
    srt_mut_genes = sorted(srt_mut_genes, key=lambda x: len(x[1]))
    for row in srt_mut_genes:
        print row[0], len(row[1]), row[1]
    return srt_mut_genes

'''
def gds_gene_rank(filenames, params):
    #input: mutations, from parse_gdfiles
    #given user input #x, list x most mutated genes across all gdfiles
    #(#muts/gene across all lines)
    mut_genes = {}
    rank_mut_types = ['nonsynonymous']
    if params['synonymous'] == 1:
        rank_mut_types.append('synonymous')
    if params['noncoding'] == 1:
        rank_mut_types.append('noncoding')
    if params['pseudogene'] == 1:
        rank_mut_types.append('pseudogene')
    if params['intergenic'] == 1:
        rank_mut_types.append('intergenic')
    for fname in filenames: 
        for mut in filenames[fname]:
            if mut['mut_type'] == 'SNP' and mut['snp_type'] in rank_mut_types:
                for tag in mut['locus_tag']:
                    mut_genes[tag] = mut_genes.get(tag, 0) + 1
    print mut_genes
    return mut_genes
    #sorted_mut_genes = sorted(mut_genes.iteritems(), key=operator.itemgetter(1), reverse = True)
    #mut_genes_number = int(len(sorted_mut_genes))
    #print 'The', params['number_of_top_genes'], 'most mutated genes of all', mut_genes_number, 'mutated genes:'
    #print sorted_mut_genes[:params['number_of_top_genes']]
    #return sorted_mut_genes
    ##diff btwn intergenic and noncoding (has no new base)? pseudogene? all exclusive?
'''

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
        for start, end, name, tag in genecoords:
            row = [tag]
            for filename in mut_sites:
                line_muts = 0
                for origbase in mut_sites[filename]:
                    line_muts += ((mut_sites[filename][origbase] >= start) & (mut_sites[filename][origbase] < end)).sum()
                row.append(line_muts)
            outfp.write(', '.join([str(c) for c in row]) + '\n')

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

def dnds_calculate(diff_dict):
    '''input: the output from parsed gd files
    goal: count the SNP mutation types for dN/dS, intergenic, etc
    output: dict of dict of a count'''
    muttypes = {}
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

def proc1(conf):
    '''simulated solution'''
    params = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    `matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords = get_genecoordinates(record)
    with open('genecoords.txt', 'wb') as fp:
        fp.write(str(genecoords))
    mut_sites = get_mut_sites(matrices, refseq, params['replicates'])
    write_gene_mut_counts(genecoords, mut_sites)
    write_gd_gene_mut_counts(genecoords, genefreqs)

'''old proc3 assumptions: Synonymous mutations are neutral, Mutations are independent of one another, No defects to DNA repair, Mutation rate is constant across the genome, There is only one chromosome'''

def proc4(conf):
    '''dnds calculation'''
    params = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    dnds_calculate(mutations)

def proc5(conf):
    '''analytical solution'''
    params = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords = get_genecoordinates(record)
    with open('genecoords.txt', 'wb') as fp:
        fp.write(str(genecoords))
    mut_sites = get_mut_sites(matrices, refseq, params['replicates'])
    write_gene_mut_counts(genecoords, mut_sites)
    write_gd_gene_mut_counts(genecoords, genefreqs)

def proc6(conf):
    '''Mike's: number of lines mutating in this particular gene'''
    params = gene_rank_and_mutate_parameters()
    record = make_record(conf['ref'])
    refseq = get_refseq(record)
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrices = snpcount(mutations)
    genefreqs = lines_gene_rank(mutations, params)
    #genefreqs = {key:value for key,value in genefreqs}
    genecoords = get_genecoordinates(record)
    with open('genecoords.txt', 'wb') as fp:
        fp.write(str(genecoords))

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

