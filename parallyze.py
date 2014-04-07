#!/usr/bin/env python   #is there a reason this was deleted before?
#usage: python parallyze.py

import argparse
import os
import sys

import config
import config_default

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

from numpy import random

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    '''returns object that holds procedure name, ref file name, and diff file names'''
    conf={}
    if config.PROCEDURE == '3':
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb') or ref.endswith('.gbk') 
        conf['procedure'] = config.PROCEDURE
	conf['ref'] = ref
	assert len(config.GENOME_DIFFS.strip())==0
        return conf
    elif config.PROCEDURE in ['1','2','4','5']:  
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

def complementary_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

def parse_ref(ref_file):
    '''input:conf['ref'] 
    #returns ref genome as a string 'refseq' '''

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

    # convert the biopython Seq object to a python string
    # refseq = list(str(record.seq))
    refseq = list(record.seq)
    length=len(refseq)    

    countA=refseq.count('A') #possibly change to 0,1,2,3
    countG=refseq.count('G')
    countC=refseq.count('C')
    countT=refseq.count('T')
    print '\n', 'Base distribution in reference: '
    print 'A:', countA, '   G:', countG, '   C:', countC, '   T:', countT

    # or do seq = ''.join(seq) to save it as a string and overwrite the list
    ##print 'Seq as string [truncated]:', ''.join(seq)[:1000], '...'
    print "Number of bases: ", length
    print 'Seq as condensed string:', ''.join(refseq)[0:100], '...', ''.join(refseq)[length-100:length] 
    print

    return refseq

def str_keyvalue(data):
    s = '\n'.join([str(key)+': '+str(data[key]) for key in data])
    return s

#IS THIS THE APPROPRIATE REFGENOME? THERE ARE MULTIPLE - maybe? REL606.6 
#NEED UPDATED BRESEQ?

def parse_gdfiles(filenames, refseq):
    '''input: conf['diffs'] and conf['ref']
    :returns dict of gdfile names, each of which contains a list of sequential #s;
    each # corresponds to a dictionary for each mutation (contains key and value)''' 
    alldiffs={}
    for fname in filenames:
        alldiffs[fname] = []
        with open(fname) as fp:
            for line in fp: 
                if line.startswith('#') or \
                    line.startswith('JC') or \
                    line.startswith('RA') or \
                    line.startswith('UN'):
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
#                        data['gene_name'] =[value.strip()]
#                    if data['snp_type'] == 'intergenic':
#                        split/partition by '/'  
#                                    '''
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
    '''input: 'mutations'  #still confused about this. what is this? what should go here?
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
    rank_and_mut_params['generations'] = input("For how many generations did your experiment run? ")
    rank_and_mut_params['replicates'] = input("How many replicates would you like? ")
    print
    return rank_and_mut_params

def gds_gene_rank(filenames, params):
    '''input: mutations, from parse_gdfiles
    given user input #x, list x most mutated genes across all gdfiles'''
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
                gene_name = mut['gene_name']
                mut_genes[gene_name] = mut_genes.get(gene_name, 0) + 1
    import operator #better way to pull top mutators than sorting? 
    sorted_mut_genes = sorted(mut_genes.iteritems(), key=operator.itemgetter(1))
    mut_genes_number = int(len(sorted_mut_genes))
    most_mutated_genes = {}  #store top mutators
    print 'The', params['number_of_top_genes'], 'most mutated genes of all', mut_genes_number, 'mutated genes:'
    print sorted_mut_genes[mut_genes_number - params['number_of_top_genes']:mut_genes_number] #put most_mutated_genes in here instead once it's working
    ##diff btwn intergenic and noncoding (has no new base)? pseudogene? all exclusive? 

def snpmutate (filename1, filename2):  #throw error if non-annotated genomediff?
    '''input: matrixdict and refseq from snpcount and parse_ref, respectively
    goal: mutate one genome once, output positions of mutations
    will later do:: for i in gdfiles: for i in # reps
    also later: # of mutations per gene in reps, etc.'''
    pass

def proc1(conf):
    params = gene_rank_and_mutate_parameters()
    refseq = parse_ref(conf['ref'])
    mutations = parse_gdfiles(conf['diffs'], refseq)
    matrix = snpcount(mutations)
    gds_gene_rank(mutations, params)
    #return mutations
    #return refseq
    #return matrix

def proc3(conf):
    print '\n', 'Assumptions:', '\n', 'Synonymous mutations are neutral'\
        '\n', 'Infinite sites model', '\n', 'Mutations are independent of one another',\
        '\n', 'No defects to DNA repair', '\n', 'Mutation rate is constant across the genome',\
        '\n', 'There is only one chromosome', '\n'
    lines=input("How many lines?  ")
    gens=input("How many generations? ")
    reps=input("How many replicates?  ")

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(dest='config', help="Config.py file should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    if conf['procedure']=='3':
        print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref']
        proc3(conf)
    else: 
        print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref'], '\n','Genome diffs: ', conf['diffs'], '\n'
	if conf['procedure']=='1':
	    proc1(conf)
	elif conf['procedure']=='2':
	    proc2(conf)
	elif conf['procedure']=='4':
	    proc4(conf)
	elif conf['procedure']=='5':
	    proc5(conf)
main()
