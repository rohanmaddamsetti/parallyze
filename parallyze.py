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
    conf={}
    if config.PROCEDURE == '3':
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb' or '.gbk') 
        conf['procedure'] = config.PROCEDURE
	conf['ref'] = ref
	assert len(config.GENOME_DIFFS.strip())==0
        return conf
    elif config.PROCEDURE in ['1','2','4','5']:  
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb' or '.gbk')
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

get_config()

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

def parse_gdfiles(filenames): #change filenames - conf['diffs']
    mutations={}
    for fname in filenames:
        mutations[fname] = {}
        with open(fname) as fp:
            for line in fp: 
                if line.startswith('#') or \
                    line.startswith('JC') or \
                    line.startswith('RA') or \
                    line.startswith('UN'):
                    continue
                line = line.split()
                mut_type = line[0]
                mut_id = line[1]
                parent_ids = line[2].split(',')
                seq_id = line[3]
                position = line[4]
                data = {}
                data['mut_type'] = mut_type
                data['mut_id'] = mut_id
                data['parent_ids'] = parent_ids
                data['seq_id'] = seq_id
                data['position'] = position
                if mut_type == 'SNP':
                    data['new_se'] = line[5]
                    key,_,value = line[12].partition('=')
                    if key == 'codon_ref_seq': #intragenic SNPs
                        data[key] = value
                        key,_,value = line[11].partition('=')
                        data['codon_position'] = int(value)-1
                        key,_,value = line[15].partition('=')
                        data['gene_name'] = value
                        key,_,value = line[17].partition('=')
                        data['gene_product'] = value
                    elif key == 'snp_type': #intergenic SNPs
                        data[key] = value
                        key,_,value = line[8].partition['=']
                        data['gene_position'] = value
                        key,_,value = line[9].partition['=']
                        data['gene_product'] = value
                elif mut_type == 'SUB':
                    data['size'] = line[5]
                    data['new_se'] = line[6]
                elif mut_type == 'DEL':
                    data['size'] = line[5]
                elif mut_type=='INS':
                    data['new_seq'] = line[5]
                elif mut_type == 'MOB':
                    data['repeat_nam'] = line[5]
                    data['strand'] = line[6]
                    data['duplication_size'] = line[7]
                elif mut_type == 'AMP':
                    data['size'] = line[5]
                    data['new_copy_number'] = line[6]
                elif mut_type == 'CON':
                    data['size'] = line[5]
                    data['region'] = line[6]
                elif mut_type == 'INV':
                    data['size'] = line[5]
                mutations[fname][mut_id] = data
    dbug_dict = mutations[mutations.keys()[0]] 
    for f_mutid in dbug_dict.keys()[:2]:
        print dbug_dict[f_mutid]
    return mutations

def snpcount(filename): #conf['ref']
    mutmatrix = {} #should be array. append dicts of dicts to the array
    for fname in mutations:
        init_base={}
        init_base={'A':'{'G':0, 'C':0, 'T':0}, 'G':{}, 'C':{}, 'T':{}}
        for mutdict in fname:
            if data['mut_type']=='SNP':
                old_base = conf['ref'][data['position']]
                new_base = data['new_se']
                d[old_base][new_base] = d[old_base].get(new_base,0) + 1
    return mutmatrix #this doesn't add anything to mutmatrix

def snpmutate (filename): #conf['ref']  #throw error if non-annotated genomediff?
    pass

def gds_gene_rank():
    pass
    mut_genes = {}
    for fname in filenames:
        if data['mut_type'] == 'SNP':
            gene_name = data['gene_name']
            d[gene_name] = d[gene_name].get() + 1

def proc1(conf):
    parse_gdfiles
    snpcount

def proc3(conf):
    print '\n', 'Assumptions:', '\n', 'Synonymous mutations are neutral' '\n', 'Infinite sites model', '\n', 'Mutations are independent of one another', '\n', 'No defects to DNA repair', '\n', 'Mutation rate is constant across the genome', '\n', 'There is only one chromosome', '\n'
    lines=input("How many lines?  ")
    gens=input("How many generations? ")
    reps=input("How many replicates?  ")

    '''
    SeqIO.parse is an iterator and so has a method named "next"
    which is called when you use it in a for loop, ie:
    for record in SeqIO.parse(stuff):
        do stuff
    We can call it explicitly, once, because we're assuming there is only
    one record in the SeqIO iterator. We will be *very* explicit and store
    the iterator itself as it, then  call next() on it like so:
    '''
    it = SeqIO.parse(conf['ref'], "genbank")
    record = it.next() 

    # convert the biopython Seq object to a python string
    seq = list(str(record.seq))
    length=len(seq)
    ## print 'Seq as list [truncated]:', seq[:5000], '...'
    
    countA=seq.count('A') #possibly change to 0,1,2,3
    countG=seq.count('G')
    countC=seq.count('C')
    countT=seq.count('T')
    print 'A:', countA, '   G:', countG, '   C:', countC, '   T:', countT

    # or do seq = ''.join(seq) to save it as a string and overwrite the list
    ##print 'Seq as string [truncated]:', ''.join(seq)[:1000], '...'
    print "Number of bases: ", length
    print 'Seq as condensed string:', ''.join(seq)[0:100], '...', ''.join(seq)[length-100:length] 

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
        print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref'], '\n','Genome diffs: ', conf['diffs']
	if conf['procedure']=='1':
	    proc1(conf)
	elif conf['procedure']=='2':
	    proc2(conf)
	elif conf['procedure']=='4':
	    proc4(conf)
	elif conf['procedure']=='5':
	    proc5(conf)
main()
