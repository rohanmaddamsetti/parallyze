#!/usr/bin/env python

#usage: python parallyze.py config.txt [genomediff1.gd genomediff2.gd ... n : referencegenome.gb(k) {genbank}]

import argparse
import os
import sys

import config
import config_default

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
        assert ref.endswith('.gb') 
        conf['procedure'] = config.PROCEDURE
	conf['ref'] = ref
	assert len(config.GENOME_DIFFS.strip())==0
        return conf
    elif config.PROCEDURE in ['1','2','4','5']:  
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb')
        diffs = file_list(config.GENOME_DIFFS)
	assert len(diffs)!=0
        for diff_file in diffs:   #annotated vs. non-annotated genomediff files
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

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(dest='config', help="Config.py file should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    if conf['procedure']=='3':
	print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref']
    else: 
        print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref'], '\n','Genome diffs: ', conf['diffs']
main()

def defreflist(filename):
    reflist=[]
#ref=''
#fp=open(filename, 'rU')
#for record in SeqIO.parse(filepath, "genbank"):
    #ref += repr(record.seq)
        #if I wanted all genome parts concatenated into one string
#determine length of ref genome. assign a position number to each base.
#create list or array (?) of the ref.gen.
#	for line in file('reference.gb,' 'r'):
#	    if line.startswith(a number) and line.isupper():
		#peel off first number, assign each character to a sequential number
		#array.count(x) or len(<seq>)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#for seq_record in SeqIO.parse("reference.gb", "genbank"):
    #print(seq_record.id)
    #refseq=repr(seq_record.seq)
    #print(len(seq_record))
    #concatenate - or not - string.join
    #do in fasta. 
    #what's the output format? - did rohan ever do a 2nd mtg w/ tracy?
    #mutable_refseq=refseq.tomutable()
def proc3(conf):
    if conf['procedure']=='3': #more elegant way to do this?
	print '\n', 'Assumptions:', '\n', 'Synonymous mutations are neutral' '\n', 'Infinite sites model', '\n', 'Mutations are independent of one another', '\n', 'No defects to DNA repair', '\n', 'Mutation rate is constant across the genome', '\n', 'There is only one chromosome', '\n'
    #if conf['procedure']=='3': #more elegant way to do this?
proc3()
