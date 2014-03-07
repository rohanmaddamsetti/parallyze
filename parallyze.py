#!/usr/bin/env python

#install BioPython and numpy (?) first (see README.md)
#usage: python parallyze.py

import argparse
import os
import sys

import config
import config_default

#import biopython and numpy

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

def reflist(filename):
    reflist=[]
#ref=''
#fp=open(filename, 'rU')
#determine length of ref genome. assign a position number to each base.
#create list or array (?) of the ref.gen.
#	for line in file('reference.gb,' 'r'):
#	    if line.startswith(a number) and line.isupper():
		#peel off first number, assign each character to a sequential number
		#array.count(x) or len(<seq>)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

from numpy import random
#for seq in SeqIO.parse(conf['ref'], "genbank"):
    #print(seq_record.id)
    #refseq=repr(seq_record.seq)
    #print(len(seq_record))
    #what's the output format? 
    #mutable_refseq=refseq.tomutable()

def proc3(conf):
    print '\n', 'Assumptions:', '\n', 'Synonymous mutations are neutral' '\n', 'Infinite sites model', '\n', 'Mutations are independent of one another', '\n', 'No defects to DNA repair', '\n', 'Mutation rate is constant across the genome', '\n', 'There is only one chromosome', '\n'
    replicates=input("How many replicates?  ")
    lines=input("How many lines?  ")
    for record in SeqIO.parse(conf['ref'], "genbank"):
        refseq=repr(record.seq) #does this do what i want?
        #mutrefseq=record.tomutable()
        length=len(record)
        print"Number of bases: ", length
        #print(record[0:9], '...', record[length-10:length])
#rohan: pass whole object or just relevant bits of configuration

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
