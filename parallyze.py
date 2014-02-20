#!/usr/bin/env python

#usage: python parallyze.py config.txt [genomediff1.gd genomediff2.gd ... n : referencegenome.gb(k) {genbank}]

#find config.txt
#find .gd files, save in list
#try to break

import argparse
import os
import sys

import config
import config_default
import config_test

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    conf = {}
    if config.PROCEDURE == 3:
        reference = config.REFERENCE_GENOME.strip()
        assert os.path.isfile(reference)
        assert reference.endswith('.gb') 
        conf['procedure'] = config.PROCEDURE
        conf['reference'] = config.reference
        return conf
    elif config.PROCEDURE in [1,2,4,5]:  
        reference = config.REFERENCE_GENOME.strip()
        assert os.path.isfile(reference)
        assert reference.endswith('.gb')
        diffs = file_list(config.GENOME_DIFF_FILES)
        for diff_file in diffs:   #annotated vs. non-annotated genomediff files
            assert diff_file.endswith('.gd')
        conf['reference'] = reference
        conf['procedure'] = config.PROCEDURE
        conf['diffs'] = diffs
        return conf
    else:
        print >>sys.stderr, 'Invalid procedure {p}! Using test data'.format(p=config.PROCEDURE)
        conf['reference'] = config_test.REFERENCE_GENOME
        conf['procedure'] = config_test.PROCEDURE
        conf['diffs'] = config_test.GENOME_DIFF_FILES
        return conf

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(dest='config', help="First entry should be a config.py file. File should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    print >>sys.stderr, 'configuration: ', conf['procedure'], conf['reference']

main()
