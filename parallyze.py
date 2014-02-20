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

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    conf = {}
    if config.PROCEDURE == 3:
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb') 
        conf['procedure'] = config.PROCEDURE
        conf['ref'] = config.ref
        return conf
    elif config.PROCEDURE in [1,2,4,5]:  
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb')
        diffs = file_list(config.GENOME_DIFFS)
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
    #parser.add_argument(dest='config', help="First entry should be a config.py file. File should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    print >>sys.stderr, 'configuration: ', conf['procedure'], conf['ref'], conf['diffs']

main()
