#!/usr/bin/python

#usage: python parallyze.py config.txt [genomediff1.gd genomediff2.gd ... n : referencegenome.gb(k) {genbank}]

#find config.txt
#find .gd files, save in list
#try to break

import argparse
import config
import os
import sys

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    conf = {}
    if PROCEDURE == 1:
        reference = REFERENCE_GENOME.strip()
        assert os.path.isfile(reference)
        assert reference.endswith('.gb')
        conf['procedure'] = PROCEDURE
        conf['reference'] = reference
        return conf
    elif PROCEDURE == 2:
        reference = REFERENCE_GENOME.strip()
        assert os.path.isfile(reference)
        assert reference.endswith('.gb')
        diffs = file_list(GENOME_DIFFS)
        for diff_file in diffs:
            assert diff_file.endswith('.gd')
        conf['reference'] = reference
        conf['procedure'] = PROCEDURE
        conf['diffs'] = diffs
        return conf
    else:
        print >>sys.stderr, 'Invalid procedure {p}! Exiting...'.format(p=PROCEDURE)
        sys.exit()

def main():
	parser = argparse.ArgumentParser()
#	parser.parse_args()
	parser.add_argument(dest='config', help="First entry should be a config.py file. File should be in working folder.")
	parser.add_argument('-r', dest='reference', default="NONE")
    args = parse.parse_args()

    conf = get_config()
    print conf['procedure'], conf['reference']


main()
