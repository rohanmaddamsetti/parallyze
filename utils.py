'''

Part of parallyze
by Elizabeth Baird
+ Camille Scott

Contains parallyze general utility functions.

'''

from Bio import SeqIO

dna = ['A', 'T', 'C', 'G']

def complementary_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

def str_keyvalue(data):
    s = '\n'.join([str(key)+': '+str(data[key]) for key in data])
    return s

def parse_genbank(gb_file):
    '''
    Parse a single-locus genbank file with SeqIO, advance the iterator
    to the first item and return that record.
    '''
    record = SeqIO.parse(gb_file, 'genbank')
    return record

def print_genbank_summary(gb_record):
    '''
    Print a summary of the given genbank record
    '''
    print '[Genbank record object]'
    print '#', record.description
    print '#', record.id, record.name
    if 'date' in record.annotations:
        print '#', record.annotations['date']
    if 'sequence_version' in record.annotations:
        print '# version', record.annotations['sequence_version']
    if 'taxonomy' in record.annotations:
        print '#', '; '.join([t for t in record.annotations.taxonomy])
    print '#\n#\tLength:', len(record.seq)
    print '#\tA:', record.seq.count('A'), '\tT:', record.seq.count('T')
    print '#\tG:', record.seq.count('G'), '\tC:', record.seq.count('C')
    gc = record.seq.count('G') + record.seq.count('C')
    print '#\t%GC:', 100.0 * gc / len(record.seq)
