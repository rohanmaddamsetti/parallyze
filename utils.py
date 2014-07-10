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
    # TODO: Look up how to quiet SeqIO.parse
    record = SeqIO.parse(gb_file, 'genbank').next()
    return record

def print_genbank_summary(gb_record):
    '''
    Print a summary of the given genbank record
    '''
    print '\n[Genbank record object]'
    print '#', gb_record.description
    print '#', gb_record.id, gb_record.name
    if 'date' in gb_record.annotations:
        print '#', gb_record.annotations['date']
    if 'sequence_version' in gb_record.annotations:
        print '# version', gb_record.annotations['sequence_version']
    if 'taxonomy' in gb_record.annotations:
        print '#', '; '.join([t for t in gb_record.annotations['taxonomy']])
    print '#\n#\tLength:', len(gb_record.seq)
    print '#\tA:', gb_record.seq.count('A'), '\tT:', gb_record.seq.count('T')
    print '#\tG:', gb_record.seq.count('G'), '\tC:', gb_record.seq.count('C')
    gc = gb_record.seq.count('G') + gb_record.seq.count('C')
    print '#\t%GC:', 100.0 * gc / len(gb_record.seq), '\n'
