import numpy as np
import operator

import utils

def calculate_dNdS(genomediffs):
    '''input: the output from parsed gd files
    goal: count the SNP mutation types for dN/dS, intergenic, etc
    output: dict of dict of a count'''

    counts = {}

    for mut_id, gd in genomediffs.iteritems():
        if gd.mut_type == 'SNP' and \
            gd.snp_type in ['synonymous', 'nonsynonymous']:
            
            # These all only have one locus_tag, and we can't
            # use a list as a key, so just get the value
            locus_tag = gd.locus_tag[0]

            if locus_tag in counts:
                dN, dS = counts[locus_tag]
            else:
                dN, dS = (0, 0)

            if gd.snp_type == 'nonsynonymous':
                dN += 1
            else:
                dS += 1
            
            counts[locus_tag] = (dN, dS)

    return counts

# TODO: Update for refactor
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

# NOTE: Updated for refactor
def mutated_lines_per_gene(genomediffs, snp_types):
    '''input: mutations, from parse_gdfiles
    given user input #x, list # mutated lines/gene'''

    '''
    Store the number of lines mutated for each gene,
    keyed by locus_tag
    '''
    tag_lines = {}

    for key, gd in genomediffs.iteritems():
        if gd.mut_type == 'SNP' and gd.snp_type in snp_types:
            for tag in gd.locus_tag:
                if tag not in tag_lines:
                    tag_lines[tag] = set()
                tag_lines[tag].add(gd.line)

    sorted_tag_lines = zip(tag_lines.keys(), tag_lines.values())
    sorted_tag_lines = sorted(sorted_tag_lines, key=lambda x: len(x[1]), reverse = True)
    
    return sorted_tag_lines

# TODO: Update for refactor
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

# TODO: Refactor
def get_mut_sites(matrices, refseq, num_replicates):
    mut_sites = {}
    refseq_arr = np.array([c for c in refseq]) 
    for filename in matrices:
        print '\n** generating mutation sites for', filename
        mut_sites[filename] = snpmutate(matrices[filename], num_replicates, refseq_arr)
    print
    '''
    mut_sites = { 'filename1': { 'A': array(row for reps, columns for mut position indices), 
                                 'T': array(...) ... },
                  'filename2': {...},
                  ...
                }
    '''
    return mut_sites

