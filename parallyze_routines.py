import numpy as np
import operator
import os

from genomediff import parse_genomediff, GenomeDiff
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO


import utils

def calculate_dNdS(genomediffs):
    '''input: the output from parsed gd files
    goal: count the SNP mutation types for dN/dS, intergenic, etc.
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

# NOTE: Updated for refactor
def snpcount(genomediffs, lines, snp_types):
    '''input: 'mutations' i.e., parsed gd files
    returns dictionary of gdfiles, each containing a matrix
    of SNP mutations - to and from base'''

    file_matrices = {}
    for line in lines:
        file_matrices[line] = np.zeros((4,4), dtype = int)

    for key, gd in genomediffs.iteritems():
        if gd.mut_type=='SNP' and gd.snp_type in snp_types:
            old_base = utils.seq_to_int(gd.old_base)
            new_base = utils.seq_to_int(gd.new_base)
            line = gd.line
            file_matrices[line][old_base, new_base] += 1

    for line, mat in file_matrices.iteritems():
        print 'file:', line
        print 'from (row) / to (column) :', '\n', mat, '\n'

    return file_matrices

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

# TODO: Update for refactor? low priority
# TODO: if position is intergenic (or etc.), choose new position
def snpmutate(mat, num_replicates, refseq_arr): 
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
#below was recently commented out for purposes of functionality
#    mut_sites = {}
#    for origbase in mat:
#        num_muts = #numpy.zeros sum across matrix row

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
    '''

# TODO: Refactor? low priority. do for comparison btwn analytic and simulation
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

def SNPsToAlignment(conf):
    ''' rows are lexicographically sorted conf.GENOMEDIFF_FILES, and
    the last row is the reference sequence. columns are all positions 
    that evolved in the set of genomes.
    '''
    
    ref_record = utils.parse_genbank(conf.REF_GENOME)
    #utils.print_genbank_summary(ref_record)

    snps = []
    ## each elt in snps is a tuple: (position, old_base, new_base, locus_tag, label)
    conf.GENOMEDIFF_FILES.sort() ## So I can assume the diffs are sorted.
    for gd_file in conf.GENOMEDIFF_FILES:
        gd_dict = parse_genomediff(gd_file, ref_record)
        for k, v in gd_dict.iteritems():
            if v.mut_type != 'SNP':  ## only consider SNPs
                continue
            old_base = ref_record[v.position]
            snps.append((v.position, old_base, v.new_base, v.locus_tag[0], gd_file))
    snps.sort(key=lambda elt: elt[0]) #sort by position.
    cols = [x for x in set([elt[0] for elt in snps])]
    cols.sort()
    ## The LAST row of the alignment is the reference.
    alignment = [['']*len(cols) for x in range(len(conf.GENOMEDIFF_FILES)+1)]
    ## make an index of positions to columns.
    for elt in snps:
        i = conf.GENOMEDIFF_FILES.index(elt[4])
        j = cols.index(elt[0])
        alignment[-1][j] = elt[1] ## the reference sequence.
        alignment[i][j] = elt[2]
    # now fill in the empty entries in the matrix with the reference seq value.
    ref = alignment[-1]
    print ref
    for i in range(len(alignment)):
        for j in range(len(cols)):
            if alignment[i][j] == '':
                alignment[i][j] = ref[j]

    str_alignment = [''.join(x) for x in alignment]
    aln_ids = [os.path.splitext(gd)[0] for gd in conf.GENOMEDIFF_FILES]
    aln_ids = aln_ids + [ref_record.id] # add the reference. 
    site_recs = [SeqRecord(Seq(x), id=y) for x,y in zip(str_alignment, aln_ids)] 
    ## turn into a Biopython Alignment object.
    msa = MultipleSeqAlignment(site_recs)
    return msa

def AlignmentToPhylogeny(aln):

    ## write alignment to a temporary file.
    out_handle = open("temp/aln.phy", "w")
    AlignIO.write(aln, out_handle, "phylip")
    out_handle.close()
    ## raxml needs phylip formatted data.       
    raxml_cline = RaxmlCommandline(sequences="temp/aln.phy", threads=2, model="GTRGAMMA", name="test")
    raxml_cline()
    return None


