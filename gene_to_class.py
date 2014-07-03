#!/usr/bin/env python
#attempt to switch gene info over to a Class object

class Gene(object):

    def __init__(self, my_start, my_end, gene_name, locus_tag, A, G, C, T):
        self.my_start = my_start
        self.my_end = my_end
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.A = 0
        self.G = 0
        self.C = 0
        self.T = 0

    def __str__(self):
        return '{} ({}): ({},{}) A={} G={} C={} T={}'.format(self.gene_name, self.locus_tag, self.my_start, self.my_end, self.A, self.G, self.C, self.T)

def get_genecoordinates(record):
    geneinfo = {}
    refseq = list(record.seq)
    for i in record.features:
        if i.type == 'CDS':
            Gene.my_start = int(i.location.start)
            Gene.my_end = int(i.location.end)
            Gene.locus_tag = i.qualifiers['locus_tag'][0]
            try:
                Gene.gene_name = i.qualifiers['gene'][0]
            except KeyError:
                Gene.gene_name = 'none'
            for base in refseq[Gene.my_start:Gene.my_end]:
                Gene._ = Gene._.get(base,0) +  1
            geneinfo.append(Gene)
    print '\n', 'Reference genome list (1st 10): gene name, locus_tag, start, stop, A, G, C, T', '\n', geneinfo[:10]
    return geneinfo
   

