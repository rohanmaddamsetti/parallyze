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
    pass

