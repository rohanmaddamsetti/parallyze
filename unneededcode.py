def base_to_int(i):
    if i == 'A':
        return 0
    elif i == 'G':
        return 1
    elif i == 'C':
        return 2
    else:
        return 3

def int_to_base(i):
    if i == 0:
        return 'A'
    elif i == 1:
        return 'G'
    elif i == 2:
        return 'C'
    else:
        return 'T'

def seq_to_int(seq):
    converted_to_int=[base_to_int(b) for b in seq]
    return converted_to_int

def int_to_seq(seq):
    converted_to_seq=[int_to_base(b) for b in seq]
    return converted_to_seq

###################

def get_genecoordinates(record):
    #input: record file from reference (i.e., ref as record file)
    #get all IDed genes and their beginning and end position;
    #store in list of tuples
    geneinfo = []
    for i in record.features:
        if i.type == 'CDS': #or i[type] if dictionary - tab everything below
        #add non-CDS records later? (tRNA, lonely 'gene', repeat_region, misc_feature
            my_start = int(i.location.start)
            my_end = int(i.location.end)
            locus_tag = i.qualifiers['locus_tag'][0]
            try:
                 gene_name = i.qualifiers['gene'][0]
            except KeyError: 
                gene_name = 'none'
            mytuple = (my_start, my_end, gene_name, locus_tag)
            geneinfo.append(mytuple)
    print '\n', 'Reference genome gene list (1st 10): start position, end pos, gene, locus_tag', '\n', geneinfo[:10]
    return geneinfo

#####################

def gds_gene_rank(filenames, params):
    #input: mutations, from parse_gdfiles
    #given user input #x, list x most mutated genes across all gdfiles
    #(#muts/gene across all lines)
    mut_genes = {}
    rank_mut_types = ['nonsynonymous']
    if params['synonymous'] == 1:
        rank_mut_types.append('synonymous')
    if params['noncoding'] == 1:
        rank_mut_types.append('noncoding')
    if params['pseudogene'] == 1:
        rank_mut_types.append('pseudogene')
    if params['intergenic'] == 1:
        rank_mut_types.append('intergenic')
    for fname in filenames: 
        for mut in filenames[fname]:
            if mut['mut_type'] == 'SNP' and mut['snp_type'] in rank_mut_types:
                for tag in mut['locus_tag']:
                    mut_genes[tag] = mut_genes.get(tag, 0) + 1
    print mut_genes
    return mut_genes
    #sorted_mut_genes = sorted(mut_genes.iteritems(), key=operator.itemgetter(1), reverse = True)
    #mut_genes_number = int(len(sorted_mut_genes))
    #print 'The', params['number_of_top_genes'], 'most mutated genes of all', mut_genes_number, 'mutated genes:'
    #print sorted_mut_genes[:params['number_of_top_genes']]
    #return sorted_mut_genes
    ##diff btwn intergenic and noncoding (has no new base)? pseudogene? all exclusive?

####################
