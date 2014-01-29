parallyze
=========

software to analyze parallel genome evolution by generating null distributions.

## Algorithm design

In all cases, the user should be able to specify only examining coding regions,
or just nonsynonymous and synonymous mutations, or just synonymous mutations.

###Case 1: Genomes from evolution experiment. Assume all independent, i.e. star phylogeny.

Count all point mutations (x_1 + x_2 + ... + x_n), and turn into 4x4 matrix.

"Sort" reference genome positions by base. represent as (base, position), e.g. ('A', 3461000)

    for 1 to N replicates:
        for 1 to n genomes:
          draw x_i mutations from the mutation matrix, and drop onto reference genome.

###Case 2: Dispersion Test. Genomes from evolution experiment. Assume star phylogeny.

    for 1 to N replicates:
        Shuffle all mutations across n genomes.
    Calculate how often a certain dispersion pattern occurs 
    (e.g., 12 mutations in nadR; all mutations occur in separate genomes).
    

###Case 3: Simulating an evolution experiment (constant u, no mutators)

Start with a 4x4 matrix and set lambda = mean number of mutations per genome, and n = number of genomes.

    for 1 to N replicates:
      for 1 to n:
        draw x_i from Poisson(lambda)
        drop x_i onto reference genome
        
###Case 4: Genomes from clinical isolates/epidemics. Also applies to multiple
isolates from the same LTEE population.

1) Infer phylogeny
2) Infer genotypes of all internal nodes by "using parsimonious assumptions" -- or better.
3) Count the number x of inferred mutations, and generate a 4x4 matrix of mutation probabilities.

        for 1 to N replicates:
            drop x mutations onto reference genome, and count number of independent mutations per gene.
        
        average the results to calculate the null distribution.
        
Also, the user should be able to supply a 4x4 matrix of mutation probabilities as an option for this case.