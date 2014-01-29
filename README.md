parallyze
=========

software to analyze parallel genome evolution

## Algorithm design

###Case 1: Genomes from evolution experiment. Assume all independent, i.e. star phylogeny.

Count all mutations (x_1 + x_2 + ... + x_n), and turn into 4x4 matrix.

"Sort" reference genome positions by base. represent as (base, position), e.g. ('A', 3461000)

    for 1 to N replicates:
        for 1 to n genomes:
          draw x_i mutations from the mutation matrix, and drop onto reference genome.


###Case 2: Simulating an evolution experiment (constant u, i.e. no mutators)

Start with a 4x4 matrix and set lambda = mean number of mutations per genome, and n = number of genomes.

    for 1 to N replicates:
      for 1 to n:
        draw x_i from Poisson(lambda)
        drop x_i onto reference genome

###Datasets for Testing


* 40K Clones from LTEE
* Bennett temperature-evolved genomes (doi:10.1126/science.1212986)
* Lieberman et al. Burkholderia outbreak isolates (doi:10.1038/ng.997)

##NOTES

Rich is thinking about a set of tools where
one idea is doing a statistical test for parallelism.
Basically take a set of mutations that happened, either
within a population or across populations, and throw these
mutations down on a reference genome. This kind of work has
been done before--see Tami Lieberman's 2011 Nature Genetics
paper. Rich thinks Jean-Baptiste Michel might be the point man
on this analysis, however. Rich is thinking this method (or a toolbox
of related methods) might make a nice paper, and might be shareable
among research groups.

This is related to my idea for looking for a signature of compensatory
adapatation. Roy's immediate criticism when I talked to him about this,
was that 12 lines do no provide enough statistical power to do this kind
of analysis. Perhaps another function of these simulations could be a
power analysis for the number of lines needed in the experiment to look
for parallelism, given the mutation rate and spectrum, and the length 
of the experiment.

Might this be equivalent to simulating different trees (different
branching patterns corresponding to different times/placement of
mutations on the tree), and calculating likelihoods?

write the parallizer program: given data or a mutation matrix, generate a null
distribution of where mutations fall on a reference genome. Include both coding and
non-coding regions. For the time being, concentrate on SNPs.

Different mutational processes can be superimposed on each other.

Look at Tenaillon paper for hierarchies of parallelism.

for parallel evolution inference: for X lines, after N generations, and rate u,
throw X*N*u mutations onto reference (or do X times onto reference to simulate X
genomes)

This assumes a star phylogeny. Figure out how to weight if non-independent mutations.
