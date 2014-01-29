parallyze
=========

software to analyze parallel genome evolution


NOTES

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
