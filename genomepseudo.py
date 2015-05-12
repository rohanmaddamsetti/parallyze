Parallyze Procedure 1 pseudocode
'''to determine if number of populations (of evol. exp.) in which a particular gene mutates is significantly different from neutral model'''

user inputs 1) reference genome and 2) genome diff files through config file
error check files
(create list of dictionaries?) convert files to matrix (experiment/population/treatment/gene/generation/typemut/position/length)
for each population, create (or search for) mutation type and length (position yet?) file (esp. tally A->T, A->G, etc)
    for i in x  replicates:
    	for i in y populations, randomly distribute:
        	(deletions)
        	(insertions)
        	SNPs
        	(translocations/recombination?)
        	-record 
from user-provided diff files, provide ranking for genes with // mutation hits
let user specify which genes to run a statistical test on ("enter numbers corresponding to gene of interest to test")
if, for the pops of the user, x of the y lines exhibit a mutation in gene z, is the number of replicates with x or more lines exhibiting the same mutation significant?
    (give user some measure of obtained/desired statistical power?) 
