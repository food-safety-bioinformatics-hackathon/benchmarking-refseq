# workflow for testing reference sequence assignment with simulations.
This is a project arising from the **Bioinformatics Food Safety Hackaton 2019**, and the main collaborators are
Boas van der Putten, Martin Lott, Thanh Le Viet and Leonardo de Oliveira Martins. 

The idea is to

1. Simulate reads from a reference genome
2. Use `refseq_masher` to find best candidates

`refseq_masher` has two modes, "contains" and "matches", where the former gives hits over all reads while the later
gives best hits per read. The natural choice for us was to use "contains" (although it is originally designed for
metagenomics), while "matches" assumes input are contigs and not reads. In both cases it uses `mash screen` behind the
scenes.

## current status

We have an instance running on CLIMB using nextflow. A singularity container was also created, to allow it to be reproducible into an HPC environment (at QIB).


