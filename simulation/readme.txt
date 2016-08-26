We randomly selectd 200 genes from Homo sapiens chromosome 1.
We select first 100 genes as DE and rest are NON-DE

In the script ./simulate.sh,
we simluate 10 groups of reads.

To use:
please download the fasta file for chromosome 1 to this folder
and rename it to 1.fa.
We were using ensembl grch37 release-84
make sure the first line of 1.fa begins with:
>1 
as in the 1.fa.example
because we are using the Ensembl convention.

Then run
./simulate.sh
The default number of reads per condition per DE/NON-DE per subject
is 700,000. It may takes a while to generate the reads.


After generation, you can use tophat to align th reads back to genome.
Sort and index the bam files generated.

Then you can use our code to model the data and compute the KL divergence of
each gene between the N and T conditions



References:
In simluation, we use the code RNASeqReadSimulator by Wei Li (li.david.wei AT gmail.com)
downloaded from https://github.com/davidliwei/RNASeqReadSimulator

I modified the code a little so we can specify the random seed 
