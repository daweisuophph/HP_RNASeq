***DEISoM*** is software using hierarchical Bayesian model for ***D***ifferentially ***E***xpressed ***Iso***form detection from ***M***ultiple biological replicates for identifying DE (Differentially Expressed) isoforms from multiple biological replicates representing two conditions, e.g., multiple samples from healthy and diseased subjects.


#Introduction

DEIsoM first estimates isoform expression within each condition by (1) capturing common patterns from sample replicates while allowing individual differences, and (2) modeling the uncertainty introduced by ambiguous read mapping in each replicate. Specifically, we introduce a Dirichlet prior distribution to capture the common expression pattern of replicates from the same condition,  and treat the isoform expression of individual replicates as samples from this  distribution. Ambiguous read mapping is modeled as a multinomial distribution, and ambiguous reads are assigned  to the most probable isoform in each replicate. Second, DEIsoM couples an efficient variational inference and a post-analysis method to improve the accuracy and speed of identification of DE isoforms over  alternative methods. 

# Installation
### Prerequisites
- A C++ compiler e.g. g++ 4.4.7
- [boost v1.53.0+](http://www.boost.org/)
  - For previous versions, please remove -lboost_system in makefile
  - Remember to set BOOST_HOME: `export BOOST_HOME=[path to boost home folder]`

### Third-party software
In the folder third_party, we have prepared the some third_party libraries:
- [liblbfgs](http://www.chokkan.org/software/liblbfgs/)
  - Used for optimization
  - `./configure` and `make`. See its README for more installation instructions.
- [samtools](http://samtools.sourceforge.net/)
  - Used for reading bam files
  - ./configure` and `make`. See its README for more installation instructions.`

### Build
After correctly installing and configure boost, liblbfgs and samtools, we can install DEIsoM now.

First return to the main directory of DEIsoM and then build DEIsom by type `make` in command line for unix-like os.

### Test
After successfully compilation of the code, one can check its installation by testing the demo.

- Change the current working directory to `demo`
- `./demo.sh`
- Check the output folder and kl.txt file. The kl.txt file will have the output like
```
    chr22 ENSG00000100065  1 1 0.0912106
```

# Instructions
### Prepare indexed reference from GFF3 file
For fast computing, DEIsoM first build indexed reference. We build a GFF3 file for every gene seperately from a single gff3 file. We now only support gff3 format. To build the indexed reference, one can run:
```
indexGFF [path to gff3 file] [output indexed reference folder]
```

### Prepare BAM files from FASTQ data
- First, we can use alignment software to map the RNAseq data to our referece. We have tested tophat and RUM. And here is an example of command to run tophat:
```
tophat -p 1 -o [output folder] [reference folder] [pair-end read 1] [pair-end read 2]
```
- Then, we need to sort the accepted_hits.bam (the output of tophat) and make index for them.
```
samtools sort [path to accepted_hits.bam] -o accepted_hits.sorted.bam
samtools index accepted_hits.sorted.bam
```

### Run a job
Suppose we have M replicates and we want to build a model for N genes: gene 1, ..., gene N. The corresponding built indexed references are: gff 1, ..., gff N. We can run the model in one script:
```
run --gene-ids [gene 1], [gene 2], ..., [gene N]   \
    --gffs [gff 1], [gff 2], ..., [gff N]   \
    --bams [bam 1],[bam 2],...,[bam M]   \
    --read-len [read length]   \
    --paired-end [mean length of insert length] [standard derivation of insert length]   \  
    --out-iter [number of iterations for updating alpha]   \
    --in-iter [number of iterations for updating variational parameters]   \
    --output [output folder]
```
The default output is in binary format. If you want to have a human readable format, you can add this option:
```
--human-readable
```
For computing the FPKM for the transcripts, we need you to pass an addition option:
```
--read-count [read count file]
```
read count file is a file of M integers seperated by spaces. Each integer represents the total number of effective reads for each replicate. This can be computed using samtools. If you dont't want to compute the FPKM, this option can be omitted.

### Split jobs for computing on a cluster
Submitting a large number of genes using the above script is not efficient. And we sometimes want to run these jobs in parallel on a large cluster. We have provided a helpder program to automatically deivide the jobs in small chunks:
```
split --trunk-size [trunk size]   \
      --gff-dir [indexed gff directory] \
      --bams [bam 1],[bam 2],...,[bam M]   \
      --read-len [read length]   \
      --paired-end [mean length of insert length] [standard derivation of insert length]   \
      --out-iter [number of iterations for updating alpha]   \
      --in-iter [number of iterations for updating variational parameters]   \
      --output [output folder]
```
You can also pass other options of `run` to it.

### Computing the scores (KL divergences) for identifying DE genes
To evaluate the results, we can use:
```
kl [indexed GFF] [output of group 1] [output of group 2]  
```
for general evaluations (either paired or unpaired between groups)
or
```
kl --beta [indexed GFF] [output of group 1] [output of group 2] 
```
for paired evalutions only.

# Simulations
We have upload the scripts and code for creating sythetic data in our paper. We also include all scripts for all other models that we used in the experiments. Please check `simulation` folder for details.

# License
MIT License

Copyright (c) 2016 Hao Peng (pengh@purude.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Reference Paper

Hao Peng, Yifan Yang, Shandian Zhe, Jian Wang, Michael Gribskov, Yuan Qi, DEIsoM: A hierarchical Bayesian model for identifying differentially expressed isoforms using  biological replicates, in submission.
