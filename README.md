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
- `.\demo.sh`
- Check the output folder and kl.txt file. The kl.txt file will have the output like
```
    chr22 ENSG00000100065  1 1 0.0912106
```

# Instructions

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
