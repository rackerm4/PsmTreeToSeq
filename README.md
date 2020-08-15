# pspTSGen
Pipeline using DendroPy & Seq-Gen.

Simulates pyhlogentic trees with protracted speciation model, using Seq-Gen to simulate the evolution of nucleotide or amino acid sequences along those phylogenies.

Status: Currently not working.

## Background


### Installation & Docker

```sh
$ git clone https://github.com/rackerm4/pspTSGen
$ docker build -t pspTSGen .
$ docker run pspTSGen --profile default --num_runs <N> --schema newick --output data
```
## Requirements

* DendroPy==4.4.0
* future==0.18.2
* PyYAML==5.3.1
* Seq-Gen==xxxx

(Will be installed by pip in Docker)

### Arguments
[`main.py`]
Arg | Notes
------- | --------
--profile/-p    | Choosing parameters file (config): no specification of parameters will result in randomized values
--num_runs/-n   | Enter number of trees you want simulate
--schema/-s | tree schema (newick, nexus..)
--output/-o | specify output directory

## Work to be done
- Seq-Gen installation via Dockerfile
- Passing Seq-Gen parameter from default config
- randomized Seq-Gen parameters
- several testing:
    - behavior with high tree counts & full size

## Open questions
- Big data runtime?
- Works better with implementing with nextflow.io ?
    - parallelization of runs
- number of digits after the decimal point for random values

### How does it work

## Resources
### psp-related publications
Estimating the duration of speciation from phylogenies  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4262007/

Prolonging the past counteracts the pull of the present: protracted speciation can explain observed slowdowns in diversification - https://pubmed.ncbi.nlm.nih.gov/21873376/

The reconstructed tree in the lineage-based model of protracted speciation - https://pubmed.ncbi.nlm.nih.gov/24615006/

### seq-gen-related publications
Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees - https://pubmed.ncbi.nlm.nih.gov/9183526/
    - Full article for free: https://academic.oup.com/bioinformatics/article-pdf/13/3/235/1170463/13-3-235.pdf﻿
### Biology
https://en.wikipedia.org/wiki/Phylogenetics

### Bioinformatics
https://en.wikipedia.org/wiki/Bioinformatics


