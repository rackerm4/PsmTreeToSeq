# PsmTreeToSeq
Pipeline using DendroPy & Seq-Gen.

Simulates pyhlogentic trees with protracted speciation model, using Seq-Gen to simulate the evolution of nucleotide sequences along those phylogenies.

Status: Currently not working. Errors while parsing parameters. Check "Work to be done"

## :thought_balloon: Background


### Installation & Docker
Install as Module:
```
$ python setup.py install

# Install Seq-Gen
$ sudo apt-get install seq-gen
```

Clone from GitHub:
```sh
$ git clone https://github.com/rackerm4/PsmTreeToSeq.git

$ pip install dendropy (or install module, see above)

$ cd PsmTreeToSeq
$ python main.py --schema nexus -config default --output Dir --num_runs <N>
```
## :wrench: Requirements

* DendroPy==4.4.0
* future==0.18.2
* numpy==1.19.1
* PyYAML==5.3.1
* Seq-Gen

### Arguments
[`main.py`]

Arg | Notes
------- | --------
--config/-c | Choosing parameters file (config): 1 = random parameters, None = None, False = False, or specify your value
--num_runs/-n   | Enter number of trees & sequence files you want simulate
--schema/-s | Tree schema (newick, nexus..)
--output/-o | Specify output directory

## :bulb: Work to be done
- Parsing ProtractedSpeciationProcess parameters from config file not working correctly.
- Dockerfile/images
- Passing Seq-Gen parameter from default config [Done]
- randomized Seq-Gen parameters [Done]
- several testing:
    - behavior with high tree counts & full size
    
## Known issues
- Newick trees raise an error when read with Seq-Gen. Every newick tree is converted into a nexus tree. All newick trees will be deleted.
- Following error:
```
"_Maximum number of runs to execute in the event of prematurely-terminated simulations due to all 
lineages going extinct. Once this number or re-runs is exceed, then TreeSimTotalExtinctionException 
is raised. Defaults to 1000. Set to None to never quit trying._"
```
- Permission denied error? Execute permission for the files in /src/: chmod u+rwx *.py
- If you find yourself with an error like that: 
```
$ docker: Got permission denied while trying to connect to the Docker daemon socket 
try
$ sudo chmod 666 /var/run/docker.sock
```
## :question: Open questions
- Big data runtime
- Number of digits after the decimal point for random values
- typical DNA specific rates

### How does it work

## :books: Resources
- DendroPy - Phylogenetic Computing Library:
    - https://dendropy.org/
- Seq-Gen - Program simulating the evolution of nucleotide or amino acid sequences along a phylogeny
    - http://tree.bio.ed.ac.uk/software/seqgen/

### Protracted speciation model-related publications
Estimating the duration of speciation from phylogenies  
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4262007/

Prolonging the past counteracts the pull of the present: protracted speciation can explain observed slowdowns in diversification 
- https://pubmed.ncbi.nlm.nih.gov/21873376/

The reconstructed tree in the lineage-based model of protracted speciation 
- https://pubmed.ncbi.nlm.nih.gov/24615006/

### Seq-gen-related publications
Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees 
- https://pubmed.ncbi.nlm.nih.gov/9183526/
- Full article for free: 
    - https://academic.oup.com/bioinformatics/article-pdf/13/3/235/1170463/13-3-235.pdf
    
    

##### :microscope: Biology
https://en.wikipedia.org/wiki/Phylogenetics
##### :computer: Bioinformatics
https://en.wikipedia.org/wiki/Bioinformatics


