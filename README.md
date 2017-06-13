# riboWaltz

R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
------------------------------------------------------------------------

### Overview

Ribosome profiling is a powerful technique used to study translation at the genome-wide level, generating unique information concerning ribosome positions along RNAs. Optimal localization of ribosomes requires the proper identification of the ribosome P-site in each ribosome protected fragment, a crucial step to determine trinucleotide periodicity of translating ribosomes and draw correct conclusions concerning where ribosomes are located. To determine the P-site within ribosome footprints at nucleotide resolution, the precise estimation of its offset with respect to the protected fragment is necessary. 

riboWaltz is an R package for calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data. Taking advantage of a two-step algorithm where offset information is passed through populations of reads with different length in order to maximize offset coherence, riboWaltz computes with high precision the P-site offset. riboWaltz also provides a variety of graphical representations, laying the foundations for further accurate RiboSeq analyses and better interpretation of positional information.

------------------------------------------------------------------------

### License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

### Installation

To install `riboWaltz`, directly from github: 
    
    install.packages("devtools")
    library("devtools")
    install_github("LabTranslationalArchitectomics/RiboWaltz")

------------------------------------------------------------------------

### Contact

fabio.lauria@gmail.com
t.tebaldi@unitn.it
gabriella.viero@gmail.com
