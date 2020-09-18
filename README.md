[![Travis Build Status](https://travis-ci.org/mmollina/MAPpoly.svg?branch=master)](https://travis-ci.org/mmollina/MAPpoly) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mmollina/mappoly?branch=master&svg=true)](https://ci.appveyor.com/project/mmollina/mappoly)
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/github/mmollina/MAPpoly/branch/master/graphs/badge.svg)](https://codecov.io/github/mmollina/MAPpoly)

![](https://raw.githubusercontent.com/mmollina/MAPpoly/master/mappoly_hexsticker.png)

MAPpoly (v. 0.2.0) is an R package to construct genetic maps in autopolyploids with even ploidy levels. In its current version, MAPpoly can handle ploidy levels up to 8 when using hidden Markov models (HMM), and up to 12 when using the two-point simplification. All the two-point based functions are fast enough to run on standard computers. However, we strongly recommend to use high-performance computation for HMM-based analysis, especially for ploidy levels higher than 4. 

![](https://raw.githubusercontent.com/mmollina/MAPpoly/master/mappoly.gif)


In its current version, MAPpoly can handle three different types of datasets:

1. CSV files 
2. MAPpoly files
  - Dosage based
  - Probability based
3. VCF files

The derivation of the HMM used in MAPpoly can be found in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). Recently, we used MAPpoly to build an ultra-dense multilocus integrated genetic map containing ~30k SNPs and characterized the inheritance system in a sweetpotato full-sib family ([Mollinari et al., 2020](https://doi.org/10.1534/g3.119.400620)). See the resulting map [here](https://gt4sp-genetic-map.shinyapps.io/bt_map/) and the haplotype composition of all individuals in the full-sib population [here](https://gt4sp-genetic-map.shinyapps.io/offspring_haplotype_BT_population/).

# Installation

## From CRAN (stable version - not available yet)

To install MAPpoly from the The Comprehensive R Archive Network (CRAN) use

```R
install.packages("mappoly")
```

## From GitHub (development version)

You can install the development version from Git Hub. Within R, you need to install `devtools`:

```R
install.packages("devtools")
```

If you are using Windows, you must install the the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from Git Hub use

```R
devtools::install_github("mmollina/mappoly", dependencies=TRUE)
```

For further QTL analysis, we recommend our [QTLpoly](https://github.com/guilherme-pereira/QTLpoly) package. QTLpoly is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model [Pereira et al. 2020](https://doi.org/10.1534/genetics.120.303080). 

# Vignettes
* [Building a genetic map using potato genotype data from SolCAP](https://mmollina.github.io/MAPpoly_vignettes/vignette_tetraploid/vignette_tetraploid.html)
* [Building a genetic map in an hexaploid full-sib population using MAPpoly](https://mmollina.github.io/tutorials/hexa_fake/haxaploid_map_construction.html)
* Dataset examples
  * [Hexaploid sweetpotato VCF dataset (Beauregard x Tanzania) obtained using VCF2SM](https://github.com/mmollina/MAPpoly_vignettes/tree/master/data/BT)
  * [Hexaploid simulation with dosage call in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/hexafake)
  * [Hexaploid simulation with dosage probabilities in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/hexafake_geno_dist)
  * [Tetraploid potato with dosage call in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/SolCAP_dosage)
  * [Tetraploid potato with dosage call in CSV format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/tetra_solcap.csv)
  * [Tetraploid potato with dosage probabilities in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/SolCAP)

# Related software

* Simulations
  * [PedigreeSim: Simulation of genetic marker data in diploid and polyploid pedigreed populations.](https://www.wur.nl/en/show/Software-PedigreeSim.htm)

* Genotype calling
  * [ClusterCall: Automated tetraploid genotype calling by hierarchical clustering](https://potatobreeding.cals.wisc.edu/software/)
  * [fitPoly: Genotype Calling for Bi-Allelic Marker Assays](https://CRAN.R-project.org/package=fitPoly)
  * [polyRAD: Genotype Calling with Uncertainty from Sequencing Data in Polyploids and Diploids](https://CRAN.R-project.org/package=polyRAD)
  * [SuperMASSA: Graphical Bayesian inference tool for genotyping polyploids](https://bitbucket.org/orserang/supermassa)
  * [updog: Flexible Genotyping for Polyploids](https://CRAN.R-project.org/package=updog)
  * [VCF2SM: Python script that integrates VCF files and SuperMASSA](https://github.com/guilherme-pereira/vcf2sm)
 
* Genetic mapping in polyploids
  * [MDSMap: High Density Genetic Linkage Mapping using Multidimensional Scaling](https://CRAN.R-project.org/package=MDSMap)
  * [polymapR: Linkage Analysis in Outcrossing Polyploids](https://CRAN.R-project.org/package=polymapR)
  * [TetraploidSNPMap: Linkage maps and mapping QTLs for autotetraploid species, using SNP dosage data.](https://www.bioss.ac.uk/knowledge/tetraploidmap/)
  
  
* Haplotype reconstruction
  * [TetraOrigin:haplotype reconstruction in a full-sib tetraploid family](https://github.com/chaozhi/TetraOrigin)
  * [PolyOriginR:haplotype reconstruction in polyploid multiparental populations](https://github.com/chaozhi/PolyOriginR)

* QTL mapping
  * [QTLpoly: QTL mapping in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model](https://github.com/guilherme-pereira/QTLpoly)
  * [diaQTL: QTL analysis of diploid and autotetraploid diallel populations](https://github.com/jendelman/diaQTL)

# Miscellaneous
* [Supplementary scripts for Mollinari and Garcia (2019)](https://github.com/mmollina/Autopolyploid_Linkage)
* [Workshop: Polyploid Genetic Data Analysis: From Dosage Calling to Linkage and QTL Analysis](http://152.1.45.19/esalq_2019.html)
* [Miscellaneous scripts](https://github.com/mmollina/MAPpoly_vignettes/blob/master/README.md)

# Acknowledgment 

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838), both funded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/).
