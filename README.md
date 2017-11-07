# CHASM2

## About

Sequencing studies have statistically implicated genetic drivers of human cancers by distinguishing these from the expected random accumulation of somatic mutations. However, prior work has coarsely focused on driver genes or regions, largely avoiding prediction of individual mutations. Here, we develop and rigorously validate CHASM2 to predict individual driver somatic missense mutations and show it exceeds state-of-the-art performance. Applied to 32 cancer types, CHASM2 identifies 3,527 unique drivers with four times higher prevalence of rare drivers than previously calculated. Our results indicate a complex relationship between the driver landscape of somatic missense mutations and each cancer type, some reveal a prominent role for rare drivers while others rely on common drivers and are already saturating discovery. We show experimentally that CHASM2 discriminates radiosensitivity mutations within the ATM gene, a phenotype of potential clinical relevance. The prevalence of rare cancer drivers has implications for future interpretation of cancer genomes and genome-driven oncology.

## Releases

You can download [releases](https://github.com/KarchinLab/CHASM2/releases) on github.

## Installation

### conda

We recommend using conda to install the CHASM2 dependencies.

```bash
$ conda env create -f environment.yml 
$ source activate CHASM2
```

You will need to activate the CHASM2 environment each time you want to run CHASM2.

### pip

CHASM2 needs R and the randomForest R package. Since pip can not do this, you will need to do this by yourself.

CHASM2 also requires the following python packages:

* numpy
* scipy
* pandas>=0.17.0
* scikit-learn
* rpy2
* probabilistic2020
* snakemake

To install these packages via `pip` you can use the following command:

```bash
$ pip install -r requirements.txt
```

## Issues

* snvGetGenomic does not work in python 3, which is what most of my code is intended to run on
* 20/20+ only works in hg19, but CHASM2 is in hg38
