# CHASMplus

## About

Large-scale cancer sequencing studies of patient cohorts have statistically implicated many cancer driver genes, with a long-tail of infrequently mutated genes. Here we present CHASMplus, a computational method to predict driver missense mutations, which is uniquely powered to identify rare driver mutations within the long-tail. We show that it substantially outperforms comparable methods across a wide variety of benchmark sets. Applied to 8,657 samples across 32 cancer types, CHASMplus identifies over 4,000 unique driver mutations in 240 genes, further distinguished by their specific cancer types. Our results support a prominent emerging role for rare driver mutations, with substantial variability in the frequency spectrum of drivers across cancer types. The trajectory of driver discovery may already be effectively saturated for certain cancer types, a finding with policy implications for future sequencing. As a resource to handle newly observed driver mutations, we systematically score every possible missense mutation across the genome and provide access to those scores through [OpenCRAVAT](https://github.com/KarchinLab/open-cravat/wiki).

## Scoring mutations with CHASMplus

We recommend that most users who just want to obtain CHASMplus predictions use OpenCRAVAT to obtain scores. A detailed walk though is available on the [quick start page](https://chasmplus.readthedocs.io/en/latest/quickstart_opencravat.html) of the [CHASMplus website](https://chasmplus.readthedocs.io/).

## Source code releases

You can download [releases](https://github.com/KarchinLab/CHASMplus/releases) on github.

## Installation

### conda

We recommend using conda to install the CHASMplus dependencies.

```bash
$ conda env create -f environment.yml 
$ source activate CHASMplus
```

Make sure the CHASMplus environment is activated when you want to run CHASMplus.

### pip

CHASMplus needs R and the randomForest R package. Since pip can not do this, you will need to do this by yourself.

CHASMplus also requires the following python packages:

* numpy
* scipy
* pandas>=0.17.0
* scikit-learn
* rpy2
* probabilistic2020
* snakemake
* pyyaml

To install these packages via `pip` you can use the following command:

```bash
$ pip install -r requirements.txt
```

## Platform

CHASMplus is only intended to run on *linux* operating systems and on a compute server.

## Issues

* snvGetGenomic does not work in python 3 so you will also need to have python2.7 installed (check `which python2.7`). This will be alleviated in a later release.
