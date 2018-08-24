.. CHASMplus documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CHASMplus: predicting driver somatic missense mutations in human cancers
========================================================================

.. figure:: /images/karchinlablogo.png
   :align: right
   :width: 225

:Author: Collin Tokheim, Rachel Karchin
:Contact: ctokheim # jhu DOT edu
:Lab: `Karchin Lab <http://karchinlab.org/>`_ 
:Source code: `GitHub <https://github.com/KarchinLab/CHASMplus>`_
:Q&A: `Biostars (tag: CHASMplus) <https://www.biostars.org/t/CHASMplus/>`_ 

Large-scale cancer sequencing studies of patient cohorts have statistically implicated many cancer driver genes, with a long-tail of infrequently mutated genes. Here we present CHASMplus, a computational method to predict driver missense mutations, which is uniquely powered to identify rare driver mutations within the long-tail. We show that it substantially outperforms comparable methods across a wide variety of benchmark sets.  Applied to 8,657 samples across 32 cancer types, CHASMplus identifies over 4,000 unique driver mutations in 240 genes, further distinguished by their specific cancer types.  Our results support a prominent emerging role for rare driver mutations, with substantial variability in the frequency spectrum of drivers across cancer types.  The trajectory of driver discovery may already be effectively saturated for certain cancer types, a finding with policy implications for future sequencing.  As a resource to handle newly observed driver mutations, we systematically score every possible missense mutation across the genome and provide access to those scores through `OpenCRAVAT <https://github.com/KarchinLab/open-cravat/wiki>`_.

Contents:

.. toctree::
   :maxdepth: 3

   quickstart_opencravat
   models
   download
   installation
   tutorial
   faq

Releases
--------

* `CHASMplus v1.0.0 <https://github.com/KarchinLab/CHASMplus/archive/v1.0.0.tar.gz>`_ - 8/17/2018 - Initial release

Citation
--------

The manuscript is currently submitted. Please cite the biorXiv paper for now:

Tokheim, C., & Karchin, R. (2018). Enhanced context reveals the scope of somatic missense mutations driving human cancers. bioRxiv, 313296. 
