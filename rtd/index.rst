.. CHASMplus documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CHASMplus: predicting which missense mutations drive human cancers
================================================================

.. figure:: /images/karchinlablogo.png
   :align: right
   :width: 225

:Author: Collin Tokheim, Rachel Karchin
:Contact: ctokhei1 # alumni DOT jh DOT edu
:Lab: `Karchin Lab <http://karchinlab.org/>`_ 
:Source code: `GitHub <https://github.com/KarchinLab/CHASMplus>`_
:Q&A: `Biostars (tag: CHASMplus) <https://www.biostars.org/t/CHASMplus/>`_ 

Large-scale DNA sequencing studies of patients' tumors have revealed that most driver mutations occur only in a few patients, which presents a challenge for precision medicine. CHASMplus is a machine learning method that accurately distinguishes between driver and passenger missense mutations, even for those found at low frequencies or are cancer type-specific. Unlike previous approaches that focus on identifying driver genes, CHASMplus identifies whether individual mutations are cancer drivers. CHASMplus can be used by both bioinformaticians and biolgists by using a graphical user interface or a command line tool.


.. note:: CHASMplus is available through a graphical user interface [see the :ref:`quickstart-ref`]

Prominent papers using CHASMplus:

* Reiter *et al.*, Minimal functional driver gene heterogeneity among untreated metastases. **Science**
* Anagnostou, Niknafs *et al.*, Multimodal genomic features predict outcome of immune checkpoint blockade in non-small-cell lung cancer. **Nature Cancer**
* Saito *et al.*, Landscape and function of multiple mutations within individual oncogenes. **Nature**
* *Reiter et al.*, An analysis of genetic heterogeneity in untreated cancers. **Nature Reviews Cancer**

Contents:

.. toctree::
   :maxdepth: 3

   quickstart_opencravat
   models
   download
   installation
   faq

Citation
--------

Please cite our paper:

Tokheim and Karchin, CHASMplus Reveals the Scope of Somatic Missense Mutations Driving Human Cancers, Cell Systems (2019), https://doi.org/10.1016/j.cels.2019.05.005
