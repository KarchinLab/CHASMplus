Quick start (OpenCRAVAT & CHASMplus)
------------------------------------

The easiest way to obtain CHASMplus scores is by using OpenCRAVAT to fetch precomputed scores. You will need python 3.5 or newer to use OpenCRAVAT.

Install OpenCRAVAT
++++++++++++++++++

You will first need to install the OpenCRAVAT python package, please follow the instructions on the OpenCRAVAT wiki: 

`Installation Instructions <https://github.com/KarchinLab/open-cravat/wiki/1.-Installation-Instructions>`_

Install chasmplus annotator
+++++++++++++++++++++++++++

OpenCRAVAT has a modular architecture to perform genomic variant interpretation including variant impact, annotation, and scoring. CHASMplus is one module available in the CRAVAT store. To install the CHASMplus module within OpenCRAVAT, please execute the following command:

.. code-block:: bash

   $ cravat-admin install chasmplus

The above command may take a couple minutes.

Running CHASMplus
+++++++++++++++++

OpenCRAVAT takes as input either a VCF file or a simple tab-delimited text file. I will describe a simple example that uses the latter. The simple tab-delimited text file should contain a variant ID, chromosome (with "chr"), start position (1-based), strand, reference allele, alternate allele, and optional sample ID.::

    var1	chr10	122050517	+	C	T
    var2	chr11	124619643	+	G	A
    var3	chr11	47358961	+	G	T
    var4	chr11	90135669	+	C	T
    var5	chr12	106978077	+	A	G

You can download an example input file `here <https://raw.githubusercontent.com/KarchinLab/CHASMplus/master/doc/input.txt>`_.

.. note:: By default, OpenCRAVAT processes variants on the hg38 reference genome. If you are using hg19 or hg18, please specify with the "-l" parameter your specific reference genome so that OpenCRAVAT will know to lift over your variants.
   
You can run CHASMplus by using the `cravat` command. For information about command line options, please see the command line help:

.. code-block:: bash

   $ cravat -h

I recommend using the "--sm --ea" options when running CHASMplus, so that large input files will finish much more quickly. To obtain CHASMplus scores for the example input file, run the following command:

.. code-block:: bash

   $ cravat --sm --ea -a chasmplus -d output_directory input.txt

The above command will run the chasmplus annotator (specified by the -a flag) and save results to the directory named "output_directory". Scores and p-values from CHASMplus are found in the "input.txt.chasmplus.var" file. You will also need the "input.txt.crm" file to merge the user provided variant IDs into the CHASMplus results. The .var file should look like this::

    #name=chasmplus
    #displayname=CHASMplus
    #column=0,UID,uid,int
    #column=1,P-value,pval,float
    #column=2,Score,score,float
    #column=3,Transcript,transcript,string
    #column=4,All results,results,string
    #no_aggregate=
    #UID	P-value	Score	Transcript	All results
    1	0.399	0.048	ENST00000453444.6	ENST00000334433.7:(0.025:0.59),ENST00000358010.5:(0.049:0.393),*ENST00000453444.6:(0.048:0.399),NM_001291876.1:(0.046:0.412),NM_001291877.1:(0.045:0.418),NM_206861.2:(0.048:0.399),NM_206862.3:(0.025:0.59)
    2	0.99	0.001	NM_052959.2	*NM_052959.2:(0.001:0.99)
    3	0.446	0.041	NM_001080547.1	ENST00000533968.1:(0.053:0.369),*NM_001080547.1:(0.041:0.446),NM_003120.2:(0.049:0.393)

Interpretation
++++++++++++++

CHASMplus scores range from 0 to 1, with higher scores meaning more likely to be a cancer driver mutation. If you are looking to identify a discrete set of putative driver mutations, then we suggest that you correct for multiple hypothesis testing. We recommend using the Benjamini-Hochberg (BH) procedure for controling the false discovery rate. You will need to use an external package to do this, e.g., the `p.adjust` function in R. False discovery rate adjustments will likely be added in the future.
