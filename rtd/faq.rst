FAQ
===

**Who should I contact if I encounter a problem?**

If you believe your problem may be encountered by other users,
please post the question on `biostars <https://www.biostars.org/>`_.
Check to make sure your question has not been already answered 
by looking at posts with the tag `CHASMplus <https://www.biostars.org/t/CHASMplus>`_.
Otherwise, create a new post with the CHASMplus tag. We will be checking
biostars for questions. You may also contact me directly at
ctokheim AT jhu dot edu.

**Are the p-values by CHASMplus valid for targeted gene panels?**

The p-values reported from CHASMplus are based on whole-exome
sequencing studies. If your mutations comes from a target gene panel, CHASMplus
cannot capture ahead of time what your specific genes are. To get an accurate
estimate of statistical significance, you will need to use the source code version
of CHASMplus to perform a customized analysis. Documentation on how to do this will be added in the future.

**Where can I obtain the training data for CHASMplus?**

You can obtain the set of mutations used for training from `here <http://karchinlab.org/data/CHASMplus/formatted_training_list.txt.gz>`_.

**I want to compare my method to CHASMplus. How should I do it?**

I recommend using the precomputed scores available through OpenCRAVAT [see quickstart-ref_]. Scores in the precompute were generatured using gene-hold out cross-validation, so there is no issue when evaluating performance about training set overlap leading to overfitting. However, the scores do reflect training based on data from The Cancer Genome Atlas (TCGA). If a new method is trained using more data than is available from the TCGA, then it is recommended to create a new CHASMplus model based on the larger data set by using the CHASMplus source code.
