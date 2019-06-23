FAQ
===

**Who should I contact if I encounter a problem?**

If you believe your problem may be encountered by other users,
please post the question on `biostars <https://www.biostars.org/>`_.
Check to make sure your question has not been already answered 
by looking at posts with the tag `CHASMplus <https://www.biostars.org/t/CHASMplus>`_.
Otherwise, create a new post with the CHASMplus tag. We will be checking
biostars for questions. You may also contact me directly at
ctokhei1 AT alumni DOT jh DOT edu.

**Are the p-values by CHASMplus valid for targeted gene panels?**

The p-values reported from CHASMplus are based on whole-exome
sequencing studies. If your mutations comes from a targeted gene panel, CHASMplus
cannot capture ahead of time what are the specific genes being assessed. To get an accurate
estimate of statistical significance, you will need to use the source code version
of CHASMplus to perform a customized analysis. Documentation on how to do this will be added in the future.

**Will CHASMplus support calibrated p-values for targeted gene panels in the future?**

Yes! We are currently trying to include a calibrated p-value calculation for the MSK-IMPACT gene panel into the easy to use OpenCRAVAT software. Since this depends on the particular genes in the targeted gene panel, we would also be willing to include other panels by popular request.

**Where can I obtain the training data for CHASMplus?**

You can obtain the set of mutations used for training from `here <http://karchinlab.org/data/CHASMplus/formatted_training_list.txt.gz>`_.

**I want to compare my method to CHASMplus. How should I do it?**

I recommend using the precomputed scores available through OpenCRAVAT [see :ref:`quickstart-ref`]. Scores in the precompute were generatured using gene-hold out cross-validation, so there is no issue when evaluating performance about training set overlap leading to overfitting. However, the scores do reflect training based on data from The Cancer Genome Atlas (TCGA). If a new method is trained using more data than is available from the TCGA, then it is recommended to create a new CHASMplus model based on the larger data set by using the CHASMplus source code.

**I want to apply CHASMplus to new data. How should I do it?**

For small datasets it is recommended that pre-computed scores obtained from OpenCRAVAT are used. Large datasets may also use the pre-computed scores, but won't benefit from predictions that are customized to your data. If a cancer type you are interested in does not have precomputed scores or you have collected a large number of cancer samples, than it is recommended to use the CHASMplus source code to perform tailored predictions for your data.
