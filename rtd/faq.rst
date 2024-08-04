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

**Does CHASMplus support targeted gene panels?**

Yes! We have added CHASMplus modules into `OpenCRAVAT <https://opencravat.org/>`_ that support targeted gene panel sequencing. The first iteration was done for the MSK-IMPACT gene panel, but others may be supported upon request. Please see the OpenCRAVAT instructions for how to install "annotators" (`GUI <https://github.com/KarchinLab/open-cravat/wiki/5.-GUI-usage#managing-modules>`_, `command line <https://github.com/KarchinLab/open-cravat/wiki/1.-Installation-Instructions#install-annotators>`_). The MSK-IMPACT version of CHASMplus will have "MSK-IMPACT" (GUI version) or "mski" (command line version) in the name.

**Can I get custom scores based on my own data from targeted gene panels?**

Yes. However, you will need to run the source code version of CHASMplus to customize the predictions to your data.

**Where can I obtain the training data for CHASMplus?**

You can obtain the set of mutations used for training from `here <https://www.dropbox.com/scl/fi/i0018dqgtcjq3n1ebkbuq/formatted_training_list.txt.gz?rlkey=f4fuf4l58cajrwov9n7gcnwry&st=nd23nrdq&dl=1>`_.

**I want to compare my method to CHASMplus. How should I do it?**

I recommend using the precomputed scores available through OpenCRAVAT [see :ref:`quickstart-ref`]. Scores in the precompute were generated using gene-hold out cross-validation, so there is no issue when evaluating performance about training set overlap leading to overfitting. However, the scores do reflect training based on data from The Cancer Genome Atlas (TCGA). If a new method is trained using more data than is available from the TCGA, then it is recommended to create a new CHASMplus model based on the larger data set by using the CHASMplus source code.

**I want to apply CHASMplus to new data. How should I do it?**

For small datasets it is recommended that pre-computed scores obtained from OpenCRAVAT are used. Large datasets may also use the pre-computed scores, but won't benefit from predictions that are customized to your data. If a cancer type you are interested in does not have precomputed scores or you have collected a large number of cancer samples, than it is recommended to use the CHASMplus source code to perform tailored predictions for your data.
