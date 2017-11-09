Installation
------------

.. image:: https://travis-ci.org/KarchinLab/CHASM2.svg?branch=master
    :target: https://travis-ci.org/KarchinLab/CHASM2

Releases
~~~~~~~~

CHASM2 can be downloaded on `github <https://github.com/KarchinLab/CHASM2/releases>`_.

Package requirements
~~~~~~~~~~~~~~~~~~~~

CHASM2 Environment
++++++++++++++++++

We recommend using `conda <https://conda.io/docs/>`_ to install the CHASM2 dependencies.

.. code-block:: bash

   $ conda env create -f environment.yml  # create environment for CHASM2
   $ source activate CHASM2  # activate environment for CHASM2

Make sure the CHASM2 environment is activated when you want to run CHASM2.

20/20+
++++++

You will need to download the `2020plus github repository <https://github.com/KarchinLab/2020plus/releases>`_.

Set the directory of 20/20+ in the configuration file for CHASM2.

.. code-block:: yaml

    twentyTwentyPlus: /path/to/2020plus  # set this directory


Check your PATH variable
========================

Make sure that you have add the 20/20+ directory to your `PATH` variable. If you have done this correctly, the following command should print the location of the 2020plus.py script.

.. code-block:: bash

   $ which 2020plus.py

MySQL
+++++

Features for mutations CHASM2 are obtained  can also be prepared by directly using the MuPIT MySQL database. A MySQL dump of the SNVBox database contains features used for our study. The SNVBox database has a fairly large file size, you may want to directly download and upload to MySQL.

.. code-block:: bash

   $ wget http://karchinlab.org/data/CHASM2/SNVBox_chasm2.sql.gz
   $ gunzip SNVBox_chasm2.sql.gz
   $ mysql [options] < SNVBox_chasm2.sql

This will create a database named mupit_modbase, where [options] is the necessary MySQL parameters to login. You will need sufficient privileges on your MySQL database to CREATE a new database. If everything worked properly, you should see a database named "SNVBox_20161028_sandbox".

SNVBox code
+++++++++++

TODO
