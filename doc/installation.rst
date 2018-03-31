Installation
------------

.. image:: https://travis-ci.org/KarchinLab/CHASMplus.svg?branch=master
    :target: https://travis-ci.org/KarchinLab/CHASMplus

Releases
~~~~~~~~

CHASMplus can be downloaded on `github <https://github.com/KarchinLab/CHASMplus/releases>`_.

Package requirements
~~~~~~~~~~~~~~~~~~~~

CHASMplus Environment
+++++++++++++++++++++

We recommend using `conda <https://conda.io/docs/>`_ to install the CHASMplus dependencies.

.. code-block:: bash

   $ conda env create -f environment.yml  # create environment for CHASMplus
   $ source activate CHASMplus  # activate environment for CHASMplus

Make sure the CHASMplus environment is activated when you want to run CHASMplus.

20/20+
++++++

You will need to download the `2020plus github repository <https://github.com/KarchinLab/2020plus/releases>`_. Please follow the installation instructions from the `20/20+ website <http://2020plus.readthedocs.io/>`_.

Set the directory of 20/20+ in the configuration file for CHASMplus. You can find this configuration file within the CHASMplus directory at chasm2/data/config.yaml.

.. code-block:: yaml

    twentyTwentyPlus: /path/to/2020plus  # set this directory


Check your PATH variable
========================

Make sure that you have add the 20/20+ directory to your `PATH` variable. If you have done this correctly, the following command should print the location of the 2020plus.py script.

.. code-block:: bash

   $ which 2020plus.py

SNVBox database (MySQL)
+++++++++++++++++++++++

Features for mutations CHASMplus are obtained  can also be prepared by directly using a MySQL database. A MySQL dump of the SNVBox database contains features used for our study. The SNVBox database has a fairly large file size, you may want to directly download and upload to MySQL.

.. code-block:: bash

   $ wget http://karchinlab.org/data/CHASMplus/SNVBox_chasmplus.sql.gz
   $ gunzip SNVBox_chasmplus.sql.gz
   $ mysql [options] < SNVBox_chasm2.sql

This will create a database named mupit_modbase, where [options] is the necessary MySQL parameters to login. You will need sufficient privileges on your MySQL database to CREATE a new database. If everything worked properly, you should see a database named "SNVBox_20161028_sandbox".

SNVBox code
+++++++++++

The next step is to download the code that fetches features from the SNVBox database. Please download the code from `here <http://karchinlab.org/data/CHASMplus/SNVBox.tar.gz>`_, or use wget:

.. code-block:: bash

   $ wget http://karchinlab.org/data/CHASMplus/SNVBox.tar.gz

The next step is to set the configuration file (snv_box.conf) to point towards the established database in the previous section. Specifically, change the db.user, db.password, and db.host to point towards your own mysql user name, mysql password, and mysql host.

The last step is to set the CHASMplus configuration file to point towards the path of the snvGetGenomic command within the SNVBox code. The yaml configuration file is found within the CHASMplus directory at chasm2/data/config.yaml.

.. code-block:: yaml

    snvGetGenomic: /path/to/SNVBox/snvGetGenomic  # set this path
