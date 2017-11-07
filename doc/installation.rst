Installation
------------

.. image:: https://travis-ci.org/KarchinLab/CHASM2.svg?branch=master
    :target: https://travis-ci.org/KarchinLab/CHASM2

Releases
~~~~~~~~

CHASM2 can be downloaded on `github <https://github.com/KarchinLab/CHASM2/releases>`_.

Package requirements
~~~~~~~~~~~~~~~~~~~~

conda
+++++

We recommend using conda to install the CHASM2 dependencies.

.. code-block:: bash

   $ conda env create -f environment.yml  # create environment for CHASM2
   $ source activate CHASM2  # activate environment for CHASM2

You will need to activate the CHASM2 environment each time you want to run CHASM2.


Check your PATH variable
~~~~~~~~~~~~~~~~~~~~~~~~

Make sure that you have add the 20/20+ directory to your `PATH` variable. If you have done this correctly, the following command should print the location of the 2020plus.py script.

.. code-block:: bash

    $ which 2020plus.py
