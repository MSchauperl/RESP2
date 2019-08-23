.. resp2 documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RESP2's documentation!
=========================================================

Installation
=========================================================
To install all necessary packages we recommend to follow these steps.

1.) Download the RESP2 package::

   git clone git@github.com:MSchauperl/RESP2.git

2.) Install and activate the conda environment 'RESP2' ::

   cd RESP2
   conda env create -f devtools/conda-envs/RESP2_environment.yaml
   conda activate RESP2

3.) Install the RESP2 package ::

   python setup.py develop

4.) Download the respyte package ::

   cd ..
   git clone https://github.com/lpwgroup/respyte.git

5.) Install respyte package ::

   cd respyte
   python setup.py develop
   cd ..


Example
========================================================

See Jupyter notebook: https://github.com/MSchauperl/RESP2/blob/master/example/Ethanol_Example.ipynb


Functions
========================================================

RESP2 charge generation and ForceBalance Input generation
----------------------------------------------------------
.. automodule:: resp2.resp2
   :members:

Other Modules
==========================

.. automodule:: resp2.create_mol2_pdb
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
