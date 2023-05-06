.. Installation and operation of GEOS-Chem TOMAS on Stetson documentation master file, created by
   sphinx-quickstart on Mon Apr 17 15:35:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Installation and operation of GEOS-Chem TOMAS on Stetson's documentation!
====================================================================================

**GEOS-Chem TOMAS** is a global chemical transport model including the aerosol
microphysics package used to study the emission and climate effect of the aerosols 
in the atmosphere. The following documentation will provide guidance for Dalhousie 
researcher to install and operate the model on the local computationnal server 
of the mathematic and physics department (`Stetson`_).

This documentation is meant to describe the parametrization of GEOS-Chem simulations
that are specific to the stetson server environment and the use of `TOMAS`_ microphysic 
package. This documentation is complementary to the official `Geos-Chem wiki`_ and 
should be used in combination with this original documentation for a full understanding
of the Geos-Chem operation.

.. note::
   This project is under active development.

.. toctree::
   :caption: Stetson User guide
   :maxdepth: 2
   
   guide/Stetson/Stetson.rst

.. toctree::
   :caption: GEOS-Chem TOMAS user guide
   :maxdepth: 2
   
   guide/Gcc_TOMAS/geos-chem.rst

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. URL link reference
.. _Geos-Chem wiki: https://geos-chem.readthedocs.io/en/stable/index.html
.. _TOMAS: http://wiki.seas.harvard.edu/geos-chem/index.php/TOMAS_aerosol_microphysics
.. _Stetson: https://www.mathstat.dal.ca/cluster/doku.php
