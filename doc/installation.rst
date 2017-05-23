.. GeneLearn documentation master file, created by
   sphinx-quickstart on Thu Jul  9 13:38:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GeneLearn's documentation!
=====================================

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`




# Dependencies
* [GenemarkS version 4.29](http://exon.gatech.edu/GeneMark/)
* [RefTree](https://bitbucket.org/berkeleylab/jgi_reftree)
* [Task Farmer](http://jgi.goe.gov)
* [Python v2.74](https://www.python.org/)
* [Scikit-learn](https://scikits.appspot.com/scikit-learn)
* [Biopython](http://biopython.org)
* [simplejson](https://github.com/simplejson/simplejson)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [khmer v1.4](https://pypi.python.org/pypi/khmer/1.4/) pip install khmer==1.4
* [Pfam/Vfam HMMER]
* [Spark 2.1.0]

# Preparation
## Spark:

You may come cross warnings like below if you run it on your laptop
```
17/05/04 14:34:34 WARN BLAS: Failed to load implementation from: com.github.fommil.netlib.NativeSystemBLAS
17/05/04 14:34:34 WARN BLAS: Failed to load implementation from: com.github.fommil.netlib.NativeRefBLAS
```
This is because the native BLAS is not used. This may affect the speed of
model training and prediction.

If you want to avoid this, you may want to build the Spark from source code.

Please refer to this page:
http://www.spark.tc/blas-libraries-in-mllib/


HMMER installation
=======
Tutorial:
http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf

on genepool:

module load hmmer/3.1b2


Pfam database
========
You can download Pfam database from
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz

Vfam database
======
http://derisilab.ucsf.edu/software/vFam/vFam-B_2014.hmm

Prepare Pfam/Vfam database for HMMER run
====
hmmpress Pfam-A.hmm
hmmpress vFam-B_2014.hmm
