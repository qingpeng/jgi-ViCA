
#ViCA
##   Classifying virus from metagenomic and metatransciptomic contigs
 


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
* [khmer v1.4](https://pypi.python.org/pypi/khmer/1.4/)
* Pfam/Vfam HMMER
pip install khmer==1.4
* [Spark]

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


# User Case

With this package, a model is offered with training using simulated data from
RefSeq genomes. Tools are provided if the users want to train the model 
themselves with their own data. 


## Model Tuning
split into training and testing
```angular2html
scripts/5_create_training_testing_with_seq_name.py
```
model training
```angular2html
spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```

model evaluation
```
spark-submit ./scripts/spark_evaluating_model.py testing.vect training.vect_model/ testing.vect.prediction testing.vect.report testing.vect.prc.png
```


## Training

model training (small data set)
```angular2html
spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```

## Prediction
There are three use cases for doing the prediction:

- Large scale prediction - pipeline (in NextFlow) used for prediction on large 
number of sequences using HPC or Cloud system

a. feature extraction using nextflow workflow management
```angular2html
scripts/feature_extraction.nf
```
b. using spark to do prediction on the vectors
```angular2html
scripts/spark_prediction.py
```

- Small scale prediction  - downloadable package used for prediction on small
number of sequences running locally (like a laptop)
```angular2html
~/scripts/prediction_pipeline_lite.py
usage: prediction_pipeline_lite.py [-h]
                                   input_file output_file genemark_path
                                   hmmer_path hmmer_db spark_path feature_file
                                   model_directory

```


- Web Application - a web interface where the users can submit small number of
sequences for prediction
```angular2html
~/web/server.py
```

### Other helper scripts
select vectors from specific features:
```angular2html
scripts/pick_vectors_by_feature.py
```

## Diagram of program

![Diagram](./doc/images/vica.png)

## Description of scripts

scripts/sklearn_training_model.py : training model using sklearn, with weighted option
scripts/sklearn_evaluate_model.py: evaluate model using sklearn

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


