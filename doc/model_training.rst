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


## Model Tuning
split into training and testing
```angular2html
scripts/5_create_training_testing_with_seq_name.py
```
model training
```angular2html
$SPARKPATH/bin/spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```

model evaluation
```
$SPARKPATH/bin/spark-submit ./scripts/spark_evaluating_model.py testing.vect training.vect_model/ testing.vect.prediction testing.vect.report testing.vect.prc.png
```


## Training

model training (small data set)
```angular2html
$SPARKPATH/bin/spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```