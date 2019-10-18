# ViCA

Classifying virus from metagenomic and metatransciptomic contigs
 

An updated version of ViCA using Deep Learning approach is hosted at:
https://github.com/USDA-ARS-GBRU/vica

## Usage

With this package, a model is offered with training using simulated data from
RefSeq genomes. 


Tools are provided if the users want to train the model 
themselves with their own data. Please refer to documentation.


There are three use cases for doing the prediction:

### 1. Prediction for large number of sequences 
pipeline (in NextFlow) used for prediction on large 
number of sequences using HPC or Cloud system

#### Step 1. feature extraction using Nextflow workflow management
```angular2html
scripts/feature_extraction.nf
```
#### Step 2. using spark to do prediction on the vectors
```angular2html
$SPARK_PATH/bin/spark-submit spark_prediction.py
usage: spark_prediction.py [-h] libsvm model scaler outfile
```

### 2. Prediction for small number of sequences   
downloadable package used for prediction on small
number of sequences running locally (like a laptop)
```angular2html
~/scripts/prediction_pipeline_lite.py
usage: prediction_pipeline_lite.py [-h]
                                   input_file output_file genemark_path
                                   hmmer_path hmmer_db spark_path feature_file
                                   model_directory scaler_directory
```

### 3. Web Application 
a web interface where the users can submit small number of
sequences for prediction
```angular2html
~/web/server.py
```

## Installation
Please refer to documentation for dependency and installation details

## Diagram of program

![Diagram](./doc/images/vica.png)


