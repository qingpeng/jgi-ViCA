Application
=======
- 5_create_training_testing.py is too slow for large vector file, need improving
- finish nextflow pipeline for testing data sets
- migrate MLlib training from RDD-based to DF-based 20170327
- docker support
- web application
- Google Cloud support
- migrate to pure Spark application
    * full pipeline from segment to vectors on spark

web service
=====
to microservice

for now, feature extraction and prediction in the same spark instance
to be separated...

add DB support...

offline -> online prediction, keep model hot... faster speed


Extra
=======



Others
=======
- Spark Version check:
    * careful about the version of Spark   2.0 or 1.5???

Thoughts
======
																											
prediction for one segment directly....


local running... no reftree.... get segments..... from segments....

standalone apps. based on spark..... package in a package... for prediction....


Remove reliablity on reftree, or nextflow.... 

work for small segments, that's it.

it's the user's responsibility to scale......

does not interact with cluster management... except Spark


input: 1 - 10 segments/sequences

output: vectors....


convert vectors for predictions.....

pure python project,,, as much as possible....

### pipeline for large training, full scale... feature extraction..... training...


get all segments from refseq!!

use reftree ????
### pipeline for large prediction, million sequences....

### standalone application, feature extraction training...... only generate model small batch job.

### standalone application, for multiple sequences....  prediction only.... small batch job, local or hpcc.... 

extraction... 

nextflow for runnign local..... slow pipeline.....




loading model...

prediction...

### web application, prediction only, small batch job...


extraction... 

loading model...

prediction...




## docker

## web interface

## further

cassandra.. db
store, retrieve data….

web service…
build server..
scalable….
can scale….