
#Gene Pattern Classifier
## A program to classify genes and genomes based on the codon transition \
probabilities calculated from gene calling. 
 
##Use cases:
* Classifying metagenomic and metatransciptomic contigs taxonomically
* Creating data files for genome binning using emergent self organizing maps (ESOM)
* Identifying metagenomic contigs and metatransciptomic contigs that come from deeply \
  divergent forms of life

# Dependencies
[GenemarkS version 4.29](http://exon.gatech.edu/GeneMark/)
[wget](http://www.gnu.org/software/wget/)
[Python v2.74](https://www.python.org/)
[Scikit-learn](https://scikits.appspot.com/scikit-learn)
[SQLite3](https://www.sqlite.org/index.html)
[Biopython](http://biopython.org)
[simplejson](https://github.com/simplejson/simplejson)
[numpy](http://www.numpy.org/)
[scipy](http://www.scipy.org/)
[matplotlib](http://matplotlib.org/)

## Files
###Setup
 * env.sh - source environment variables
 * taxdmpsetup.sql - Creates a sqlite database from NCBI taxonomy files
 * gm_parameters/ A directory with modified GeneMark Parameter files
 8 config.json, config2.json - Config files for running the classifier 
### Creating Training data 
 * create_training_data.py - A wrapper to run all tasks necessary to create training data
 	1. download_training_data.sh - Fetches RefSeq training data
 	2. genome_parserdb.py - Splits RefSeq data into genomes
 	3. genemarksub.sh - A script to submit all fastas in a directory for gene calling by \
 	   genmark as a grid array job
 	4. create_vectors.py - A script to parse a directory with directories containing one \
 	   directory per training class and 1 metamark .mod file per genome into a into a pickle \
 	   file containing a data matrix.
 * metamarkparser.py - a script to parse metamark MOD files info a character state \
 	vector, 1 vector per file in directory. Used for  classifying the experimental data 

###Classifing via SVM
* testSVM.py - A script used to test ideas about how to perform multiple training tasks \
and analyze the perfomance using cross validation. Specifically the script:
   	1. load a pickle file containing training data for a dataset
   	2. Creates a two class and  12 class model
   	3. Normalizes the data runs a randomized grid parameter Search\
   	4. to optimize the parameters
   	5. cross validate the optimized model and print Statistics, ROC curve and Contingency \
   	   table (Confusion Matrix)
* trainSVM.py - A script to take training and test vectors as picles and predict their \
   class and plot a confusion matrix
 	
