
#Gene Structure Classifier
## A program to classify genes and genomes based on the codon transition probabilities calculated from gene calling 
 
##Use cases:
* Classifying metagenomic and metatransciptomic contigs taxonomically
* Creating data files for genome binning using emergent self organizing maps (ESOM)
* Identifying metagenomic contigs and metatransciptomic contigs that come from deeply 
  divergent forms of life

# Dependencies
* [GenemarkS version 4.29](http://exon.gatech.edu/GeneMark/)
* [RefTree](http://jgi.goe.gov)
* [Python v2.74](https://www.python.org/)
* [Scikit-learn](https://scikits.appspot.com/scikit-learn)
* [Biopython](http://biopython.org)
* [simplejson](https://github.com/simplejson/simplejson)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)

## Files
###Setup
 * env.sh - source environment variables
 * gm_parameters A directory with modified GeneMark Parameter files
 * config.json, training.json - config files for running the classifier 

### Creating Training data
 *  taskfarmer
 * single_taxon_training.py - A wrapper to run all tasks necessary to create training data for a single genome
 	1. shred.py - shreds the genome into smaller pieces of fixed length or conforming to a gamma distribution
 	2. Feature extraction script
 	3. feature_formatter.py - a utility script  to format csv feature data into json data
 * add_to_reftree.py - a utility to add all the individual jsons to a local Reftree directory


###model training and cross validation via SVM
* svm_training_cross_validation.py
	1. training_data_formatter.py - a script to take taxonomic nodes or levels and create
	 a multi feature JSON file for classification 
	2. svm_training.py - A script to take training data, do a parameter search, create a model and cross-validate it 

### Sample classification
* svm_classification.py - a wrapper to run SVM classification
	1. feature extraction scripts - see below
	2. feature_formatter.py - a utility script  to format csv feature data into json data
	3. svm_classifier.py a script to load a SVM model and classify sequences

### Producing data for genome binning
* metagenome_binning_data.py -a a wrapper for producing genome binning data
	1. Feature extraction script
	2. feature formatter
	

### Feature extraction ecripts 
* feature_extraction_metamark.py - a wrapper to run metamark and extract the genomic feature
* feature_extraction_kmer.py an alternate method of extracting feature information

### Utility scripts 
* feature_formatter.py

## Diagram of program

![Diagram](GenePatternClassifier.png)


