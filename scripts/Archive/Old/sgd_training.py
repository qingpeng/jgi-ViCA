#!/usr/bin/env python
from sklearn import preprocessing
#from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
#from sklearn.cross_validation import train_test_split
import numpy as np
#import cPickle as pickle
import argparse
import os
import scipy.stats
import pylab as pl
import matplotlib.pyplot as plt
import simplejson as json

from sklearn.linear_model import SGDClassifier

parser = argparse.ArgumentParser(description='A script to train a SGD classifier, optimizing paramerers and estimating accuracy with a test dataset')
parser.add_argument('-i','--train', help='Input training matrix', type=argparse.FileType('r'), default='-' )
parser.add_argument('-t','--test', help = 'test data matrix', type = argparse.FileType('r'), required=True)
parser.add_argument('-r,','--reportdir',help= 'Directory for reports on the model' ,  default='sgdresults' )
parser.add_argument('-c,','--config',help= 'config file', default='config.json' )
args = parser.parse_args()

## Read the configuration file
config = json.load( open(args.config, 'r'))["svm"]
reportdir = os.path.abspath(args.reportdir)

def import_matrix(filehandle):
	classvect = [] # Classes for training
	for i,line in enumerate(filehandle):
		vect = line.split('\t') 
		vect = vect[:0]+ [float(k) for k in vect[1:]] # convert numbers to float
		classvect.append(line[0])
		if i == 0:
			modt = np.array(vect[1:])
		else:
			modt = np.vstack([modt, vect[1:]])
	return (classvect, modt)
	
def plot_confusion_matrix(cm,plotname, title='Confusion matrix', cmap=plt.cm.Blues ):
    """Plots a confusion matrix (cross tabulation) of classification assignments""" 
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    pl.savefig(plotname, format="pdf")
    	
# Open Data 
trainclass, traindat = import_matrix(args.train)
testclass, testdat = import_matrix(args.test)

# Normalize and scale training data
scaler = preprocessing.StandardScaler().fit(traindat)
trainscaled = scaler.transform(traindat)
testscaled = scaler.transform(testdat)

clf = SGDClassifier(loss="hinge", penalty="l2")
clf.fit(trainscaled, trainclass)

testpredictions = clf.predict(testscaled)

cm = confusion_matrix(testclass, testpredictions)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
np.set_printoptions(precision=3)
print(cm_normalized)
plot_confusion_matrix(cm_normalized, os.path.join(reportdir,"cm_norm.pdf") )
print(cm)
plot_confusion_matrix(cm, os.path.join(reportdir,"cm.pdf") )