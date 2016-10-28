#!/usr/bin/env python
from sklearn import svm, preprocessing, cross_validation, grid_search
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn.cross_validation import train_test_split
import numpy as np
import cPickle as pickle
import argparse
import os
import scipy.stats
import pylab as pl
import matplotlib.pyplot as plt
import simplejson as json


parser = argparse.ArgumentParser(description='A script to train a SVM model')
parser.add_argument('-i','--train', help='Input training matrix', type=argparse.FileType('r'), default='-' )
parser.add_argument('-t','--test', help = 'test data matrix', type = argparse.FileType('r'), required=True)
parser.add_argument('-r,','--reportdir',help= 'Directory for reports of he model' ,  default='svmresults' )
parser.add_argument('-c,','--config',help= 'config file', default='config.json' )
args = parser.parse_args()

## Read the configuration file
config = json.load( open(args.config, 'r'))["svm"]
reportdir = os.path.abspath(args.reportdir)


## Functions

def plot_confusion_matrix(cm,plotname, title='Confusion matrix', cmap=plt.cm.Blues ):
    """Plots a confusion matrix (cross tabulation) of classification assignments""" 
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    pl.savefig(plotname, format="pdf")

def randomized_grid_search_rbf(datamat, classvect, C=4.6, gamma=0.01, n=20, n_jobs=1 ):
    """Retun best parameters from randomized grid search"""
    clf = svm.SVC(kernel = 'rbf',cache_size=4000, class_weight='auto')
    param_dist = {'C': scipy.stats.expon(scale=C), 'gamma': scipy.stats.expon(scale=gamma)}
    # run randomized search
    random_search = grid_search.RandomizedSearchCV(clf, param_distributions=param_dist,n_iter=n, n_jobs = n_jobs)
    random_search.fit(datamat, classvect)
    return random_search

def randomized_grid_search_linear(datamat, classvect, C=4.6, n=20, n_jobs=1 ):
    """Retun best parameters from randomized grid search"""
    clf = svm.LinearSVC(class_weight='auto')
    param_dist = {'C': scipy.stats.expon(scale=C)}
    # run randomized search
    random_search = grid_search.RandomizedSearchCV(clf, param_distributions=param_dist,n_iter=n, n_jobs = n_jobs)
    random_search.fit(datamat, classvect)
    return random_search
    
def plot_roc_curve(datamat, classvect, C, gamma, plotname, reporthandle):
    """Write ROC parameters to a file and plot a ROC curve"""
    #shuffle and split training and test sets
    random_state = np.random.RandomState(0)
    n_samples, n_features = datamat.shape
    X, y = shuffle(datamat, classvect, random_state=random_state)
    half = int(n_samples / 2)
    X_train, X_test = X[:half], X[half:]
    y_train, y_test = y[:half], y[half:]
    #Run classifier
    classifier = svm.SVC(kernel='rbf',C=C,gamma=gamma, probability=True, class_weight='auto')
    probas_ = classifier.fit(X_train, y_train).predict_proba(X_test)
    #Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
    roc_auc = auc(fpr, tpr)
    reporthandle.write("Area under the ROC curve : %f" % roc_auc)
    #Plot ROC curve
    pl.clf()
    pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('Receiver operating characteristic, Viral vs. Non-viral')
    pl.legend(loc="lower right")
    pl.savefig(plotname, format="pdf")

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

## Run Parameter Search cross-validation and write the model

# Open Data 
classvect, modt = import_matrix(args.train)
testclassvect, testt = import_matrix(args.test)

# Normalize and scale training data
scaler = preprocessing.StandardScaler().fit(modt)
modtscaled = scaler.transform(modt)
testtscaled = scaler.transform(testt)

# Open report document
reportname = os.path.join(reportdir, "report.txt")
reporthandle = open(reportname, "w")


if config["model"] == "rbf":
	# Perform parameter search
	random_search = randomized_grid_search_rbf(datamat=modtscaled, classvect=classvect,C=4.6,gamma=0.01,n=int(config["grid_runs"]),n_jobs=int(config["cores"]))
	params =random_search.best_params_
	reporthandle.write("Randomized Grid Search was run")
	reporthandle.write("\n")
	# Fit the model
	clf = svm.SVC(kernel='rbf', C=params["C"], gamma = params["gamma"], class_weight='auto', probability=True, cache_size=5000)
	clf.fit(modtscaled, classvect)
	
if config["model"] == "linear":
	# Perform parameter search
	random_search = randomized_grid_search_linear(datamat=modtscaled, classvect=classvect,n=int(config["grid_runs"]),n_jobs=int(config["cores"]))
	params =random_search.best_params_
	reporthandle.write("Randomized Grid Search was run")
	reporthandle.write("\n")
	# Fit the model
	clf = svm.LinearSVC(C=params["C"], class_weight='auto')
	clf.fit(modtscaled, classvect)
	
#  Fivefold Cross Validation
reporthandle.write("Cross validation Results")
reporthandle.write("\n")
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='precision_weighted')
# reporthandle.write(scores)
reporthandle.write("Weighted Precision: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                                               
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='recall_weighted')
# reporthandle.write(scores)
reporthandle.write("Weighted Recall: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='f1_weighted')
#reporthandle.write(scores)
reporthandle.write("Weighted f1: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2)) 

#Generate ROC curve
# rocfile=os.path.join(reportdir,"roc.pdf")
# reporthandle.write("\n")
# reporthandle.write("Writing ROC Curve to the file" )
# plot_roc_curve(datamat=modtscaled, classvect=classvect, C=params["C"], \
#    gamma=params["gamma"], plotname=rocfile, reporthandle=report)
    
# Test with real test data

testpredict= clf.predict(testtscaled)
cm = confusion_matrix(testclassvect, testpredict)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
np.set_printoptions(precision=3)
reporthandle.write('Normalized confusion matrix\n')
print(cm_normalized)
plot_confusion_matrix(cm_normalized, os.path.join(reportdir,"cm_norm.pdf") )
print(cm)
plot_confusion_matrix(cm, os.path.join(reportdir,"cm.pdf") )
#write model to file
modelfile = os.path.join(reportdir,"svmmodel.p")
reporthandle.write("\n")
reporthandle.write("Writing SVM model to the file %s" % modelfile)
with open(modelfile, 'wb') as modelhandle:
	pickle.dump(clf,modelhandle)
modelhandle.close()

#close report document
reporthandle.close()

