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


parser = argparse.ArgumentParser(description='A script to train a radial basis function kernel SVM model')
parser.add_argument('-i','--input', help='Input training matrix', type=argparse.FileType('r'), default='-' )
parser.add_argument('-o,','--output',help= 'Output model',type=argparse.FileType('w'), default='-' )
parser.add_argument('-r,','--reportdir',help= 'Directory for reports of he model' ,  default='svmresults' )
parser.add_argument('-c,','--config',help= 'config file', default='config.json' )
args = parser.parse_args()

## Read the configuration file
config = json.load( open(args.config, 'r'))
reportdir = os.path.abspath(args.reportdir)


## Functions
def report(grid_scores, n_top=3):
    """Return best scores after parameter grid search"""
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")

def plot_confusion_matrix(cm, title='Confusion matrix', cmap=plt.cm.Blues, plotname):
    """Plots a confusion matrix (cross tabulation) of classification assignments""" 
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    pl.savefig(plotname, format="pdf")

def randomized_grid_search_rbf(C=4.6, gamma=0.01, n=20, datamat, classvect, n_jobs=1 ):
    """Retun best parameters from randomized grid search"""
    clf = svm.SVC(kernel = 'rbf',cache_size=4000, class_weight='auto')
    param_dist = {'C': scipy.stats.expon(scale=C), 'gamma': scipy.stats.expon(scale=gamma),}
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

## Run Parameter Search cross-validation and write the model

# Open Data 
classvect = [] # Classes for training
modvect = [] # list of training vectors
for line in args.output:
    vect = line.split()
    classvect.append(line[0])
    modvect.append(vect[1:])
modt = np.array(modtrain) # Convert vector list to numpy array

# Normalize and scale data
scaler = preprocessing.StandardScaler().fit(modt)
modtscaled = scaler.transform(modt)

# Open report document
reportname = os.path.join(reportdir, "report.txt")
report = open(reportname, "w")

# Perform parameter search
random_search = randomized_grid_search_rbf(4.6,0.01,20, modtscaled, classvect,1)
params = get_params(random_search)
report.write("Randomized Grid Results")
report.write("\n")
report.write(report(random_search.grid_scores_))

# Fit the model
clf = svm.SVC(kernel='rbf', C=params["C"], gamma = params["gamma"], class_weight='auto', probability=True, cache_size=4000)
clf.fit(modtscaled, classvect)

#  Fivefold Cross Validation
report.write("Cross validation Results")
report.write("\n")
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='precision')
report.write(scores)
report.write("Precision: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='accuracy')
report.write(scores)
report.write("Accuracy: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='recall')
report.write(scores)
report.write("Recall: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
scores = cross_validation.cross_val_score(clf, modtscaled, classvect, cv=5, scoring='f1')
report.write(scores)
report.write("f1: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2)) 

#Generate ROC curve
rocfile=os.path.join(reportdir,"roc.pdf")
report.write("\n")
report.write("Writing SVM model to the file %s" % modelfile)
plot_roc_curve(datamat=modtscaled, classvect=classvect, C=params["C"], \
    gamma=params["gamma"], plotname=rocfile, reporthandle=report)
    


#write model to file
modelfile = os.path.join(reportdir,"svmmodel.p")
report.write("\n")
report.write("Writing SVM model to the file %s" % modelfile)
pickle.dump(clf,os.pathjoin(reportdir,open(modelfile, 'wb'))
modelfile.close()

#close report document
report.close()

