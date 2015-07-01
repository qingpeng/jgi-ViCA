#trainSVM.py
#this script loads genemark data from RefSeq database, scales, optimizes and trains the SVM
#from sklearn import preprosessing
from sklearn import svm
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn import grid_search
import numpy as np
import cPickle as pickle
import argparse
import os
import scipy.stats
from collections import Counter

from time import time
from operator import itemgetter

#for ROC Curve
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc
import pylab as pl

parser = argparse.ArgumentParser(description='A script to train an SVM model')
parser.add_argument('-r','--root', help='a root input directory', required=True)
#parser.add_argument('-c,','--config',help= 'a json formated config file', default='config.json')
args = parser.parse_args()

svmdir = os.path.join(os.path.abspath(args.root),"SVMmatrix.p")
#load objects classvect, modvect
with open(svmdir, 'rb') as pk:
	modvect = pickle.load(pk)
	classvect = pickle.load(pk)
print("classes and elements in class vector:")
print(Counter(classvect))


mod = np.array(modvect)


#create class vector with euks lumped 
classvect_euk = classvect
for n,i in enumerate(classvect_euk):
	if i in ["vertebrate_other","vertebrate_mammalian","invertebrate","plant","fungi","protozoa" ]:
		classvect_euk[n]="eukaryote"
cls_euk = np.array(classvect_euk)
print("classes and elements in eukaryote class vector:")
print(Counter(classvect_euk))


#create viral/non-viral two class model 
classvect_viral = classvect
for n,i in enumerate(classvect_viral):
	if i in ['dsDNAviruses', 'viruses_other', 'ssRNAviruses','dsRNAviruses', 'retroviruses','viruses_other']:
		classvect_viral[n]=True
	else:
		classvect_viral[n]=False
cls_viral = np.array(classvect_viral)
print("classes and elements in viral/non-viral class vector:")
print(Counter(classvect_viral))


#verify that the arrays are of the correct size
def verify_array(mod, cls):
	"""verify that the dimensions of the array are correct and that the data array and the description array match"""
	modshape = mod.shape
	assert modshape[1] == 256, "There does not appear to be 256 elements in the data matrix"
	assert modshape[0] == len(cls), "The data matrix and classification vector have different lengths"
	print("Data array dimensions are " + str(modshape))
	print("Class vector length is " + str(len(cls)))
	print("Data type of data matrix is " + str(mod.dtype))
	print ("Data type of class vector is " + str(cls.dtype))
verify_array(mod, cls_viral)



##normalize
print("Scaling Data")
#modscaled = mod
modscaled = preprocessing.scale(mod)
##validate data


#############Parameter optimization ####################
# Utility function to report best scores
def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")
        

# specify parameters and distributions to sample from
#param_dist = {"C": scipy.stats.expon(scale=4.6), "class_weight": ['auto']}
# clf = svm.SVC(kernel = 'rbf',cache_size=4000, class_weight='auto')             
# param_dist = {'C': scipy.stats.expon(scale=4.6), 'gamma': scipy.stats.expon(scale=.01),}
# # run randomized search
# n_iter_search = 20
# random_search = grid_search.RandomizedSearchCV(clf, param_distributions=param_dist,n_iter=n_iter_search, n_jobs = 4)
# random_search.fit(modscaled, cls_viral)
# report(random_search.grid_scores_)

############# Cross-validation of optimized model ##############
# print("Beginning cross-validation of the optimized model")
# clf = svm.SVC(kernel='rbf', C=10, gamma = 0.008, class_weight='auto', cache_size=4000)
# #clf = svm.LinearSVC(C = 4.6, class_weight= 'auto')
# scores = cross_validation.cross_val_score(clf, modscaled, cls_viral, cv=5, scoring='precision')
# print(scores)
# print("Precision: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
# scores = cross_validation.cross_val_score(clf, modscaled, cls_viral, cv=5, scoring='accuracy')
# print(scores)
# print("Accuracy: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
# scores = cross_validation.cross_val_score(clf, modscaled, cls_viral, cv=5, scoring='recall')
# print(scores)
# print("Recall: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))                        
# scores = cross_validation.cross_val_score(clf, modscaled, cls_viral, cv=5, scoring='f1')
# print(scores)
# print("f1: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2)) 


###############Generate ROC curve ###################
random_state = np.random.RandomState(0)

#shuffle and split training and test sets
n_samples, n_features = modscaled.shape
X, y = shuffle(modscaled, cls_viral, random_state=random_state)
half = int(n_samples / 2)
X_train, X_test = X[:half], X[half:]
y_train, y_test = y[:half], y[half:]


# #Run classifier
# classifier = svm.SVC(kernel='rbf',C=10,gamma=0.008, probability=True, class_weight='auto')
# probas_ = classifier.fit(X_train, y_train).predict_proba(X_test)
# 
# #Compute ROC curve and area the curve
# fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
# roc_auc = auc(fpr, tpr)
# print "Area under the ROC curve : %f" % roc_auc
# 
# #Plot ROC curve
# pl.clf()
# pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
# pl.plot([0, 1], [0, 1], 'k--')
# pl.xlim([0.0, 1.0])
# pl.ylim([0.0, 1.0])
# pl.xlabel('False Positive Rate')
# pl.ylabel('True Positive Rate')
# pl.title('Receiver operating characteristic, Viral vs. Non-viral')
# pl.legend(loc="lower right")
# pl.show()

#####Confusion matrix ##########
import matplotlib.pyplot as plt

from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix

X = modscaled
y = cls_viral
# Split the data into a training set and a test set
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

#Run classifier, using a model that is too regularized (C too low) to see
#the impact on the results
classifier = svm.SVC(kernel='rbf', C=10, gamma = 0.008, class_weight='auto')
y_pred = classifier.fit(X_train, y_train).predict(X_test)

print(classifier.get_params())
def plot_confusion_matrix(cm, title='Confusion matrix', cmap=plt.cm.Blues):
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    #tick_marks = np.arange(2)
    #plt.xticks(tick_marks, ["Non-viral", "Viral"], rotation=45)
    #plt.yticks(tick_marks, ["Non-viral", "Viral"])
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')


# Compute confusion matrix
cm = confusion_matrix(y_test, y_pred)
np.set_printoptions(precision=2)
print('Confusion matrix, without normalization')
print(cm)
plt.figure()
plot_confusion_matrix(cm)

# Normalize the confusion matrix by row (i.e by the number of samples
# in each class)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
print('Normalized confusion matrix')
print(cm_normalized)
plt.figure()
plot_confusion_matrix(cm_normalized, title='Normalized confusion matrix')

plt.show()
