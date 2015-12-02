#!/usr/bin/env python

import os
import subprocess
import argparse
import graphlab as gl

gl.set_runtime_config('GRAPHLAB_CACHE_FILE_LOCATIONS','/scratch')
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY',100000000000)
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY_PER_FILE', 100000000000)
gl.set_runtime_config('GRAPHLAB_SFRAME_SORT_BUFFER_SIZE',100000000000)

parser = argparse.ArgumentParser(description='A script to train and validate a model')
parser.add_argument('-i','--train', help='Input training matrix', required=True )
parser.add_argument('-t','--test', help = 'test data matrix', required=True)
parser.add_argument('-r,','--reportdir',help= 'Directory for reports of he model' ,  default='svmresults' )
args = parser.parse_args()

test =  gl.SFrame.read_csv(args.test, delimiter='\t', header=False)
train =  gl.SFrame.read_csv(args.train, delimiter='\t', header=False)

test.save('test_sframe')
train.save('train_sframe')

#model = gl.classifier.create(train, target='X1')

print "boosted_trees:\n"
model = gl.boosted_trees_classifier.create(train, target='X1', class_weights='auto', max_iterations= 15, verbose=True)
model.save("boosted_trees.model")

results = model.evaluate(test)
print(results)
print(results["confusion_matrix"].print_rows(200,3))


print "logistic:\n"
model = gl.logistic_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=True, max_iterations= 15, verbose=True)
model.save("logistic.model")

results = model.evaluate(test)
print(results)
print(results["confusion_matrix"].print_rows(200,3))

print "SVM:\n"
model = gl.svm_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=True, max_iterations= 15, verbose=True)
model.save("SVM.model")

results = model.evaluate(test)
print(results)
print(results["confusion_matrix"].print_rows(200,3))



