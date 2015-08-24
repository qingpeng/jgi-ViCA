#!/usr/bin/env python

import os
import subprocess
import argparse
import graphlab as gl

parser = argparse.ArgumentParser(description='A script to train and validate a model')
parser.add_argument('-i','--train', help='Input training matrix', required=True )
parser.add_argument('-t','--test', help = 'test data matrix', required=True)
parser.add_argument('-r,','--reportdir',help= 'Directory for reports of he model' ,  default='svmresults' )
args = parser.parse_args()

gl.set_runtime_config('GRAPHLAB_CACHE_FILE_LOCATIONS','/scratch')
test =  gl.SFrame.read_csv(args.test, delimiter='\t', header=False)
train =  gl.SFrame.read_csv(args.train, delimiter='\t', header=False)

test.save('test_sframe')
train.save('train_sframe')

model = gl.classifier.create(train, target='X1')
predictions = model.classify(test)
print(predictions)
results = model.evaluate(test)
print(results)
model.save("dato.model")
print(model.summary)
