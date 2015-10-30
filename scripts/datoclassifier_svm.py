#!/usr/bin/env python

import os
import subprocess
import argparse
import graphlab as gl

parser = argparse.ArgumentParser(description='A script to train and validate a model')
parser.add_argument('-i','--train', help='Input training matrix', required=True )
parser.add_argument('-t','--test', help = 'test data matrix', required=True)
parser.add_argument('-d,','--modeldir',help= 'Directory to save the model' ,  default='svmmodel' )
parser.add_argument('-r,','--report',help= 'report file' ,  default='report.txt' )
args = parser.parse_args()

gl.set_runtime_config('GRAPHLAB_CACHE_FILE_LOCATIONS','/scratch')
test =  gl.SFrame.read_csv(args.test, delimiter='\t', header=False)
train =  gl.SFrame.read_csv(args.train, delimiter='\t', header=False)
file_report = open(args.report, 'w')


test.save('test_sframe')
train.save('train_sframe')

model = gl.svm_classifier.create(train, target='X1')
predictions = model.predict(test)
print predictions
#file_report.write(predictions)
results = model.evaluate(test)
print results
#file_report.write(results)
model.save(args.modeldir)
print model.summary
#file_report.write(model.summary)
