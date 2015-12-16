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
parser.add_argument('-r,','--reportdir',help= 'Directory to store related files for the model' ,  default='dato_model/' )
args = parser.parse_args()

dir = args.reportdir

test =  gl.SFrame.read_csv(args.test, delimiter='\t', header=False)
train =  gl.SFrame.read_csv(args.train, delimiter='\t', header=False)

#test_sframe_file = args.test+'.sframe'
#train_sframe_file = args.train+'.sframe'
#test.save(test_sframe_file)
#train.save(train_sframe_file)

#model = gl.classifier.create(train, target='X1')

def write_results(results, file_out, svm):
    file_out_fh = open(file_out,'w')
    if svm == True:
        keys = ['accuracy','precision','recall','f1_score']
    else:
        keys = ['accuracy','auc','precision','recall','f1_score','log_loss']
    for key in keys:
        file_out_fh.write(key+':'+str(results[key])+'\n')
    file_out_fh.close()
    

print "boosted_trees:\n"
model = gl.boosted_trees_classifier.create(train, target='X1', class_weights='auto', max_iterations= 50, verbose=True)
model.save(dir+"boosted_trees.model")
predictions = model.classify(test)
predictions.save(dir+"boosted_trees.predictions")
results = model.evaluate(test)
write_results(results,dir+"boosted_trees.results",False)
results['roc_curve'].save(dir+"boosted_trees.roc_curve")
print(results)
print(results["confusion_matrix"].print_rows(200,3))


print "logistic: feature_rescaling = True\n"
model = gl.logistic_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=True, max_iterations= 50, verbose=True)
model.save(dir+"logistic_rescaling.model")
predictions = model.classify(test)
predictions.save(dir+"logistic_rescaling.predictions")
results = model.evaluate(test)
write_results(results,dir+"logistic_rescaling.results",False)
results['roc_curve'].save(dir+"logistic_rescaling.roc_curve")
print(results)
print(results["confusion_matrix"].print_rows(200,3))

print "logistic: feature_rescaling = False\n"
model = gl.logistic_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=False, max_iterations= 50, verbose=True)
model.save(dir+"logistic_no_rescaling.model")
predictions = model.classify(test)
predictions.save(dir+"logistic_no_rescaling.predictions")
results = model.evaluate(test)
write_results(results,dir+"logistic_no_rescaling.results",False)
results['roc_curve'].save(dir+"logistic_no_rescaling.roc_curve")
print(results)
print(results["confusion_matrix"].print_rows(200,3))



print "SVM: feature_rescaling = True\n"
model = gl.svm_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=True, max_iterations= 50, verbose=True)
model.save(dir+"SVM_rescaling.model")
predictions = model.classify(test)
predictions.save(dir+"SVM_rescaling.predictions")
results = model.evaluate(test)
write_results(results,dir+"SVM_rescaling.results",True)
#results['roc_curve'].save(dir+"SVM_rescaling.roc_curve")
print(results)
print(results["confusion_matrix"].print_rows(200,3))

print "SVM: feature_rescaling = False\n"
model = gl.svm_classifier.create(train, target='X1', class_weights='auto', feature_rescaling=False, max_iterations= 50, verbose=True)
model.save(dir+"SVM_no_rescaling.model")
predictions = model.classify(test)
predictions.save(dir+"SVM_no_rescaling.predictions")
results = model.evaluate(test)
write_results(results,dir+"SVM_no_rescaling.results",True)
#results['roc_curve'].save(dir+"SVM_no_rescaling.roc_curve")
print(results)
print(results["confusion_matrix"].print_rows(200,3))


