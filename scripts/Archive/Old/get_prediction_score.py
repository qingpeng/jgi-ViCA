#!/usr/bin/env python

import graphlab as gl


import argparse

gl.set_runtime_config('GRAPHLAB_CACHE_FILE_LOCATIONS','/scratch')
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY',40000000000)
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY_PER_FILE', 40000000000)
gl.set_runtime_config('GRAPHLAB_SFRAME_SORT_BUFFER_SIZE',40000000000)

parser = argparse.ArgumentParser(description='A script to get P-R curve for a model')
parser.add_argument('-m','--model', help='Model to test', required=True )
parser.add_argument('-t','--test', help = 'test data matrix', required=True)
parser.add_argument('-r,','--report',help= 'predition with score', required = True )
#parser.add_argument('-f,','--figure',help= 'figure name for plotting', required = True )

args = parser.parse_args()


#test =  gl.SFrame.read_csv('/global/projectb/scratch/arrivers/geneleanrntest/20150818/test.twoclass.txt', delimiter='\t', header=False)
#train =  gl.SFrame.read_csv('/global/projectb/scratch/arrivers/geneleanrntest/20150818/train.twoclass.txt', delimiter='\t', header=False)
#test.save('test_twoclass_sframe')
#train.save('train_twoclass_sframe')

model = gl.load_model(args.model)
test =  gl.SFrame.read_csv(args.test, delimiter='\t', header=False)


#test = gl.load_sframe('../test_Both_1.vect.pro.convert.sframe')
#train = gl.load_sframe('train_twoclass_sframe')
#model = gl.svm_classifier.create(train, target='X1', class_weights='auto', max_iterations=50)
#model.save("dato.SVM.model")
#model = gl.load_model('dato.SVM.model')
#results = model.evaluate(test)
#pc = model.predict(test)
pm = model.predict(test, output_type='margin',missing_value_action='auto')


#print precision
#print recall
#print thresholds

test_in = open(args.test,'r')

file_report_obj = open(args.report,'w')
i = 0
for line in test_in:
    line = line.rstrip()
    fields = line.split()
    file_report_obj.write(fields[0]+' '+str(pm[i])+'\n')
    i = i+1

# 
