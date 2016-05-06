#!/usr/bin/env python

import graphlab as gl
from sklearn.metrics import precision_recall_curve
import numpy as np
#import matplotlib.pyplot as plt
import argparse

gl.set_runtime_config('GRAPHLAB_CACHE_FILE_LOCATIONS','/scratch')
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY',40000000000)
gl.set_runtime_config('GRAPHLAB_FILEIO_MAXIMUM_CACHE_CAPACITY_PER_FILE', 40000000000)
gl.set_runtime_config('GRAPHLAB_SFRAME_SORT_BUFFER_SIZE',40000000000)

parser = argparse.ArgumentParser(description='A script to get P-R curve for a model')
parser.add_argument('-m','--model', help='Model to test', required=True )
parser.add_argument('-t','--test', help = 'test data matrix', required=True)
parser.add_argument('-r,','--report',help= 'data for PR curve', required = True )
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

x1 = np.array(test['X1'])
#print x1[:100]

pm_np = np.array(pm)
#print pm_np[:100]
precision, recall, thresholds = precision_recall_curve(x1,pm_np, pos_label="virus")

#print precision
#print recall
#print thresholds


file_report_obj = open(args.report,'w')

for k in range(len(precision)-1):
    file_report_obj.write(str(precision[k])+'\t'+str(recall[k])+'\t'+str(thresholds[k])+'\n')

file_report_obj.close()
# 
# #print(results)
# #print(results["confusion_matrix"].print_rows(200,3))
# 
# #prdata = gl.SFrame(test['X1'], gl.SArray(pm))
# # Plot Precision-Recall curve
# #plt.clf()
# plt.plot(recall[0], precision[0], label='Precision-Recall curve')
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.ylim([0.0, 1.05])
# plt.xlim([0.0, 1.0])
# plt.title('Precision-Recall example: AUC={0:0.2f}'.format(average_precision[0]))
# plt.legend(loc="lower left")
# plt.savefig(args.figure)
# 
