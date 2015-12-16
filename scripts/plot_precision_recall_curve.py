#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
import glob
import os

dir = sys.argv[1]
switch = sys.argv[2] # 'single' or 'multiple'
title = sys.argv[3]

matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_np(file):
    file_in_obj = open(file,'r')
    precision = []
    recall = []
    for line in file_in_obj:
        line = line.rstrip()
        fields = line.split()
        precision.append(fields[0])
        recall.append(fields[1])
    p_np = np.array(precision)
    r_np = np.array(recall)
    return p_np, r_np
dir = sys.argv[1]

if switch == 'single':
    files = ['SVM_no_rescaling.model.p_r',  'SVM_rescaling.model.p_r',  'logistic_no_rescaling.model.p_r','logistic_rescaling.model.p_r']
else:
    files = ['logistic_no_rescaling.model.p_r','logistic_rescaling.model.p_r']

def plotting(file, c):
    print file
    p_np,r_np = load_np(file)
    
    idx = np.random.choice(np.arange(len(p_np)), 1000)
    plt.scatter(r_np[idx], p_np[idx], color=c, label=os.path.basename(file))

plotting(dir+'/'+files[0],'r')
plotting(dir+'/'+files[1],'b')
plotting(dir+'/'+files[2],'g')
plotting(dir+'/'+files[3],'y')

plt.legend(loc='lower left')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])

plt.title(title)
plt.savefig(dir+'/precision_recall.pdf')
