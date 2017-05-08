#!/usr/bin/env python

from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
import argparse
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix, precision_recall_curve


def evaluating(model, libsvm, report):
    x_test, y_test = load_svmlight_file(libsvm)
    clf = joblib.load(model)

    probability = clf.predict_proba(x_test)
    probability_list = [i[1] for i in probability]
    auprc = average_precision_score(y_test, probability_list)
    precision, recall, thresholds = precision_recall_curve(y_test,
                                                           probability_list)
    prediction = clf.predict(x_test)
    matrix = confusion_matrix(y_test, prediction)

    report_obj = open(report, 'w')
    report_obj.write('auprc='+str(auprc)+'\n')
    report_obj.write('confusion_matrix:\n')
    report_obj.write(str(matrix[0][0])+' '+str(matrix[0][1])+'\n')
    report_obj.write(str(matrix[1][0])+' '+str(matrix[1][1])+'\n')
    report_obj.write('precision recall thresholds\n')
    for i in range(len(precision)):
        report_obj.write(precision[i]+' '+recall[i]+' '+thresholds[i]+'\n')


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of testing data')
    parser.add_argument('model',  help='model name to save')
    parser.add_argument('report',  help='evaluation report file')

    args = parser.parse_args()
    evaluating(args.model, args.libsvm, args.report)

if __name__ == '__main__':
    main()
