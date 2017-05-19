#!/usr/bin/env python

from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
import argparse
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix, precision_recall_curve
from sklearn import preprocessing


def evaluating(model, scaler_file, libsvm, report, scaler_with_mean):
    clf = joblib.load(model)
    num_features = len(clf.coef_[0])
    x_test, y_test = load_svmlight_file(libsvm, n_features=num_features)
    scaler = joblib.load(scaler_file)
    if scaler_with_mean == 'False':

        scaler.fit(x_test)
        x_test_scaled = scaler.transform(x_test)
    else:

        x_test_dense = x_test.todense()
        scaler.fit(x_test_dense)
        x_test_scaled = scaler.transform(x_test_dense)

    probability = clf.predict_proba(x_test_scaled)

    probability_list = [i[1] for i in probability]
    auprc = average_precision_score(y_test, probability_list)
    precision, recall, thresholds = precision_recall_curve(y_test,
                                                           probability_list)
    prediction = clf.predict(x_test_scaled)
    matrix = confusion_matrix(y_test, prediction)

    report_obj = open(report, 'w')
    report_obj.write('auprc='+str(auprc)+'\n')
    report_obj.write('confusion_matrix:\n')
    report_obj.write(str(matrix[0][0])+' '+str(matrix[0][1])+'\n')
    report_obj.write(str(matrix[1][0])+' '+str(matrix[1][1])+'\n')
    report_obj.write('precision recall thresholds\n')
    for i in range(len(precision)-1):
        report_obj.write(str(precision[i])+' '+str(recall[i])+' '
                         + str(thresholds[i])+'\n')
    report_obj.write(str(precision[len(precision)-1])+' '
                     + str(recall[len(precision)-1])+'\n')


def main():
    parser = argparse.ArgumentParser(
        description='A script to evaluate model with sklearn')

    parser.add_argument('libsvm', help='libsvm file of testing data')
    parser.add_argument('model',  help='model file name to load')
    parser.add_argument('scaler',  help='scaler file name to load')
    parser.add_argument('scaler_with_mean', help='True or False')
    parser.add_argument('report',  help='evaluation report file')

    args = parser.parse_args()
    evaluating(args.model, args.scaler, args.libsvm, args.report,
               args.scaler_with_mean)

if __name__ == '__main__':
    main()
