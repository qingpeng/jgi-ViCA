#!/usr/bin/env python

from sklearn.datasets import load_svmlight_file
from sklearn import linear_model
from sklearn import preprocessing
from sklearn.externals import joblib
import argparse


class Model:
    def __init__(self):
        self.logreg = linear_model.LogisticRegression(C=0.0001,
                                                      class_weight='balanced')

    def train(self, libsvm, model, scaler_file, scaler_with_mean):
        x_train, y_train = load_svmlight_file(libsvm)
        if scaler_with_mean == 'False':
            scaler = preprocessing.StandardScaler(with_mean=False)
        else:
            scaler = preprocessing.StandardScaler(with_mean=True)
        scaler.fit(x_train)
        joblib.dump(scaler, scaler_file)
        scaled_x_train = scaler.transform(x_train)

        self.logreg.fit(scaled_x_train, y_train)
        joblib.dump(self.logreg, model)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use sklearn to train model')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('model',  help='model name to save')
    parser.add_argument('scaler',  help='standardscaler name to save')
    parser.add_argument('scaler_with_mean', help='True or False')

    args = parser.parse_args()
    model = Model()
    model.train(args.libsvm, args.model, args.scaler, args.scaler_with_mean)

if __name__ == '__main__':
    main()
