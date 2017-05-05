#!/usr/bin/env python

from sklearn.datasets import load_svmlight_file
from sklearn import linear_model, datasets

from sklearn.externals import joblib
import argparse


def training(model_dir, libsvm):
    X_train, y_train = load_svmlight_file(libsvm)
    logreg = linear_model.LogisticRegression(C=0.0001)
    logreg.fit(X_train, y_train)
    joblib.dump(logreg, model_dir)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('model_dir',  help='model directory to save')
#    parser.add_argument('scaler_dir',  help='directory to save standardscaler')

    args = parser.parse_args()
    training(args.model_dir, args.libsvm)

if __name__ == '__main__':
    main()
