#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.mllib.util import MLUtils
from pyspark.mllib.feature import StandardScaler, StandardScalerModel
from pyspark.mllib.linalg import Vectors

from pyspark.mllib.classification import LogisticRegressionWithLBFGS, LogisticRegressionModel
import argparse


def training(model_directory, libsvm, scaler):
    sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")
    training_rdd = MLUtils.loadLibSVMFile(sc, libsvm)
    training_rdd.cache()
    if scaler == '1':
        label = training_rdd.map(lambda x: x.label)
        features = training_rdd.map(lambda x: x.features)

        scaler1 = StandardScaler().fit(features)
        data1 = label.zip(scaler1.transform(features))
        model_logistic = LogisticRegressionWithLBFGS.train(data1)
    else:
        model_logistic = LogisticRegressionWithLBFGS.train(training_rdd)
    model_logistic.save(sc, model_directory)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('model',  help='model directory to save')
    parser.add_argument('scaler',  help='standard scaler, 0 or 1')

    args = parser.parse_args()
    training(args.model, args.libsvm, args.scaler)

if __name__ == '__main__':
    main()
