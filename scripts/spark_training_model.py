#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.mllib.util import MLUtils
from pyspark.mllib.classification import LogisticRegressionWithLBFGS, LogisticRegressionModel
import argparse


def training(model_directory, libsvm):
    sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")
    training_rdd = MLUtils.loadLibSVMFile(sc, libsvm)
    training_rdd.cache()
    model_logistic = LogisticRegressionWithLBFGS.train(training_rdd)
    model_logistic.save(sc, model_directory)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('model',  help='model directory to save')

    args = parser.parse_args()
    training(args.model, args.libsvm)

if __name__ == '__main__':
    main()
