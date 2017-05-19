#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.ml.linalg import Vectors
from pyspark.ml.classification import LogisticRegression
from pyspark.sql.types import *
from pyspark.ml.feature import StandardScaler
from pyspark.sql import SparkSession

import argparse


def training_model(model_dir, libsvm, scaler_dir):

    spark = SparkSession \
        .builder \
        .appName("PythonSQL") \
        .config("spark.some.config.option", "some-value") \
        .getOrCreate()

    training = spark.read.format("libsvm").option(
        "numFeatures", 19423).load(libsvm)
    scaler = StandardScaler(inputCol="features", outputCol="scaledFeatures")
    scalerModel = scaler.fit(training)
    scaledData = scalerModel.transform(training)
    scalerModel.save(scaler_dir)

    lr_scaler = LogisticRegression(featuresCol="scaledFeatures")
    lrModel_scaler = lr_scaler.fit(scaledData)
    lrModel_scaler.save(model_dir)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('model_dir',  help='model directory to save')
    parser.add_argument('scaler_dir',  help='directory to save standardscaler')

    args = parser.parse_args()
    training_model(args.model_dir, args.libsvm, args.scaler_dir)

if __name__ == '__main__':
    main()
