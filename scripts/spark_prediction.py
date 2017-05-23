#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.ml.classification import LogisticRegressionModel
from pyspark.ml.feature import StandardScaler, StandardScalerModel
from pyspark.sql import SparkSession
import argparse


class Model:
    def __init__(self, model_path, scaler_path):
        self.model = LogisticRegressionModel.load(model_path)
        self.scaler = StandardScalerModel.load(scaler_path)
        self.prediction_list = []
        self.probability_list = []

    def predict(self, libsvm_path):
        spark = SparkSession \
            .builder \
            .appName("LogisticRegressionWithElasticNet") \
            .getOrCreate()

        testing = spark.read.format(
            "libsvm").option("numFeatures", self.model.numFeatures).load(
            libsvm_path)
        testing_scaler = self.scaler.transform(testing)
        testing_scaler_prediction = self.model.transform(testing_scaler)
        # testing_scaler_prediction.show()
        self.prediction_list = [i.prediction for i in
                                testing_scaler_prediction.select('prediction')
                                .collect()]

        self.probability_list = [i.probability[1] for i in
                                 testing_scaler_prediction.select(
                                'probability').collect()]
        # print self.probability_list

    def write_output(self, output_path):
        file_output_obj = open(output_path, 'w')
        for i in range(len(self.probability_list)):
            file_output_obj.write(str(self.probability_list[i])+' ' +
                                  str(self.prediction_list[i])+'\n')
        file_output_obj.close()


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to predict label')

    parser.add_argument('libsvm', help='libsvm file')
    parser.add_argument('model',  help='model directory to use')
    parser.add_argument('scaler',  help='scaler file name to load')
    parser.add_argument('outfile', help='output file')

    args = parser.parse_args()
    sc = SparkContext(appName="prediction")

    model = Model(args.model, args.scaler)
    print 'load model done!\n'
    model.predict(args.libsvm)
    print 'prediction done!\n'
    model.write_output(args.outfile)
    print 'output result done!\n'

if __name__ == '__main__':
    main()
