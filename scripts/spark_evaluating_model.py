#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.ml.classification import LogisticRegressionModel
from pyspark.ml.feature import StandardScaler, StandardScalerModel
from pyspark.sql import SparkSession
import argparse
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt


class Model:
    def __init__(self, model_path, scaler_path):

        self.model = LogisticRegressionModel.load(model_path)
        self.scaler = StandardScalerModel.load(scaler_path)
        self.au_PRC = 0
        self.precision = []
        self.recall = []
        self.thresholds = []
        self.matrix = []

    def evaluate(self, test_path):
        spark = SparkSession \
            .builder \
            .appName("LogisticRegressionWithElasticNet") \
            .getOrCreate()

        testing = spark.read.format(
            "libsvm").option("numFeatures", self.model.numFeatures).load(
            test_path)
        testing_scaler = self.scaler.transform(testing)
        testing_scaler_prediction = self.model.transform(testing_scaler)
        label_list = [i.label for i in
                      testing_scaler_prediction.select('label').collect()]
        prediction_list = [i.prediction for i in
                           testing_scaler_prediction.select(
                               'prediction').collect()]
        probability_list = [i.probability[1] for i in
                            testing_scaler_prediction.select(
                                'probability').collect()]
        self.au_PRC = average_precision_score(label_list, probability_list)
        self.precision, self.recall, self.thresholds = precision_recall_curve(
            label_list, probability_list)
        self.matrix = confusion_matrix(label_list, prediction_list)

    def write_to_report(self, report_path):
        report_obj = open(report_path, 'w')
        report_obj.write('auprc=' + str(self.au_PRC) + '\n')
        report_obj.write('confusion_matrix:\n')
        report_obj.write(str(self.matrix[0][0]) + ' ' + str(self.matrix[0][1])
                         + '\n')
        report_obj.write(str(self.matrix[1][0]) + ' ' + str(self.matrix[1][1])
                         + '\n')
        report_obj.write('precision recall thresholds\n')
        for i in range(len(self.precision) - 1):
            report_obj.write(str(self.precision[i]) + ' ' + str(self.recall[i])
                             + ' '
                             + str(self.thresholds[i]) + '\n')
        report_obj.write(str(self.precision[len(self.precision) - 1]) + ' '
                         + str(self.recall[len(self.precision) - 1]) + '\n')

    def draw_prc(self, figure_path):
        lw = 2
        plt.figure(figsize=(6, 4), dpi=80)
        plt.clf()
        plt.plot(self.recall, self.precision, lw=lw, color='teal',
                 label='all, AUC={0:0.4f}'.format(self.au_PRC))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve')
        plt.legend(loc="lower left")
        plt.savefig(figure_path)


def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to train model')

    parser.add_argument('libsvm', help='libsvm file of testing data')
    parser.add_argument('model',  help='model directory to save')
    parser.add_argument('scaler',  help='scaler file name to load')
    parser.add_argument('report', help='report file')
    parser.add_argument('PRC', help='PRC figure file')

    args = parser.parse_args()
    sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")

    model = Model(args.model, args.scaler)
    #print 'load model done!\n'
    model.evaluate(args.libsvm)
    #print 'evaluating_model done!\n'
    model.write_to_report(args.report)
    #print 'write_to_report done!\n'
    model.draw_prc(args.PRC)
    #print 'draw_prc done!\n'

if __name__ == '__main__':
    main()
