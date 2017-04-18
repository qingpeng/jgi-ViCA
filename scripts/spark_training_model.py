#!/usr/bin/env python

from pyspark import SparkContext

from pyspark.mllib.util import MLUtils
from pyspark.mllib.classification import LogisticRegressionWithLBFGS, LogisticRegressionModel
from pyspark.mllib.regression import LabeledPoint
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

def evaluate_model(testing, model):
    model.clearThreshold()
    labelsAndScores = testing.map(lambda p: (p.label, model.predict(p.features)))
    LS = labelsAndScores.collect()
    y_true = []
    y_scores = []
    for s in LS:
        y_true.append(s[0])
        y_scores.append(s[1])
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    au_PRC = average_precision_score(y_true, y_scores)
    return au_PRC, precision, recall, thresholds



sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")
# training = MLUtils.loadLibSVMFile(sc, "training.vect")
# training.cache()
testing = MLUtils.loadLibSVMFile(sc, "testing.vect",19421)
testing.cache()
# model_logistic = LogisticRegressionWithLBFGS.train(training)
# model_logistic.save(sc, "./model_logistic")
sameModel = LogisticRegressionModel.load(sc, "./model_logistic")
result = evaluate_model(testing, sameModel)
print result
# file_out_obj = open("spark_evaluate_result.out",'w')
# file_out_obj.write(result[0]+'\n')

