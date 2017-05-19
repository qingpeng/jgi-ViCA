#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.mllib.util import MLUtils
from pyspark.mllib.classification \
    import LogisticRegressionWithLBFGS, LogisticRegressionModel
import argparse
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt


def evaluate_model(testing, model):
    model.clearThreshold()
    labels_and_scores_rdd = testing.map(
        lambda p: (p.label, model.predict(p.features)))
    labels_and_scores = labels_and_scores_rdd.collect()
    y_true = []
    y_scores = []
    for s in labels_and_scores:
        y_true.append(s[0])
        y_scores.append(s[1])
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    au_prc = average_precision_score(y_true, y_scores)
    return au_prc, precision, recall, thresholds, y_true, y_scores


def write_to_report(au_prc, precision, recall, thresholds, report):
    file_report_obj = open(report, 'w')
    file_report_obj.write('#au_prc:\n')
    file_report_obj.write(str(au_prc)+'\n')
    file_report_obj.write('#precision recall thresholds\n')
    # print precision
    # print 'length='+str(len(precision))+'\n'
    # print recall
    # print 'length='+str(len(recall))+'\n'
    # print thresholds
    # print 'length='+str(len(thresholds))+'\n'

    for i in xrange(len(precision)-1):
        file_report_obj.write(
            str(precision[i])+' '+str(recall[i])+' '+str(thresholds[i])+'\n')
    file_report_obj.write(
            str(precision[-1])+' '+str(recall[-1])+' '+'N/A'+'\n')
    file_report_obj.close()


def write_to_prediction(y_true, y_scores, prediction):
    file_prediction_obj = open(prediction, 'w')
    file_prediction_obj.write('#true_label score\n')

    for i in xrange(len(y_true)):
        file_prediction_obj.write(str(y_true[i])+' '+str(y_scores[i])+'\n')
    file_prediction_obj.close()


def draw_prc(precision, recall, prc_file, au_prc):
    lw = 2
    plt.figure(figsize=(6, 4), dpi=80)
    plt.clf()
    plt.plot(recall, precision, lw=lw, color='teal',
             label='all, AUC={0:0.4f}'.format(au_prc))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve')
    plt.legend(loc="lower left")
    plt.savefig(prc_file)


def testing_model(model_directory, libsvm, prediction, report, prc_file):
    sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")
    model = LogisticRegressionModel.load(sc, model_directory)
    testing_rdd = MLUtils.loadLibSVMFile(sc, libsvm,
                                         numFeatures=model.numFeatures)
    testing_rdd.cache()
    au_prc, precision, recall, thresholds, y_true, y_scores = evaluate_model(
        testing_rdd, model)
    print 'evaluating_model done!\n'
    write_to_report(au_prc, precision, recall, thresholds, report)
    print 'write_to_report done!\n'
    write_to_prediction(y_true, y_scores, prediction)
    print 'write_to_prediction done!\n'
    draw_prc(precision, recall, prc_file, au_prc)
    print 'draw_prc done!\n'

def evaluation_model(model_dir, libsvm, scaler_dir):
    spark = SparkSession \
        .builder \
        .appName("PythonSQL") \
        .config("spark.some.config.option", "some-value") \
        .getOrCreate()

    testing = spark.read.format("libsvm").load(libsvm)
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

    parser.add_argument('libsvm', help='libsvm file of testing data')
    parser.add_argument('model',  help='model directory to save')
    parser.add_argument('prediction',
                        help='prediction result, label and score')
    parser.add_argument('report', help='report file')
    parser.add_argument('PRC', help='PRC figure file')

    args = parser.parse_args()
    testing_model(args.model, args.libsvm, args.prediction, args.report,
                  args.PRC)

if __name__ == '__main__':
    main()
