#!/usr/bin/env python

from pyspark import SparkContext
from pyspark.mllib.util import MLUtils
from pyspark.mllib.classification import LogisticRegressionWithLBFGS, LogisticRegressionModel
from pyspark.mllib.linalg import Vectors
import argparse


def prediction(model_directory, libsvm_file, outputfile):
    sc = SparkContext(appName="PythonLinearRegressionWithSGDExample")

    model = LogisticRegressionModel.load(sc, model_directory)
    #print "numfeature",model.numFeatures
    #print "aaaaaaaa"
    vectors = MLUtils.loadLibSVMFile(sc, libsvm_file, numFeatures=model.numFeatures)
    vectors.cache()
    model.clearThreshold()
    # vector = vectors.collect()
    # for v in vector:
    #
    #     features = v.features
    #     print features
    #     print "bbbb",len(features),model.predict(Vectors.dense(features))
    # exit()
    scores = vectors.map(
        lambda p: (model.predict(Vectors.dense(p.features))))
     #   lambda p: (p.label, model.predict(p.features)))
    scores_list = scores.collect()
    file_out_obj = open(outputfile, 'w')
    for score in scores_list:
        #print '----->',score
        file_out_obj.write(str(score)+'\n')
    file_out_obj.close()



def main():
    parser = argparse.ArgumentParser(
        description='A script to use spark to predict label')

    parser.add_argument('libsvm', help='libsvm file')
    parser.add_argument('model',  help='model directory to use')
    parser.add_argument('outfile', help='output file')

    args = parser.parse_args()
    prediction(args.model, args.libsvm, args.outfile)

if __name__ == '__main__':
    main()
