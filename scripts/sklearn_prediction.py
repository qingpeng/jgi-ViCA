#!/usr/bin/env python

from sklearn.datasets import load_svmlight_file
from sklearn.externals import joblib
import argparse
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix, precision_recall_curve


class Model:
    def __init__(self, model_path, scaler_path):
        self.model = joblib.load(model_path)
        self.scaler = joblib.load(scaler_path)

        self.prediction_list = []
        self.probability_list = []

    def predict(self, libsvm_path, scaler_with_mean):
        num_features = len(self.model.coef_[0])
        x_test, y_test = load_svmlight_file(libsvm_path, n_features=num_features)

        if scaler_with_mean == 'False':

            self.scaler.fit(x_test)
            x_test_scaled = self.scaler.transform(x_test)
        else:

            x_test_dense = x_test.todense()
            self.scaler.fit(x_test_dense)
            x_test_scaled = self.scaler.transform(x_test_dense)

        probability = self.model.predict_proba(x_test_scaled)
        self.prediction_list = self.model.predict(x_test_scaled)

        self.probability_list = [i[1] for i in probability]

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
    parser.add_argument('scaler_with_mean', help='True or False')
    parser.add_argument('outfile', help='output file')

    args = parser.parse_args()

    model = Model(args.model, args.scaler)
    print 'load model done!\n'
    model.predict(args.libsvm, args.scaler_with_mean)
    print 'prediction done!\n'
    model.write_output(args.outfile)
    print 'output result done!\n'

if __name__ == '__main__':
    main()

