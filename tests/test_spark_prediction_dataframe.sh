#!/usr/bin/env bash
~/Downloads/spark-2.1.0-bin-hadoop2.7/bin/spark-submit ../scripts/spark_prediction_dataframe.py test-data/testing.vect.200 training.vect.200.model training.vect.200.scaler test-data/testing.vect.200.prediction
