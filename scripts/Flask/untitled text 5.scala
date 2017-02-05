
        val training_5k = MLUtils.loadLibSVMFile(sc,"/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/training.svmlib")
        val test_5k = MLUtils.loadLibSVMFile(sc, "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/testing.svmlib")
        training_5k.cache()
        test_5k.cache()
        val scaler1 = new StandardScaler().fit(training_5k.map(x => x.features))
        val data3 = training_5k.map(x => LabeledPoint(x.label, scaler1.transform(x.features)))
        val scaler2 = new StandardScaler().fit(test_5k.map(x => x.features))
        val data2 = test_5k.map(x => LabeledPoint(x.label, scaler2.transform(x.features)))
        val model = new LogisticRegressionWithLBFGS().setNumClasses(2).run(data3)
        val predictionAndLabels = data2.map { case LabeledPoint(label, features) =>
                val prediction = model.predict(features)
                (prediction, label)
            }

        predictionAndLabels.saveAsTextFile("/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/predictionAndLabels_test")
        model.save(sc, "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/predictionAndLabels_logistics_model_test")

        sc.stop()
    }
}
