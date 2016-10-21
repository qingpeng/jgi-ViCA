import org.apache.spark.mllib.evaluation.MulticlassMetrics
import org.apache.spark.mllib.regression.LabeledPoint
import org.apache.spark.mllib.linalg.Vectors
import org.apache.spark.mllib.util.MLUtils
import org.apache.spark.mllib.feature.StandardScaler


import org.apache.spark.mllib.classification.{LogisticRegressionWithLBFGS, LogisticRegressionModel}


import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf


object Simple {

    def main(args: Array[String]) {
        val conf = new SparkConf().setAppName("Building model and testing")
        val sc = new SparkContext(conf)
    


        val training_5k = MLUtils.loadLibSVMFile(sc,"/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam_run/training.vect.svmlib")
        val test_5k = MLUtils.loadLibSVMFile(sc, "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam_run/testing.vect.svmlib")
        training_5k.cache()
        val scaler1 = new StandardScaler().fit(training_5k.map(x => x.features))
        val data3 = training_5k.map(x => LabeledPoint(x.label, scaler1.transform(x.features)))
        val scaler2 = new StandardScaler().fit(test_5k.map(x => x.features))
        val data2 = test_5k.map(x => LabeledPoint(x.label, scaler2.transform(x.features)))
        val model = new LogisticRegressionWithLBFGS().setNumClasses(2).run(data3)
        val predictionAndLabels = data2.map { case LabeledPoint(label, features) =>
                val prediction = model.predict(features)
                (prediction, label)
            }

        predictionAndLabels.saveAsTextFile("predictionAndLabels_binary")
        model.save(sc, "predictionAndLabels_logistics_model_binary")


    }
}
