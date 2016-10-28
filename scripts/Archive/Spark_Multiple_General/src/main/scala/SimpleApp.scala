import org.apache.spark.mllib.evaluation.MulticlassMetrics
import org.apache.spark.mllib.regression.LabeledPoint
import org.apache.spark.mllib.linalg.Vectors
import org.apache.spark.mllib.util.MLUtils
import org.apache.spark.mllib.feature.StandardScaler


import org.apache.spark.mllib.classification.{LogisticRegressionWithLBFGS, LogisticRegressionModel}


import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import java.io._

object Simple {

    def main(args: Array[String]) {


        if (args.length != 4) {
            println("Usage: multiple_class <training_file> <testing_file> <output_dir> <class_num>") 
            System.exit(1)
        }

        val training_file = args(0)
        val testing_file = args(1)
        val output_dir = args(2)
        val class_num = args(3).toInt

    


        val conf = new SparkConf().setAppName("Building model and testing")
        val sc = new SparkContext(conf)
    


        val training_5k = MLUtils.loadLibSVMFile(sc,training_file)
        val test_5k = MLUtils.loadLibSVMFile(sc, testing_file)
        training_5k.cache()
        val scaler1 = new StandardScaler().fit(training_5k.map(x => x.features))
        val data3 = training_5k.map(x => LabeledPoint(x.label, scaler1.transform(x.features)))
        data3.cache()
        val scaler2 = new StandardScaler().fit(test_5k.map(x => x.features))
        val data2 = test_5k.map(x => LabeledPoint(x.label, scaler2.transform(x.features)))
        data2.cache()
        val model = new LogisticRegressionWithLBFGS().setNumClasses(class_num).run(data3)
        val predictionAndLabels = data2.map { case LabeledPoint(label, features) =>
                val prediction = model.predict(features)
                (prediction, label)
            }

        val dir_predictionAndLabels = output_dir + "/predictionAndLabels"
        predictionAndLabels.saveAsTextFile(dir_predictionAndLabels)

        val metrics = new MulticlassMetrics(predictionAndLabels)

        val cMatrix = metrics.confusionMatrix

        val file_matrix_obj = new File(output_dir + "matrix.txt")
        var pw = new PrintWriter(file_matrix_obj)
        pw.write(cMatrix.toString(class_num,class_num * 10))
        pw.close()

        val file_statistics_obj = new File(output_dir + "statistics.txt")
        pw = new PrintWriter(file_statistics_obj)
        val precision = metrics.precision
        val recall = metrics.recall // same as true positive rate
        val f1Score = metrics.fMeasure

        pw.write("training_file: %s \n".format(training_file))
        pw.write("testing_file: %s \n".format(testing_file))
        pw.write("class_num: %s \n".format(class_num))

        pw.write("Precision = %s\n".format(precision))
        pw.write("Recall = %s\n".format(recall))
        pw.write("F1 Score = %s\n".format(f1Score))

        pw.write("Weighted precision = %s\n".format(metrics.weightedPrecision))
        pw.write("Weighted recall = %s\n".format(metrics.weightedRecall))
        pw.write("Weighted F1 score = %s\n".format(metrics.weightedFMeasure))
        pw.write("Weighted false positive rate = %s\n".format(metrics.weightedFalsePositiveRate))

        pw.close()

        val file_precision_obj = new File(output_dir + "precision_by_label.txt")
        pw = new PrintWriter(file_precision_obj)
        val labels = metrics.labels
        labels.foreach { l =>
          pw.write("%s = %s\n".format(l,metrics.precision(l)))
        }
        pw.close()

        val file_recall_obj = new File(output_dir + "recall_by_label.txt")
        pw = new PrintWriter(file_recall_obj)
        labels.foreach { l =>
          pw.write("%s = %s\n".format(l,metrics.recall(l)))
        }
        pw.close()

        val file_fpr_obj = new File(output_dir + "fpr_by_label.txt")
        pw = new PrintWriter(file_fpr_obj)
        labels.foreach { l =>
          pw.write("%s = %s\n".format(l,metrics.falsePositiveRate(l)))
        }
        pw.close()

        val file_f1_obj = new File(output_dir + "f1_by_label.txt")
        pw = new PrintWriter(file_f1_obj)
        labels.foreach { l =>
          pw.write("%s = %s\n".format(l,metrics.fMeasure(l)))
        }
        pw.close()


    }
}
