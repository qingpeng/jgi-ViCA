/* 5k_5mer.scala */
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf

import org.apache.spark.mllib.classification.LogisticRegressionWithLBFGS
import org.apache.spark.mllib.evaluation.BinaryClassificationMetrics
import org.apache.spark.mllib.regression.LabeledPoint
import org.apache.spark.mllib.util.MLUtils
import org.apache.spark.mllib.classification.SVMModel

import org.apache.spark.mllib.classification.SVMWithSGD

import org.apache.spark.mllib.feature.StandardScaler
import org.apache.spark.mllib.linalg.Vectors



import java.io._

object Simple {
  
  def main(args: Array[String]) {
  
    if (args.length != 4) {
        println("Usage: SimpleApp <training_file> <testing_file> <output_dir> <algorithm>") 
        System.exit(1)
    }

    val training_file = args(0)
    val testing_file = args(1)
    val output_dir = args(2)
    val algorithm = args(3)
    
    
    val conf = new SparkConf().setAppName("Simple Application")
    val sc = new SparkContext(conf)
    
    
    val training = MLUtils.loadLibSVMFile(sc,training_file)
    val testing = MLUtils.loadLibSVMFile(sc, testing_file)
    training.cache()

// do scaling
    val scaler_training = new StandardScaler().fit(training.map(x => x.features))
    val training_scaled = training.map(x => LabeledPoint(x.label, scaler_training.transform(x.features)))
    val scaler_testing = new StandardScaler().fit(testing.map(x => x.features))
    val testing_scaled = testing.map(x => LabeledPoint(x.label, scaler_testing.transform(x.features)))
    val model = algorithm match {
        // training using logistic
        case "logistic" =>  new LogisticRegressionWithLBFGS().setNumClasses(2).run(training_scaled)
        case "svm" => {
            val numIterations = 100
            SVMWithSGD.train(training_scaled, numIterations).clearThreshold()
        }
        case "tree" => {
            val numIterations = 100
            SVMWithSGD.train(training_scaled, numIterations)
        }
    }
        
// testing the model
    val scoreAndLabels = testing_scaled.map { point =>
      val score = model.predict(point.features)
      (score, point.label)
      }
      
      
// evaluate the performance
    val metrics = new BinaryClassificationMetrics(scoreAndLabels)
    val auROC = metrics.areaUnderROC()
    val auPRC = metrics.areaUnderPR
    val PRC = metrics.pr
    
// output file into directory
    val file_report = output_dir + "/report.txt"
    val file_report_obj = new File(file_report)
    file_report_obj.getParentFile().mkdirs()
    val pw = new PrintWriter(file_report_obj)
    pw.write("training_file: %s \n".format(training_file))
    pw.write("testing_file: %s \n".format(testing_file))
    pw.write("algorithm: %s \n".format(algorithm))
    pw.write("auROC: %s \n".format(auROC))
    pw.write("auPRC: %s \n".format(auPRC))
    pw.close
    
    val dir_model = output_dir + "/model"
    model.save(sc, dir_model)
    
    val dir_PRC = output_dir + "/PRC"
    PRC.saveAsTextFile(dir_PRC)
    
    val dir_scoreAndLabels = output_dir + "/scoreAndLabels"
    scoreAndLabels.saveAsTextFile(dir_scoreAndLabels)



/*



    val numIterations = 100
    val model = SVMWithSGD.train(data3, numIterations)

    model.clearThreshold()
    val scoreAndLabels = data2.map { point =>
      val score = model.predict(point.features)
      (score, point.label)
      }
    val metrics = new BinaryClassificationMetrics(scoreAndLabels)
    val auROC = metrics.areaUnderROC()
    val auPRC = metrics.areaUnderPR


    val PRC = metrics.pr
    model.save(sc, "all_svm_model_log_scaling")
    PRC.saveAsTextFile("PRC_all_svm_log_scaling")
    scoreAndLabels.saveAsTextFile("scoreAndLabels_all_svm_log_scaling")


    val conf = new SparkConf().setAppName("Simple Application")
    val sc = new SparkContext(conf)

    val training = MLUtils.loadLibSVMFile(sc, training_file)
    val test = MLUtils.loadLibSVMFile(sc, testing_file)
    
    training.cache()
    
    val model = new LogisticRegressionWithLBFGS().setNumClasses(2).run(training)
    
    
    model.clearThreshold
    
    model.save(sc, "/global/projectb/scratch/qpzhang/Full_Training/Pfam/test_all_log_model")
    
    val predictionAndLabels = test.map { case LabeledPoint(label, features) =>
        val prediction = model.predict(features)
        (prediction, label)
    }
    predictionAndLabels.saveAsTextFile("/global/projectb/scratch/qpzhang/Full_Training/Pfam/test_all_log_predictionAndLabels")
    val metrics = new BinaryClassificationMetrics(predictionAndLabels)

    val PRC = metrics.pr

    val auPRC = metrics.areaUnderPR

    PRC.saveAsTextFile("/global/projectb/scratch/qpzhang/Full_Training/Pfam/test_all_log_PRC.txt")

    val pw = new PrintWriter(new File("/global/projectb/scratch/qpzhang/Full_Training/Pfam/test_all_log_auPRC.txt" ))
    pw.write("auPRC: %s".format(auPRC))
    pw.close
    

    val pw = new PrintWriter(new File(output))
    pw.write("training: %s".format(training))
    pw.write("testing: %s".format(testing))
    pw.close
*/

  }

  
}



