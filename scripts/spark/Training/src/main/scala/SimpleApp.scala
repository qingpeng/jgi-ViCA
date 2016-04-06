/* model_building.scala */
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

import org.apache.spark.mllib.tree.RandomForest
import org.apache.spark.mllib.tree.model.RandomForestModel
import org.apache.spark.mllib.tree.DecisionTree
import org.apache.spark.mllib.tree.model.DecisionTreeModel

import org.apache.spark.mllib.tree.GradientBoostedTrees
import org.apache.spark.mllib.tree.configuration.BoostingStrategy
import org.apache.spark.mllib.tree.model.GradientBoostedTreesModel


import java.io._

object Simple {

    def testing_model(model, testing_scaled, output_dir, algorithm) = {

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
    }

    def testing_GBT_model(model, testing_scaled, output_dir, algorithm) = {
    // Evaluate model on test instances and compute test error
        val labelAndPreds = testing_scaled.map { point =>
        val prediction = model.predict(point.features)
        (point.label, prediction)
        }
        val testErr = labelAndPreds.filter(r => r._1 != r._2).count.toDouble / testing_scaled.count()

    // output file into directory
        val file_report = output_dir + "/report.txt"
        val file_report_obj = new File(file_report)
        file_report_obj.getParentFile().mkdirs()
        val pw = new PrintWriter(file_report_obj)
        pw.write("training_file: %s \n".format(training_file))
        pw.write("testing_file: %s \n".format(testing_file))
        pw.write("algorithm: %s \n".format(algorithm))
        pw.write("numInteration: %s \n".format(numIteration))
        pw.write("Test Error: %s \n".format(testErr))
        pw.write("Learned classification GBT model: %s \n".format(model.toDebugString))
    }

        
    def logistic(parameters:String, output_dir:String, training_scaled, testing_scaled, algorithm) = {
            val model = new LogisticRegressionWithLBFGS().setNumClasses(2).run(training_scaled)
            testing_model(model, testing_scaled, output_dir, algorithm)
    }
        
    def SVM(parameters:String, output_dir:String, training_scaled, testing_scaled, algorithm) = {
            val parameters_array = parameters.split(" = ")
            val numIterations = parameters_array(1) 
            val model = SVMWithSGD.train(training_scaled, numIterations).clearThreshold()
            testing_model(model, testing_scaled, output_dir, algorithm)
    }
    
    def RFT(parameters:String, output_dir:String, training_scaled, testing_scaled, algorithm) = {
            val parameters_array = parameters.split(" = ")
            val numClasses = 2
            val categoricalFeaturesInfo = Map[Int, Int]()
            val numTrees = parameters_array(1)
            val featureSubsetStrategy = "auto" // Let the algorithm choose.
            val impurity = "gini"
            val maxDepth = 4
            val maxBins = 32
            val model = RandomForest.trainClassifier(training_scaled, numClasses, categoricalFeaturesInfo, numTrees, featureSubsetStrategy, impurity, maxDepth, maxBins) 
            testing_model(model, testing_scaled, output_dir, algorithm)
    }
                
    def GBT(parameters:String, output_dir:String, training_scaled, testing_scaled, algorithm) = {
            val boostingStrategy = BoostingStrategy.defaultParams("Classification")
            boostingStrategy.numIterations = numIteration // Note: Use more iterations in practice.
            boostingStrategy.treeStrategy.numClasses = 2
            boostingStrategy.treeStrategy.maxDepth = 5
        //  Empty categoricalFeaturesInfo indicates all features are continuous.
            boostingStrategy.treeStrategy.categoricalFeaturesInfo = Map[Int, Int]()

            val model = GradientBoostedTrees.train(training_scaled, boostingStrategy)  
            testing_GBT_model(model, testing_scaled, output_dir, algorithm)
    }
                
  def main(args: Array[String]) {
/* parameters support:
    logistic: "logistic", "None"
    SVM: "SVM", "numIterations = 100"
    Random Forest: "RFT", "numTrees = 3"
    GBT: "GBT", "numIteration = 30"
*/
  
    if (args.length != 5) {
        println("Usage: model_building <training_file> <testing_file> <output_dir> <algorithm> <parameters>") 
        System.exit(1)
    }

    val training_file = args(0)
    val testing_file = args(1)
    val output_dir = args(2)
    val algorithm = args(3)
    val parameters = args(4)
    
    
    val conf = new SparkConf().setAppName("Building model and testing")
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
        case "logistic" =>  logistic(parameters, output_dir, training_scaled, testing_scaled, algorithm)

        case "SVM" => SVM(parameters, output_dir, training_scaled, testing_scaled, algorithm)
        case "RFT" => RFT(parameters, output_dir, training_scaled, testing_scaled, algorithm)
        case "GBT" => GBT(parameters, output_dir, training_scaled, testing_scaled, algorithm)
    }
       


  }

  
}



