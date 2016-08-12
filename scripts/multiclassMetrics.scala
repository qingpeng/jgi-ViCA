import org.apache.spark.mllib.evaluation.MulticlassMetrics
import java.io._

val lines = sc.textFile("predictionAndLabels_family.txt.in_order")
val pattern = """\(([\d.]+),([\d.]+)\)""".r
val predictionAndLabels = lines.map{ line =>
	val pattern(s1,s2) = line
	(s1.toDouble,s2.toDouble)
}
val metrics = new MulticlassMetrics(predictionAndLabels)

val cMatrix = metrics.confusionMatrix

val file_matrix_obj = new File("matrix.txt")
val pw = new PrintWriter(file_matrix_obj)
pw.write(cMatrix.toString(811,6000))
pw.close()

val file_statistics_obj = new File("statistics.txt")
val pw = new PrintWriter(file_statistics_obj)
val precision = metrics.precision
val recall = metrics.recall // same as true positive rate
val f1Score = metrics.fMeasure
pw.write("Precision = %s\n".format(precision))
pw.write("Recall = %s\n".format(recall))
pw.write("F1 Score = %s\n".format(f1Score))

pw.write("Weighted precision = %s\n".format(metrics.weightedPrecision))
pw.write("Weighted recall = %s\n".format(metrics.weightedRecall))
pw.write("Weighted F1 score = %s\n".format(metrics.weightedFMeasure))
pw.write("Weighted false positive rate = %s\n".format(metrics.weightedFalsePositiveRate))

pw.close()

val file_precision_obj = new File("precision_by_label.txt")
val pw = new PrintWriter(file_precision_obj)
val labels = metrics.labels
labels.foreach { l =>
	pw.write("%s = %s\n".format(l,metrics.precision(l)))
}
pw.close()

val file_recall_obj = new File("recall_by_label.txt")
val pw = new PrintWriter(file_recall_obj)
val labels = metrics.labels
labels.foreach { l =>
	pw.write("%s = %s\n".format(l,metrics.recall(l)))
}
pw.close()

val file_fpr_obj = new File("fpr_by_label.txt")
val pw = new PrintWriter(file_fpr_obj)
val labels = metrics.labels
labels.foreach { l =>
	pw.write("%s = %s\n".format(l,metrics.falsePositiveRate(l)))
}
pw.close()

val file_f1_obj = new File("f1_by_label.txt")
val pw = new PrintWriter(file_f1_obj)
val labels = metrics.labels
labels.foreach { l =>
	pw.write("%s = %s\n".format(l,metrics.fMeasure(l)))
}
pw.close()


