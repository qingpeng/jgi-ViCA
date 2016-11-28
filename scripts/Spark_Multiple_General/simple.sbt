name := "Simple Project"

version := "1.0"

scalaVersion := "2.10.5"


libraryDependencies ++= Seq(
  "org.apache.spark"  % "spark-core_2.10"              % "1.5.1",
  "org.apache.spark"  % "spark-mllib_2.10"             % "1.5.1"
  )
