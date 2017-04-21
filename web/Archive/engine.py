import os

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from pyspark.mllib.classification import LogisticRegressionWithLBFGS, LogisticRegressionModel
from pyspark.mllib.regression import LabeledPoint
from pyspark.mllib.util import MLUtils
from pyspark.mllib.linalg import SparseVector

class GeneLearnEngine:
    """genelearn engine
    """



    def __init__(self, sc):
        """Init the engine and train the model
        """

        logger.info("Starting up the GeneLearn Engine: ")

        self.sc = sc


        logger.info("Loading training data...")
        dataset_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts/Flask"
        training_file_path = os.path.join(dataset_path, 'training.svmlib')
        training = MLUtils.loadLibSVMFile(sc,training_file_path)
        self.model = LogisticRegressionWithLBFGS.train(training)



def classification(self, segment):
        """classify as virus or not
        """

        s_RDD = self.sc.parallelize([segment])
        scriptPath = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts/Flask/run_feature_extraction_for_spark.sh"
        result = s_RDD.pipe(scriptPath)
        vector = result.collect()
        
        index=[]
        value=[]
        for v in vector[0].split()[1:]:
            f = v.split(":")
            index.append(int(f[0]))
            value.append(float(f[1]))
        p = SparseVector(2104,index,value)
        
        target = self.model.predict(p)

        return target

    # Attach the functions to class methods
GeneLearnEngine.classification = classification


