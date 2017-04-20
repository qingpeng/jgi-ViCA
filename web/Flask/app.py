from flask import Blueprint
main = Blueprint('main', __name__)

import json
from engine import GeneLearnEngine

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from flask import Flask, request


@main.route("/predict/<segment>", methods=["GET"])
def predict(segment):
  #  logger.debug("User %s rating requested for movie %s", user_id, movie_id)
    target = genelearn_engine.classification(segment)
    return json.dumps(target)




def create_app(spark_context):
    global genelearn_engine 

    genelearn_engine = GeneLearnEngine(spark_context)

    app = Flask(__name__)
    app.register_blueprint(main)
    return app