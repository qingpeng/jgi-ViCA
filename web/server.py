from flask import Flask, redirect, url_for, request
import subprocess
app = Flask(__name__)


@app.route('/')
def main_form():
    return '<form action="submit" id="textform" method="post"><p>Enter the fasta sequences:</p><textarea name="text" cols="40" rows="5" maxlength="5000000"></textarea><input type="submit" value="Submit"></form>'


@app.route('/submit', methods=['POST'])
def submit_textarea():
    fasta_seqs = request.form["text"]
    file_fasta_obj = open("file_fasta.fa", 'w')
    file_fasta_obj.write(fasta_seqs)
    file_fasta_obj.close()
    genelearn_path = "/Users/qingpeng/Dropbox/Development/Github/jgi-ViCA/scripts/"
    genemark_path = "/Users/qingpeng/bin/genemark_suite_macosx/gmsuite/"
    hmmer_path = "/Users/qingpeng/bin/hmmer-3.1b2-macosx-intel/"
    pfam_db = "/Users/qingpeng/Local/Pfam_DB/"
    spark_path = "/Users/qingpeng/Downloads/spark-2.1.0-bin-hadoop2.7/"
    feature_index_file = "/Users/qingpeng/Local/GeneLearn/all_segment.fasta.vect.feature_index"
    model_path = "/Users/qingpeng/Dropbox/Development/Github/jgi-ViCA/scripts/model/subsample_model/"
    prediction_pipeline_lite_command = [
        'python',
        genelearn_path+'prediction_pipeline_lite.py', "file_fasta.fa",
        "file_fasta_fa.prediction",
        genemark_path, hmmer_path, pfam_db, spark_path, feature_index_file,
        model_path]
    print prediction_pipeline_lite_command
    print "spark prediction running...\n"
    return_code = subprocess.call(prediction_pipeline_lite_command)

    if return_code != 0:
        return return_code, "prediction_pipeline_lite_command"
    print "done!\n"

    file_prediction_obj = open("file_fasta_fa.prediction", 'r')
    result = '<br>'.join(file_prediction_obj.readlines())

    return "prediction score:<br> {}".format(result)


if __name__ == '__main__':
   app.run(debug = True)