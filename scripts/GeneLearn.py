#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
from Bio import SeqIO



parser = argparse.ArgumentParser(description='A script to create relative scripts to reflect current environment variables' )
parser.add_argument('-c', '--config', help ="json formatted config file", required= True)
args = parser.parse_args()

def write_sh(configpath):
    # Variables
    #configpath = os.path.abspath(argsconfig)

    # Read the configuration file
    config = json.load(open(configpath, 'r'))["environment_variables"]

    # write to "env.sh" file
    fh_env_sh = open(config["genelearn_path"]+"scripts/env.sh", 'w')

    fh_env_sh.write("export PATH=$PATH:" + config["genemark_path"] + '\n')
    fh_env_sh.write("export PATH=$PATH:" + config["genelearn_path"] + 'scripts/\n')
    fh_env_sh.write("export PATH=" + config["reftree_script_path"] + 'perl/bin:$PATH\n')
    fh_env_sh.write("export PERL5LIB=" + config["reftree_script_path"] + 'perl/lib:$PERL5LIB\n')
    fh_env_sh.write("export PYTHONPATH=$PYTHONPATH:" + config["python_site_packages_path"] + '\n')
    fh_env_sh.write("export PYTHONPATH=$PYTHONPATH:" + config["genelearn_path"] + 'scripts/\n')
    fh_env_sh.write("export REFTREE_DIR=" + config["reftree_db_path"] + '\n')
    fh_env_sh.write("export GENELEARN_DIR=" + config["genelearn_path"] + 'scripts/\n')
    fh_env_sh.write("export TMPDIR=" + config["TMPDIR"] + '\n')
    #fh_env_sh.write("module load tfmq\n")
    fh_env_sh.write("module load sqlite3\n")
    fh_env_sh.write("module load biopython\n")

    fh_env_sh.close()

    # write to "single_taxon_training.sh" file
    fh_single_taxon_training_sh = open(config["genelearn_path"]+"scripts/single_taxon_training.sh", 'w')

    fh_single_taxon_training_sh.write("#!/bin/bash \n")
    fh_single_taxon_training_sh.write(". " + config["genelearn_path"]+"scripts/env.sh\n")
    string_block = """
    module load biopython
    module load sqlite3

    TAXON=$1
    OUT=$2
    CONFIG=$3

    single_taxon_training.py --taxid "$TAXON" --outfile "$OUT" --config "$CONFIG"
    """
    fh_single_taxon_training_sh.write(string_block)

    fh_single_taxon_training_sh.close()

    os.chmod(config["genelearn_path"]+"scripts/single_taxon_training.sh",0750)
    os.chmod(config["genelearn_path"]+"scripts/env.sh",0750)



def main():

    parser = argparse.ArgumentParser(description='A script to run GeneLearn pipeline' )
    parser.add_argument('-c', '--config', help ="json formatted config file", required= True)
    parser.add_argument('-o', '--output', help ="Reftree database directory location",required= True) 
    args = parser.parse_args()


    # Variables
    configpath = os.path.abspath(args.config)
    write_sh(configpath)
    
    
    sts = os.path.join(os.environ['GENELEARN_DIR'],"single_taxon_training.sh")

    # Read the configuration file
    config = json.load(open(configpath, 'r'))["create_training_data"]

    # If an id list is specified in config write it otherwise query everything from the root node \
    if config["trainingid"]:
        node = str(config["trainingid"])
        print("Loading subtree taxon id from config file: %s" % node)
    else:
        node = str(1)
        print("No subtree specified, training on all sequences")
    
    # Open a subprocess to create a new reftree database directory, have it execute a shell \
    # script that calls single_taxon_training and writes the results to the reftree database directory

    if config["taskfarmer"] == "True":
        numworkers = str(config["numworkers"])
        resources = config["resources"]
        reftreeopts = ["reftree.pl", "--db", "genomic","--slots",numworkers, "--resources", resources,"--keep", "--foreach", node, args.output,\
    "--",sts, "__TAXON__","__OUTFILE__", configpath]
    else:
        resources = config["resources"]
        reftreeopts = ["reftree.pl", "--db", "genomic","--resources", resources, "--keep", "--foreach", node, args.output,\
    "--",sts, "__TAXON__","__OUTFILE__", configpath]
    #print(reftreeopts)
    p1 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    reftreeout, reftreeerr= p1.communicate()
    #print(reftreeout)
    
#prepare for training and testing datasets
    if config["testing_data"] == "True":
        formatter_opts = ["training_data_formatter.py", "--config", args.config, "--trainfile", config["training_vector_file_name"], "--testfile", config["testing_vector_file_name"], "--reftree", config["local_reftree_db_path"]]
        p2 = subprocess.Popen(formatter_opts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        formatter_out, formatter_err = p2.communicate()
    
        

if __name__ == '__main__':
    main()
