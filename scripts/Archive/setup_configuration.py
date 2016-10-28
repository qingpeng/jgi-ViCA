#!/usr/bin/env python
# write two scripts: env.sh and single_taxon_training.sh
#


import argparse
import simplejson as json
import os


parser = argparse.ArgumentParser(description='A script to create relative scripts to reflect current environment variables' )
parser.add_argument('-c', '--config', help ="json formatted config file", required= True)
args = parser.parse_args()


# Variables
configpath = os.path.abspath(args.config)

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

