#!/usr/common/usg/languages/python/2.7.4/bin/python
# This runs a pipeline that 
#usage create_training_data.py <refseq file location> <root directory>
import sys
import simplejson as json
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='A script to generate SVM training data from RefSeq genomes')
parser.add_argument('-r','--root', help='a root input directory', required=True)
parser.add_argument('-p','--par', help='an alternate location for the metamark parameter files', default=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),"gm_parameters"))
parser.add_argument('-c,','--config',help= 'a json formated config file', default='config.json')
args = parser.parse_args()

#Load config settings
with open(args.config) as cfgfile:
	config = json.load(cfgfile)
	
print(args.par)
#define paths
root = os.path.abspath(args.root)
print(root)
genomes = os.path.join(root,"data","traininggenomes")
print(genomes)
mod = os.path.join(root, "data","trainingmod")


#function to create shell script			
def run_metamark_shellscripts(config,root,genomes,mod, par):
	"""qsub all taxonomic groups to genepool.nersc.gov UGE grid engine"""
	
	
	def check_dir_structure(dir,genomes,mod):
		"""check if input and output dir structure exists, if not raise exception (input) or create (output)"""
		if not os.path.exists(os.path.join(genomes,dir)):
			raise Exception("An expected folder %s is missing from the traininggenome directory" % dir)
		if not os.path.exists(os.path.join(mod,dir)):
			os.makedirs(os.path.join(mod,dir))
			
	def check_error_and_output_dir_structure(dir,root):
		"""check if batch output and batch error dir structure exists, if not create"""
		if not os.path.exists(os.path.join(root,"batch_output", dir)):
			os.makedirs(os.path.join(root,"batch_output",dir))
		if not os.path.exists(os.path.join(root,"batch_error", dir)):
			os.makedirs(os.path.join(root,"batch_error",dir))
	
				
	#return metamark parameters based on 
	def genome_type(dir):
		"""set metamark parameters"""
		if dir in["plant", "vertebrate_mammalian","vertebrate_other","protozoa","fungi", "invertebrate"]:
			type = "euk"
		elif dir in ["archaea", "bacteria","plastid","mitochondrion","plasmid"]:
			type = "prok"
		elif dir in ["ssRNAviruses","ssDNAviruses","dsRNAviruses","dsDNAviruses", "retroviruses", "viruses_other"]:
			type = "virus"
		else:
			raise Exception("Genome directory %s does not match with a known organism type, the correct parameters are not known" % dir)
		return type
		
	
	def calculate_threads(dir,genomes):
		"""Create the -t option string for qsub"""
		dirpath = os.path.join(genomes, dir)
		p1 = subprocess.Popen(["ls", dirpath], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
		p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
		maxt = p2.communicate()[0]
		p2.stdout.close()
		range = "1-" + str(maxt)
		return str(range.strip())
		
	def parameter_file(type, par):
		path = os.path.join(par,config["parfile"][type])
		return path
		
		
		
	# create and submit qsub jobs
	for dir in config["dir"]["data"]["traininggenomes"]:
		gmpath = os.path.join(os.path.dirname(__file__),"genemarksub.sh")
		check_dir_structure(dir, genomes, mod)
		check_error_and_output_dir_structure(dir,root)
		absinpath = os.path.join(genomes, dir)
		absrunpath = os.path.join(mod, dir)
		type = genome_type(dir)
		parfile = parameter_file(type,par)
		metamarkparams = config["metamarkparams"][type]
		qsubparams = config["qsubparams"][type] + " -o " + os.path.abspath(os.path.join(root,"batch_output",dir)) + " -e " + os.path.abspath(os.path.join(root,"batch_error",dir))
		threads = calculate_threads(dir,genomes)
		qsubcommand = "qsub " + qsubparams + " -t " + threads + " " + gmpath +" " + absinpath + " " + absrunpath + " \"" + metamarkparams + " --par " + parfile + "\""
		print(qsubcommand)	
		output = subprocess.check_output(qsubcommand,shell=True)
		print(output)
run_metamark_shellscripts(config=config,root=args.root,genomes=genomes,mod=mod, par=args.par)


#task 2
# for each file in all directories build a matrix of emission probabilities and build a class vector


