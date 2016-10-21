import os
import subprocess
import argparse
import tempfile

parser = argparse.ArgumentParser(description='A script to create array jobs list' )
parser.add_argument('-t', '--subtree', help ="subtree to get", required= True)
parser.add_argument('-c', '--config', help ="json formatted config file", required= True)
parser.add_argument('-d', '--tempdir', help ="temp dir", required= True)
parser.add_argument('-o', '--output', help ="output file",required= True) 
args = parser.parse_args()

subtree = args.subtree
tempdir = args.tempdir
config_file = args.config
file_out_obj = open(args.output, 'w')


command_line = "reftree.pl --db genomic --subtree " + str(subtree) + "| grep \">\""

print command_line
output = subprocess.check_output(command_line, shell=True,universal_newlines=True)
#print output
lines = output.split('\n')
tax_set = set()
#print lines
for line in lines:
    line = line.rstrip()
#    print line
    if "taxid=" in line:
        taxid = line.split("=")[-1]
        print taxid
        tax_set.add(int(taxid))
    

#print tax_set
#for i in range(1000):
#    os.makedirs(tempdir+"/"+str(i))
    
for tax in tax_set:
#    print "here"
#    print tax
    dir_name = int(tax)%1000
    command_line =  os.environ.get('GENELEARN_DIR')+'/single_taxon_training.sh '+str(tax)+' '+tempdir+'/'+str(dir_name)+'/'+str(tax)+'.vect '+config_file
    file_out_obj.write(command_line+'\n')
    
file_out_obj.close()





