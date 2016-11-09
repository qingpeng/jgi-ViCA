
import sys
file_out = open(sys.args[1],'w')
length = sys.args[2]
taxid = sys.args[3]

block1 = """
{
  "environment_variables" : {
    "genemark_path" : "/global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite/",
    "genelearn_path" : "/global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/",
    "reftree_script_path" : "/global/homes/q/qpzhang/Dropbox/Bitbucket/jgi_reftree/",
    "python_site_packages_path" : "/global/homes/q/qpzhang/anaconda/lib/python2.7/site-packages/",
    "reftree_db_path" : "/global/projectb/scratch/ekirton/RefTree/",
    "TMPDIR" : "/global/projectb/scratch/qpzhang"
  }, 

  "create_training_data" : {
    "taskfarmer" : "False",
    "numworkers" : "50000",
    "resources" : "-l ram.c=5G,h_rt=0:30:00  -pe pe_slots 2",
    "method" : "metagenemark_v2",
    "shred" : "fixed",
    "lognorm" : {
      "shape" : "1.3335",
      "loc" : "3000",
      "scale": "1140"
    },"""
    
file_out.write(block1)

line = "    \"fixed\" : \""+str(length)+"\",\n"

file_out.write(line)    
#    "fixed" : "5000",
    
    
block2 = """    "shredsamples" : "200",
    "tmpdir": "/scratch/qpzhang/","""

file_out.write(block2)

line = "    \"trainingid\" : \""+str(taxid)+"\",\n"
file_out.write(line)

block3 = """    "mmp" : "/global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/gm_parameters/par_11.modified"
  },
  "training_data_formatter" : {
    "method" : "taxids",
    "testprop" : "0.2",
    "categories" : {
      "Archaea" : {"taxcat" : "2157"},
      "Bacteria" : {"taxcat" : "2"},
      "ssRNAPhage" : {"taxcat" : "439488", "gencode" : "11"},
      "ssRNAVirus" : {"taxcat" : "439488", "gencode" : "1"},
      "dsRNAPhage" : {"taxcat" : "35325", "gencode" : "11"},
      "dsRNAVirus" : {"taxcat" : "35325", "gencode" : "1"},
      "dsDNAPhage" : {"taxcat" : "35237", "gencode" : "11"},
      "dsDNAVirus" : {"taxcat" : "35237", "gencode" : "1"},
      "ssDNAPhage" : {"taxcat" : "29258", "gencode" : "11"},
      "ssDNAVirus" : {"taxcat" : "29258", "gencode" : "1"},
      "Retroviruses" : {"taxcat" : "35268"},			
      "Eukaryota" : {"taxcat" : "2759", "organelle" : "None"},
      "Mitochondrion" : {"taxcat" : "2759", "organelle": "mitochondrion"},
      "Chloroplast" : {"taxcat" : "2759", "organelle": "plastid:chloroplast"}
    },
    "levelattr" : {
      "level" : "Phylum"
    }
  }
}"""

file_out.write(block3)

