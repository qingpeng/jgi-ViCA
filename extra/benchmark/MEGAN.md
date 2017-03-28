

```python
from IPython.display import Image
```

# blast megan to evaluate performance

- download NR from NCBI

ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

Jan 19, 2017


- diamond against NR

diamond v0.8.34.96

does not use module version, does not work


- run megan

MEGAN v 6.6.6

abin file

http://ab.inf.uni-tuebingen.de/data/software/megan6/download/prot_acc2tax-Nov2016.abin.zip


## Try MGRAST to see how it works

### get the ~5k virus and nonvirus for testing
```
Genepool:
/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow

~/Bitbucket/jgi-genelearn/scripts/random_seq.py
from virus_segment.fa to virus_segment_5k.fa
from nonvirus_segment.fa to nonvirus_segment_5k.fa


Local:
/Users/qingpeng/GeneLearn/

Get rid of bad sequences for MGRAST 


python modify_fasta_for_mgrast.py nonvirus_segment_5k.fa >nonvirus_segment_5k_for_MGRAST.fa
python modify_fasta_for_mgrast.py virus_segment_5k.fa >virus_segment_5k_for_MGRAST.fa &

./seqtk seq -N ../virus_segment_5k_for_MGRAST.fa  >../virus_segment_5k_for_MGRAST_N.fa
./seqtk/seqtk seq -N nonvirus_segment_5k_for_MGRAST.fa >nonvirus_segment_5k_for_MGRAST_N.fa

grep -c ">" virus_segment_5k_for_MGRAST.fa
5120
grep -c ">" virus_segment_5k_for_MGRAST_N.fa # no description
4916
grep -c ">" virus_segment_5k.fa  # with description
5120

get 
ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz


```






# Try Diamond/MEGAN to see the performance

## use the segments used for MGRAST

diamond to get blast result to check if it is good...

use downloaded version of diamond.....


Path:
/global/projectb/scratch/qpzhang/NCBI_NR


Megan...

to offer seperate file for accession number -> taxonomy

it is protein sequences...

Nextflow has issue to run DIAMOND BLASTX 

Run directly.





## Use all alignment result agianst nr





### For ~5000 virus segments



Most of 4888 virus segments are assigned. (393 unassigned)

most of them are assigned as virus (161 assigned to non-virus)


### For ~5000 non-virus segments


Most of the 4929 nonvirus segments are assigned. (20 unassigned)

nearly all of them are assigned as non-virus 


#### This is not surprising. Since the segments can almost always aligned to some protein sequences in nr, since they are simulated from genomes.


## Remove the hits to the reference sequences where the segment is from


convert from .daa to .m8 

15489306

qpzhang@genepool14:/global/projectb/scratch/qpzhang/NCBI_NR$ wc prot.accession2taxid
  359376600  1437506400 13756184293 prot.accession2taxid
q


```
diamond view -a nonvirus_segment_5k_for_MGRAST_N.fa.daa -o nonvirus_segment_5k_for_MGRAST_N.fa.daa.m8 -f 6
diamond view -a virus_segment_5k_for_MGRAST_N.fa.daa -o virus_segment_5k_for_MGRAST_N.fa.daa.m8 -f 6

remove hits with the species , that is in the same family.....

filter_alignment.py

python filter_alignment.py virus_segment_5k.fa prot.accession2taxid virus_segment_5k_for_MGRAST_N.fa.top.m8 virus_segment_5k_for_MGRAST_N.fa.sensitive.top.filter.m8 &

python filter_alignment.py nonvirus_segment_5k.fa prot.accession2taxid nonvirus_segment_5k_for_MGRAST_N.fa.top.m8 nonvirus_segment_5k_for_MGRAST_N.fa.top.filter.m8 &


python filter_alignment.py virus_segment_5k.fa prot.accession2taxid virus_segment_5k_for_MGRAST_N.fa.all.m8 virus_segment_5k_for_MGRAST_N.fa.all.filter.m8 &
```
### Experiment1:

diamond output default... (25 top hits)

problem, 25 top hits all to the family where the query is from, no hits after filtering

### Experiment2:

diamond output top hits until >20\% less than top...

similar to above... bad hits to reference from other families...



### Experiment3:

diamond output all hits... with sensitive option, for longer query sequences

filter out with 100 hits to different family reference




get 100 effective hits after hits to the family where the query is from


```
python filter_alignment_100output.py virus_segment_5k.fa prot.accession2taxid virus_segment_5k_for_MGRAST_N.fa.sensitive.all.m8 virus_segment_5k_for_MGRAST_N.fa.sensitive.all.filter100.m8 &



```

## Annotation by MEGAN only using alignements against reference sequences from different family






### viurus segments

Out of 4916 virus segments, 4772 segments still have hits to species in different family from where the segment is from. 

Out of the 4772 segments with hits:

2213 are annotated as virus.(46.4%)

1211 are annotated as bacteria (927) archaea (20) Eukaryota (223)

573 are annotaed as "artificial sequences"

223 are annotated as "Not assigned"






```python
Image(filename='Images/virus_segment_5k_for_MGRAST_N.fa.sensitive.all.filter100.png')
## double click the figure to enlarge
```




![png](output_11_0.png)



### non- viurus segments

Out of 4949 virus segments, 4913 segments still have hits to species in different family from where the segment is from. 

Out of the 4913 segments with hits:

48 are annotated as virus.

4738 are annotated as "cellular organisms", bacteria (3196) archaea (52) Eukaryota (1220)

15 are annotaed as "artificial sequences"

22 are annotated as "Not assigned"


```python
Image(filename='Images/nonvirus_segment_5k_for_MGRAST_N.fa.sensitive.all.filter100.png')
## double click the figure to enlarge
```




![png](output_13_0.png)



## Other tax_level... - order -family  -  genus


also the same for Genelearn....






## Two experiments

genus... fewer hits.... fewer segments with hits..... many segments don't have genus information, can't decide if they are from the same genus or not......

family... if the same genus... even without family infor,,, still ...
do we really need this????

make sure definitely does not from the same taxonomic levels....

only look at the segments with hits..... reigious... criteria....

We only keep the alignment between query and target who do not belong to the same taxonomic rank explictly and definitely. Any alignment with query or target without the taxonomic information on that rank level will be discarded. 

If target accesion number does not have tax_id information, discard such alignment.


2]+  Done                    python ~/Bitbucket/jgi-genelearn/scripts/filter_close_alignment_for_MEGAN.py -s nonvirus_segment_5k.fa -a prot.accession2taxid -i nonvirus_segment_5k_for_MGRAST_N.fa.sensitive.all.m8 -t 100 -o nonvirus_segment_5k_for_MGRAST_N.fa.sensitive.all.m8.filter_keep_ambiguous2 --filter_option 0



## Run the segments using GeneLearn, with the application

### as the testing datasets!!

Use the same sequences in the testing data sets,,, with vectors....

get the original fasta files from the all fasta file!!!!!

directly comparison!!!!!




qpzhang@genepool14:/global/projectb/scratch/qpzhang/NCBI_NR/Nextflow$ nextflow run ~/Bitbucket/jgi-genelearn/scripts/feature_extraction.nf --segment_file /global/projectb/scratch/qpzhang/NCBI_NR/Nextflow/virus_segment_5k.fa --out /global/projectb/scratch/qpzhang/NCBI_NR/Nextflow/virus_segment_5k.fa.vect --chunkSize 50 --package_path /global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/ --genemark_path /global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite/ --hmmer_path /global/homes/q/qpzhang/bin/hmmer-3.1b2-linux-intel-x86_64/ --hmmer_db /global/homes/q/qpzhang/Pfam_DB/

nextflow run ~/Bitbucket/jgi-genelearn/scripts/feature_extraction.nf --segment_file /global/projectb/scratch/qpzhang/NCBI_NR/Nextflow/nonvirus_segment_5k.fa --out /global/projectb/scratch/qpzhang/NCBI_NR/Nextflow/nonvirus_segment_5k.fa.vect --chunkSize 50 --package_path /global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/ --genemark_path /global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite/ --hmmer_path /global/homes/q/qpzhang/bin/hmmer-3.1b2-linux-intel-x86_64/ --hmmer_db /global/homes/q/qpzhang/Pfam_DB/


5000 segments... 100 processes, 8 minutes to finish for each...


python ~/GDrive/Development/Bitbucket/jgi-genelearn/scripts/vect2svmlib.py -v virus_segment_5k.fa.vect -f all_segment.fasta.vect.feature_index -o virus_segment_5k.fa.vect.svmlib -t 1

python ~/GDrive/Development/Bitbucket/jgi-genelearn/scripts/vect2svmlib.py -v nonvirus_segment_5k.fa.vect -f all_segment.fasta.vect.feature_index.new_list -o nonvirus_segment_5k.fa.vect.svmlib -t 0

Qingpen-lm:GeneLearn qingpeng$ cat virus_segment_5k.fa.vect.svmlib nonvirus_segment_5k.fa.vect.svmlib >virus_nonvirus_segment_5k.fa.vect.svmlilb

python ~/GDrive/Development/Bitbucket/jgi-genelearn/scripts/Flask/modify_label.py virus_nonvirus_segment_5k.fa.vect.svmlilb 19421 >virus_nonvirus_segment_5k.fa.vect.svmlilb2

[[4890,  228],

[1251, 3869]])
       
       

## compare speed and performance


```python

```
