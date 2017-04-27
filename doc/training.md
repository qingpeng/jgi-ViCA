
###Pipeline to do training using simulated segments from RefSeq genomes


#### number of segments
qpzhang@genepool14:/global/projectb/scratch/qpzhang/RefTree.tmp/Output_Shred_Full$ grep -c ">" all_segment.fa
8708600

8708600/200
43543

qpzhang@genepool14:/global/projectb/scratch/qpzhang/RefTree.tmp/Output_Shred_Full$ wc taxa.txt
   43568   174272 10635013 taxa.txt
   
   
qpzhang@genepool14:/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow$ mv ../../RefTree.tmp/Output_Shred_Full/all_segment.fa ./


#### use nextflow to extract features
```angular2html
nextflow run ~/Bitbucket/jgi-genelearn/scripts/genelearn.nf
```
[2]+  Running                 nextflow run ~/Bitbucket/jgi-genelearn/scripts/get_vector_lite.nf -qs 2000 &