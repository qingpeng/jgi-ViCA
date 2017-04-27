
prepare vectors file with k-mer/codon usage features only
```angular2html
python ~/Github/jgi-ViCA/scripts/pick_vectors_by_feature.py -d ../all_segment.fasta.vect.feature_index -f 0_1 -i all_segment.fasta.vect.family.training.svmlib.no4 -o all_segment.fasta.vect.family.training.svmlib.no4.0_1 &
```
run spark on vectors without pfam
```angular2html
sbatch training_no_pfam.sl
```

