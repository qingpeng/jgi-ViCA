
## Model Tuning
split into training and testing
```angular2html
scripts/5_create_training_testing_with_seq_name.py
```
model training
```angular2html
$SPARKPATH/bin/spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```

model evaluation
```
$SPARKPATH/bin/spark-submit ./scripts/spark_evaluating_model.py testing.vect training.vect_model/ testing.vect.prediction testing.vect.report testing.vect.prc.png
```


## Training

model training (small data set)
```bash
$SPARKPATH/bin/spark-submit ./scripts/spark_training_model.py training.vect training.vect_model
```

## Get kmer/codon only features for full testing set
```bash
python ~/Github/jgi-ViCA/scripts/pick_vectors_by_feature.py -d ~/Github/jgi-Vi
CA/scripts/model/all_segment.fasta.vect.feature_index -f 0_1 -i all_segment.
fasta.vect.family.testing.svmlib.no4 -o all_segment.fasta.vect.family.testing.svmlib.no4.0_1
```