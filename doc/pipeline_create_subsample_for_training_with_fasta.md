### The pipeline to pick randomly subset for training and testing, as well as the corresponding segment sequences


1. split vector file into training and testing by different taxonomic level (order, family, genus)

```
python ~/Github/jgi-ViCA/scripts/5_create_training_testing_with_seq_name.py -v all_segment.fasta.vect -r True &
```

2. split vector files into virus-only and nonvirus-only subsets

```
python ~/Github/jgi-ViCA/scripts/split_vector_into_virus_nonvirus.py all_segment.fasta.vect.family.testing.svmlib all_segment.fasta.vect.family.testing.svmlib.virus all_segment.fasta.vect.family.testing.svmlib.nonvirus &
python ~/Github/jgi-ViCA/scripts/split_vector_into_virus_nonvirus.py all_segment.fasta.vect.family.training.svmlib all_segment.fasta.vect.family.training.svmlib.virus all_segment.fasta.vect.family.training.svmlib.nonvirus &
```

3. get the number of records in each subsets.

```
wc all_segment.fasta.vect.family.testing.svmlib.nonvirus all_segment.fasta.vect.family.testing.svmlib.virus >all_segment.fasta.vect.family.testing.svmlib.wc
wc all_segment.fasta.vect.family.training.svmlib.nonvirus all_segment.fasta.vect.family.training.svmlib.virus >all_segment.fasta.vect.family.training.svmlib.wc

qpzhang@genepool14:/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/Split_with_name$ more *.wc
::::::::::::::
all_segment.fasta.vect.family.testing.svmlib.wc
::::::::::::::
   1579000  563265401 8885107856 all_segment.fasta.vect.family.testing.svmlib.nonvirus
     58000   21153543  334347475 all_segment.fasta.vect.family.testing.svmlib.virus
   1637000  584418944 9219455331 total
::::::::::::::
all_segment.fasta.vect.family.training.svmlib.wc
::::::::::::::
    5076200  1878265615 29471358010 all_segment.fasta.vect.family.training.svmlib.nonvirus
     150600    49967056   800591102 all_segment.fasta.vect.family.training.svmlib.virus
    5226800  1928232671 30271949112 total
```    
    
4. randomly select 5k for testing and 25k segments for training from each subset

```
python pick_random_line.py all_segment.fasta.vect.family.testing.svmlib.virus 58000 5000 virus_testing_5k.vect &
python pick_random_line.py all_segment.fasta.vect.family.testing.svmlib.nonvirus 1579000 5000 non_virus_testing_5k.vect &
python pick_random_line.py all_segment.fasta.vect.family.training.svmlib.virus 150600 25000 virus_training_25k.vect &
python pick_random_line.py all_segment.fasta.vect.family.training.svmlib.nonvirus 5076200 25000 non_virus_training_25k.vect &
```


5. combine vector files for training and testing, with segment name as the first column

```
cat non_virus_testing_5k.vect virus_testing_5k.vect >testing_5k.vect_with_segment_name &
cat non_virus_training_25k.vect virus_training_25k.vect >training_25k.vect_with_segment_name &
```

6. remove the first column of segment name from vector file, now the vector files can be used for training directly
```
python remove_segment_name_from_vector.py testing_5k.vect_with_segment_name testing_5k.vect
python remove_segment_name_from_vector.py training_25k.vect_with_segment_name training_25k.vect
```

7. get the corresponding fasta sequences from all_segment.fa for each segment for the randomly selected vectors file

```
python get_fasta_for_vector.py
```

output:
```
non_virus_testing_5k.vect.fasta  non_virus_training_25k.vect.fasta
virus_testing_5k.vect.fasta  virus_training_25k.vect.fasta
```


8. combine segment fasta files for training and testing together
```
cat non_virus_testing_5k.vect.fasta virus_testing_5k.vect.fasta >testing_5k.vect.fasta &
cat non_virus_training_25k.vect.fasta virus_training_25k.vect.fasta >training_25k.vect.fasta &
```

9. get 3 virus and non-virus segments sequences as testing file for pipeline

```angular2html
head -255  virus_testing_5k.vect.fasta >virus_testing_5k.vect.fasta.head255
head -255 non_virus_testing_5k.vect.fasta >non_virus_testing_5k.vect.fasta.head255
cat non_virus_testing_5k.vect.fasta.head255 virus_testing_5k.vect.fasta.head255 >virus_nonvirus_3seqs.fa
```