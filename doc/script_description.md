

### Other helper scripts
select vectors from specific features:
```angular2html
scripts/pick_vectors_by_feature.py
```


## Description of scripts
- scripts/sklearn_training_model.py : training model using sklearn, with weighted option
- scripts/sklearn_evaluate_model.py: evaluate model using sklearn
- scripts/2_feature_genemark.py: extract genemark feature only
- scripts/2_feature_genemark_pfam_vfam.py: extract genemark/ codon usage and 
pfam/vfam/IMG virus hits
- scripts/2_feature_genemark_with_extra_features.py: experimental, genemark 
feature and other new features, like codon proportion, etc.

- scripts/pick_kmer_codon_from_svmlib.py: pick kmer and codon features from svmlib with all features

python ../../scripts/pick_kmer_codon_from_svmlib.py  testing.vect.200 testing.vect.200_0_1

