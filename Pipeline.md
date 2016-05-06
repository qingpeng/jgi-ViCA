0. prepare for HMMER run


1. run create_training_data.py to generate raw vectors 
python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/create_training_data.py -c ~/Dropbox/Bitbucket/genelearn_paper/Script/config.json.shred.metagenemark_array_job_all_5k_pfam_combine -o pfam_run2 &

2. aggregate vectors into a whole file
python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/extract_vector_from_tmp_dir.py  -d ./ -o ../../Full_Training/Pfam/Pfam_run2/all.vect &

3. split whole vector file into training part and testing part
python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/training_data_formatter.py -c /global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/config.json.template -t training.vect -e testing.vect -s vector -r /global/projectb/scratch/ekirton/RefTree/ -v all.vect -g segment

4. 
  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/conert_split_pfam_vector.py -i testing.vect -o testing.vect.all -p testing.vect.pfam &
  

python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/process_pfam_vector.py -i testing.vect.all  -o testing.vect.all.sort &
  755  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/process_pfam_vector.py -i training.vect.all -o training.vect.all.sort &
  
  
    827  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/transform_log_pfam.py testing.vect.all.sort testing.vect.all.sort.log &
  828  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/transform_log_pfam.py training.vect.all.sort training.vect.all.sort.log &
  
  
  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/conert_split_pfam_vector.py -i testing.vect -o testing.vect.all -p testing.vect.pfam &
  
  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/transform_log_pfam.py testing.vect.pfam.sort testing.vect.pfam.sort.log
  
  
 Vfam Run!
====
 python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/create_training_data.py -c ~/Dropbox/Bitbucket/genelearn_paper/Script/config.json.shred.metagenemark_array_job_all_5k_pfam_combine -o pfam_combine_vfam

 python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/extract_vector_from_tmp_dir.py  -d ./ -o ../../Full_Training/Pfam/Vfam_run/all.vect &

  python ~/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/training_data_formatter.py -c /global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/config.json.template -t training.vect -e testing.vect -s vector -r /global/projectb/scratch/ekirton/RefTree/ -v all.vect -g segment &

  