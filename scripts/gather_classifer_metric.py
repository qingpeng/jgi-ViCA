#!/usr/bin/env python

import os
import glob



def parse_vect(dir):
    files = glob.glob(dir+"/*.results")
    metric = {}
    for file in files:
        file_in_obj = open(file,'r')
        metric[os.path.basename(file)] = {}
        
        for line in file_in_obj:
            fields = line.rstrip().split(":")
            metric[os.path.basename(file)][fields[0]] = fields[1]
        file_in_obj.close()
    return metric
    
file_out_o = open("all_metric.out",'w')



metrics = ['accuracy','precision','recall','f1_score','auc','log_loss']
header = '\t'.join(['file']+metrics)
file_out_o.write(header+'\n')


folders = \
["/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_genome/non-normalized_model", 
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_genome/normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_genome_multiple_labels/non-normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_genome_multiple_labels/normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_segment/non-normalized_model", 
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_segment/normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_segment_multiple_labels/non-normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Genome/Testing_on_segment_multiple_labels/normalized_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/Both_Feature/Full_data/Multiple_label/dato_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/Both_Feature/Full_data/Single_label/dato_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/Both_Feature/Pro_data/Multiple_label/dato_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/Both_Feature/Pro_data/Single_label/dato_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/GeneMark/dato_model",
"/global/projectb/scratch/qpzhang/Genelearn_Project/Segment/Kmer/dato_model"]


for folder in folders:
    file_out_o.write(folder+'\n')
    
    results = parse_vect(folder)
    for key in results.keys():
        file_out_o.write(key+'\t')
        str = ''
        for metric in metrics:
            if metric in results[key]:
                str = str+ '\t' + results[key][metric]
            else:
                str = str+'\tN/A'
        file_out_o.write(str+'\n')
file_out_o.close()

        
    
    