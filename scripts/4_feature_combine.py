#!/usr/bin/env python
import argparse

def main():

    parser = argparse.ArgumentParser(
        description='A script to combine features into one vector file')
    parser.add_argument('feature_files', help="list of feature files to combine", nargs='+')
    parser.add_argument('--output', help= "Output file, space delimited format", default='all.vect')
    parser.add_argument('--length', help= "length cutoff", default=2000)

    args = parser.parse_args()

    file_lists = args.feature_files
    cutoff = args.length
    file_output_obj = open(args.output, 'w')
    
    file_handles = []
    
    for file in file_lists:
        file_obj = open(file,'r')
        file_handles.append(file_obj)
        
    file_1st = file_handles[0]
# Output format:
# seq_id'\t'seq_length'\t'seq_des'\t'vectors, separated by " "
#



    for line in file_1st:
        line = line.rstrip()
        fields = line.split('\t')
        seq_id = fields[0]
        seq_des = fields[2]
        length = fields[1]
        file_id = 0
        vect = fields[3].split(' ')
        # integrate feature id in the vector output file
        vect_id = [ str(file_id)+"_"+vec for vec in vect]
        
        for file_obj in file_handles[1:]:
            line = file_obj.readline()
            line = line.rstrip()
            fields = line.split('\t')
            file_id = file_id + 1
            if len(fields) == 4:
                vect = fields[3].split(' ')
            
                vect_id_this = [str(file_id)+"_"+vec for vec in vect]
            
                vect_id.extend(vect_id_this)
        
        if length >= cutoff:
            print_line = seq_id+'\t'+length+'\t'+seq_des+'\t'+ ' '.join(vect_id)
            file_output_obj.write(print_line+'\n')
            
            
    
    file_output_obj.close()
    
    
if __name__ == '__main__':
    main()