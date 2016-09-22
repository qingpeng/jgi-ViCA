#!/usr/bin/env python
import argparse

def main():

    parser = argparse.ArgumentParser(description='A script to generate k-mer coposition frequency')
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
        vect = fields[3].split(' ')
        
        for file_obj in file_handles[1:]:
            line = file_obj.readline()
            line = line.rstrip()
            fields = line.split('\t')
            if len(fields) == 4:
                vect.extend(fields[3].split(' '))
        
        if length >= cutoff:
            print_line = seq_id+'\t'+length+'\t'+seq_des+'\t'+ ' '.join(vect)
            file_output_obj.write(print_line+'\n')
            
            
    
    file_output_obj.close()
    
    
if __name__ == '__main__':
    main()