import sys

file_train_obj = open(sys.argv[1], 'r')
file_train_out_obj = open(sys.argv[1]+'.num', 'w')
file_test_obj = open(sys.argv[2], 'r')
file_test_out_obj = open(sys.argv[2]+'.num', 'w')

file_index_obj = open(sys.argv[3], 'w') # correlate number 0....1,2, to real tax id number

tax_id_list = []

for line in file_train_obj:
    line = line.rstrip()
    fields = line.split()
    tax_id = fields[0]
    try:
        tax_id_num = tax_id_list.index(tax_id)
    except ValueError:
        tax_id_num = len(tax_id_list)
        tax_id_list.append(tax_id)
#    print tax_id_num
    new_line = str(tax_id_num) + ' ' + ' '.join(fields[1:])
    file_train_out_obj.write(new_line+'\n')


for line in file_test_obj:
    line = line.rstrip()
    fields = line.split()
    tax_id = fields[0]
    try:
        tax_id_num = tax_id_list.index(tax_id)
    except ValueError:
        tax_id_num = len(tax_id_list)
        tax_id_list.append(tax_id)
    print tax_id_num
    new_line = str(tax_id_num) + ' ' + ' '.join(fields[1:])
    file_test_out_obj.write(new_line+'\n')
    
for i in range(len(tax_id_list)):
    file_index_obj.write(str(i) + ' ' + tax_id_list[i] + '\n')

    
        