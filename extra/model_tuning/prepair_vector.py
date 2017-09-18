
file_in_obj = open("/Users/qingpeng/all_segment.fasta.vect.family.training.svmlib.no4.0_1.1x_nonvirus", 'r')

virus_out = open("virus.vector", 'w')
nonvirus_out = open("nonvirus.vector", 'w')

for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    output = ''
    for field in fields[1:]:
        s = field.split(":")
        output = output+' '+s[1]
    if fields[0] == '0':
        nonvirus_out.write(output+'\n')
    else:
        virus_out.write(output+'\n')

