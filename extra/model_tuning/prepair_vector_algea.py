
file_in_obj = open("/Users/qingpeng/Dropbox/Development/Github/jgi-ViCA/extra/model_tuning/Micromonas/tmpL1fNBb/shred.fa.vect.libsvm", 'r')

virus_out = open("algea.vector", 'w')


for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    output = ''
    for field in fields[1:]:
        s = field.split(":")
        output = output+' '+s[1]

    virus_out.write(output+'\n')

