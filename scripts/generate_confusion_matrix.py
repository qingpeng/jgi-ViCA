import sys

fi =open(sys.argv[1],'r')

wordDic = [
"Archaea" ,"Bacteria" ,"ssRNAPhage" ,"ssRNAVirus" ,"dsRNAPhage","dsRNAVirus","dsDNAPhage","dsDNAVirus" ,"ssDNAPhage" ,"ssDNAVirus" ,"Retroviruses" ,"Eukaryota" ,"Mitochondrion","Chloroplast"]

labels = wordDic


print 'N/A'+' '+' '.join(labels)


list = []
i = 0
for line in fi:
    line = line.rstrip()
    if len(list) == 14:
        print_line = ' '.join(list)
        print labels[i]+' '+print_line
        list = [line]
        i +=1
    else:
        list.append(line)

print_line = ' '.join(list)
print labels[i]+' '+print_line
list = [line]
