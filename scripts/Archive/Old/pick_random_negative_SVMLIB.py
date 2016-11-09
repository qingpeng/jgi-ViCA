import sys

fo = open(sys.argv[2],'w')
k=0
for line in open(sys.argv[1],'r'):
    line = line.rstrip()
    if k%100==0:
        if line[0] == '0':
            fo.write(line+"\n")
        k=0
    k = k+1

fo.close()
