import sys

file_in = open(sys.argv[1],'r')
file_out = open(sys.argv[2],'w')
#to_pick = int(sys.argv[3])

k=0

for line in file_in:
    if k==10:
        file_out.write(line)
        k=0
    else:
        k+=1
        

