import sys

f_in_obj = open(sys.argv[1],'r')
f_out_obj = open(sys.argv[2],'w')

num = 1
count= 0
for line in f_in_obj:
    if line[0] == '>':
    
        line = ">"+str(num)+"_"+line[1:]
        f_out_obj.write(line)
        count +=1
        if count == 6:
            num +=1
            count = 0
    else:
        f_out_obj.write(line)
