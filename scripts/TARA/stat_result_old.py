import sys
import random


vfam_obj = open("all_vfam.result", 'r')
img_obj = open("all_img.result", 'r')
hit = set()

evaluedict = {}

hit_vfam = set()
for line in vfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit.add(fields[2])
        hit_vfam.add(fields[2])
        evaluedict[fields[2]] = fields[5]

hit_img = set()

for line in img_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit.add(fields[2])
        hit_img.add(fields[2])
        evaluedict[fields[2]] = fields[5]
prediction_obj = open("all_4k_prediction.out", 'r')

tp = 0
fn = 0
tp_sum = 0
fn_sum = 0
for line in prediction_obj:
    line = line.rstrip()

    fields = line.split()
    if fields[0] in hit:
        if fields[3] == "1.0":
            tp += 1
            tp_sum += float(evaluedict[fields[0]])
        else:
            fn += 1
            fn_sum += float(evaluedict[fields[0]])
        if fields[0] in hit_vfam:
            vfam_mark = 1
        else:
            vfam_mark = 0
        if fields[0] in hit_img:
            img_mark = 1
        else:
            img_mark = 0
        print line+' '+evaluedict[fields[0]]+' '+str(vfam_mark)+' '+str(img_mark)

print tp, tp_sum, tp_sum/tp
print fn, fn_sum, fn_sum/fn



