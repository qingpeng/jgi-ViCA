import sys
import random


vfam_obj = open("all_2k_contigs.fa.vfam", 'r')
pfam_obj = open("all_2k_contigs.fa.pfam", 'r')
hit_vfam = set()
hit_pfam = set()


for line in vfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit_vfam.add(fields[2])

for line in pfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit_pfam.add(fields[2])

prediction_obj = open("all_2k_prediction.out", 'r')

n_2 = 0
n_3 = 0
n_4 = 0
n_5 = 0

s_2 = 0
s_3 = 0
s_4 = 0
s_5 = 0

p_2 = 0
p_3 = 0
p_4 = 0
p_5 = 0

v_2 = 0
v_3 = 0
v_4 = 0
v_5 = 0


tp = 0
fp = 0
tn = 0
fn = 0

score = {}

all_score = 0

for line in prediction_obj:
    line = line.rstrip()

    fields = line.split()
    score[fields[0]] = float(fields[2])
    all_score += float(fields[2])
    if int(fields[1]) >= 5000:
        n_5 += 1
        s_5 += float(fields[2])
        if fields[0] in hit_pfam:
            p_5 += 1
        if fields[0] in hit_vfam:
            v_5 += 1
    elif int(fields[1]) >= 4000:
        n_4 += 1
        s_4 += float(fields[2])
        if fields[0] in hit_pfam:
            p_4 += 1
        if fields[0] in hit_vfam:
            v_4 += 1
    elif int(fields[1]) >= 3000:
        n_3 += 1
        s_3 += float(fields[2])
        if fields[0] in hit_pfam:
            p_3 += 1
        if fields[0] in hit_vfam:
            v_3 += 1
    else:
        n_2 += 1
        s_2 += float(fields[2])
        if fields[0] in hit_pfam:
            p_2 += 1
        if fields[0] in hit_vfam:
            v_2 += 1
    if float(fields[2]) >=0.9:
        if fields[0] in hit_vfam:
            tp += 1
        else:
            fp += 1
    else:
        if fields[0] in hit_vfam:
            fn += 1
        else:
            tn += 1


print 'tp', str(tp)
print 'fp', str(fp)
print 'fn', str(fn)
print 'tn', str(tn)

print 'pv5',p_5,v_5,s_5,n_5,s_5/n_5
print 'pv4',p_4,v_4,s_4,n_4,s_4/n_4
print 'pv3',p_3,v_3,s_3,n_3,s_3/n_3
print 'pv2',p_2,v_2,s_2,n_2,s_2/n_2



pfam_score_sum = 0

for i in hit_pfam:
    pfam_score_sum += score[i]

vfam_score_sum = 0

for i in hit_vfam:
    vfam_score_sum += score[i]

pfam_score_ave = pfam_score_sum/len(hit_pfam)
vfam_score_ave = vfam_score_sum/len(hit_vfam)


print 'pfam_score_ave', pfam_score_ave, len(hit_pfam)
print 'vfam_score_ave', vfam_score_ave, len(hit_vfam)
print 'ave_score', all_score/len(score), len(score)
