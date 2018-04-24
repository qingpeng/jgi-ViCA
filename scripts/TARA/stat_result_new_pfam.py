import sys
import random


vfam_obj = open("all_2k_contigs.fa.vfam", 'r')
pfam_obj = open("all_2k_contigs.fa.protein.pfam", 'r')
hit_vfam = set()
hit_pfam = set()

diamond_obj = open("all_2k_contigs.fa.diamond.virus", 'r')

hit_diamond = set()

for line in diamond_obj:
    line = line.rstrip()
    fields = line.split()
    hit_diamond.add(fields[0])
    print fields


for line in vfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit_vfam.add(fields[2])

for line in pfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        # print line
        fields = line.split()
        #print fields[4], fields[-1]
        # print fields
        if float(fields[4]) <= 0.00001:

            if ('RNA dependent RNA polymerase' in line or
                'Minor capsid protein' in line or
                    'Major capsid protein' in line):

                hit_pfam.add(fields[2][:-2])

#print hit_pfam
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

tp_diamond = 0
fp_diamond = 0
tn_diamond = 0
fn_diamond = 0

tp_all = 0
fp_all = 0
tn_all = 0
fn_all = 0


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
    if float(fields[2]) >0.98:

        if fields[0] in hit_vfam or fields[0] in hit_pfam:
            tp += 1
        else:
            fp += 1
        if fields[0] in hit_diamond:
            tp_diamond += 1
        else:
            fp_diamond += 1

        if fields[0] in hit_vfam or fields[0] in hit_pfam or fields[0] in hit_diamond:
            tp_all += 1
        else:
            fp_all += 1


    else:
        if fields[0] in hit_vfam or fields[0] in hit_pfam:
            fn += 1
        else:
            tn += 1
        if fields[0] in hit_diamond:
            fn_diamond += 1
        else:
            tn_diamond += 1

        if fields[0] in hit_vfam or fields[0] in hit_pfam or fields[0] in hit_diamond:
            fn_all += 1
        else:
            tn_all += 1

print 'tp', str(tp)
print 'fp', str(fp)
print 'fn', str(fn)
print 'tn', str(tn)

print "diamond:\n"
print 'tp', str(tp_diamond)
print 'fp', str(fp_diamond)
print 'fn', str(fn_diamond)
print 'tn', str(tn_diamond)

print "all:\n"
print 'tp', str(tp_all)
print 'fp', str(fp_all)
print 'fn', str(fn_all)
print 'tn', str(tn_all)


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
