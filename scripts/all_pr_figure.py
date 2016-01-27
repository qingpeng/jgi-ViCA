import matplotlib

matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# 
# 
# file_in = open("all_metric.out",'r')
# 
# file_in.readline()
# 
# p={}
# r={}
# 
# for line in file_in:
#     line = line.rstrip()
#     if "logistic" in line:
#         method = "logistic"
#         c = 'ro'
#     elif "boosted" in line:
#         method = "boosted"
#         c = 'bo' 
#     elif "SVM" in line:
#         method = "SVM"
#         c = 'go'
#     
#     fields = line.split('\t')
#     if len(fields)<4:
#         continue
#     if method not in p:
#         p[method] = [float(fields[3])]
#         r[method] = [float(fields[4])]
#     else:
#         p[method].append(float(fields[3]))
#         r[method].append(float(fields[4]))
# print r
# print p
# plt.scatter(r['logistic'],p['logistic'],label = 'logistic',color='r')
# plt.scatter(r['boosted'],p['boosted'],label = 'boosted',color='g')
# plt.scatter(r['SVM'],p['SVM'],label = 'SVM',color='b')
# 
#     
# plt.legend(loc='lower left')
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.ylim([0.0, 1.05])
# plt.xlim([0.0, 1.0])
# plt.title("Recall-Precision by method")
# plt.savefig('method.pdf')

# 
# import matplotlib.pyplot as plt
# file_in = open("all_metric.out",'r')
# 
# file_in.readline()
# 
# p={}
# r={}
# 
# for line in file_in:
#     line = line.rstrip()
#     if "non-normalized" in line:
#         method = "non-normalized"
#         c = 'ro'
#     elif "normalized" in line:
#         method = "normalzied"
#         c = 'bo' 
# 
#     fields = line.split('\t')
#     if len(fields)<4:
#         continue
#     if method not in p:
#         p[method] = [float(fields[3])]
#         r[method] = [float(fields[4])]
#     else:
#         p[method].append(float(fields[3]))
#         r[method].append(float(fields[4]))
# print r
# print p
# plt.scatter(r['non-normalized'],p['non-normalized'],label = 'non-normalized',color='r')
# plt.scatter(r['normalzied'],p['normalzied'],label = 'normalzied',color='g')
# 
#     
# plt.legend(loc='lower left')
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.ylim([0.0, 1.05])
# plt.xlim([0.0, 1.0])
# plt.title("Recall-Precision by normalized or not")
# plt.savefig('normalized.pdf')


# 
# import matplotlib.pyplot as plt
# file_in = open("all_metric.out",'r')
# 
# file_in.readline()
# 
# p={}
# r={}
# 
# for line in file_in:
#     line = line.rstrip()
#     if "GeneMark" in line:
#         method = "GeneMark"
#         c = 'ro'
#     elif "Kmer" in line:
#         method = "Kmer"
#         c = 'bo' 
#     elif "projectb" in line:
#         method = "Both"
#         
#     fields = line.split('\t')
#     if len(fields)<4:
#         continue
#     if method not in p:
#         p[method] = [float(fields[3])]
#         r[method] = [float(fields[4])]
#     else:
#         p[method].append(float(fields[3]))
#         r[method].append(float(fields[4]))
# print r
# print p
# plt.scatter(r['GeneMark'],p['GeneMark'],label = 'GeneMark',color='r')
# plt.scatter(r['Kmer'],p['Kmer'],label = 'Kmer',color='g')
# plt.scatter(r['Both'],p['Both'],label = 'Both',color='b')
#     
# plt.legend(loc='lower left')
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.ylim([0.0, 1.05])
# plt.xlim([0.0, 1.0])
# plt.title("Recall-Precision by features")
# plt.savefig('features.pdf')



import matplotlib.pyplot as plt
file_in = open("all_metric.out",'r')

file_in.readline()

p={}
r={}

for line in file_in:
    line = line.rstrip()
    if "Testing_on_genome" in line:
        method = "genome"
        c = 'ro'

    elif "projectb" in line:
        method = "segment"
        
    fields = line.split('\t')
    if len(fields)<4:
        continue
    if method not in p:
        p[method] = [float(fields[3])]
        r[method] = [float(fields[4])]
    else:
        p[method].append(float(fields[3]))
        r[method].append(float(fields[4]))
print r
print p
plt.scatter(r['genome'],p['genome'],label = 'testing on genome',color='r')
plt.scatter(r['segment'],p['segment'],label = 'testing on segment',color='g')
#plt.scatter(r['Both'],p['Both'],label = 'Both',color='b')
    
plt.legend(loc='lower left')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title("Recall-Precision by testing")
plt.savefig('testing.pdf')


