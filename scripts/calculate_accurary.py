# get precision from the multiple class prediction file
# this can be done in Spark run directly.
# cd ./predictionAndLabels_family/
# cat part-00* >all.txt
#

file_in_obj = open("all.txt", 'r')
file_out_obj = open("all.txt.out", 'w')

correct = {} 
incorrect = {} 
for line in file_in_obj:
	line = line.rstrip()
	fields = line[1:-1].split(',')
#	print fields
	if fields[0] == fields[1]:
		try:
			correct[fields[1]] += 1
		except KeyError:
			correct[fields[1]] = 1		
	else:
		try:
			incorrect[fields[1]] += 1
		except KeyError:
			incorrect[fields[1]] = 1

correct_total = 0
incorrect_total = 0

label_list = set(correct.keys() +  incorrect.keys())

for label in label_list:
	if label not in correct:
		correct[label] = 0
	if label not in incorrect:
		incorrect[label] = 0
	file_out_obj.write(label+' '+str(correct[label])+ ' ' + str(incorrect[label])+
		' ' + str(float(correct[label])/(correct[label] + incorrect[label])) +'\n')
	correct_total += correct[label]
	incorrect_total += incorrect[label]

file_out_obj.write(str(correct_total)+' '+str(incorrect_total)+' '+ 
	str(float(correct_total)/(correct_total+incorrect_total))+'\n')
	
