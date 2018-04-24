import sys

file_old_result_obj = open("all_2k_prediction.out.label", 'r')

euk_non_euk_obj = open("all_2k.libsvm.prediction_euk_non_euk", 'r')
non_euk_virus_obj = open(
    "all_2k.libsvm.prediction_non_euk_non_virus_and_virus", 'r')

full_training_obj = open("all_2k.libsvm.prediction_virus_non_virus", 'r')

output_obj = open("all_2k_prediction.out.label.prediction_combined", 'w')


score_euk_non_euk = []
pred_euk_non_euk = []
for line in euk_non_euk_obj:
    line = line.rstrip()
    fields = line.split()
    score_euk_non_euk.append(fields[0])
    pred_euk_non_euk.append(fields[1])

score_non_euk_virus = []
pred_non_euk_virus = []
for line in non_euk_virus_obj:
    line = line.rstrip()
    fields = line.split()
    score_non_euk_virus.append(fields[0])
    pred_non_euk_virus.append(fields[1])

score_full_training = []
pred_full_training = []
for line in full_training_obj:
    line = line.rstrip()
    fields = line.split()
    score_full_training.append(fields[0])
    pred_full_training.append(fields[1])

count = 0
for line in file_old_result_obj:
    line = line.rstrip()
    print_line = line+' '+str(score_euk_non_euk[count])+' '+pred_euk_non_euk[
                 count]+' '+score_non_euk_virus[count]+' '+pred_non_euk_virus[
                 count]+' '+score_full_training[count]+' '+pred_full_training[
                 count]+'\n'

    output_obj.write(print_line)
    count += 1

