file_rank_obj = open("/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam_run/all.vect.rank", 'r')
file_out_obj = open("predictionAndLabels_family.txt.in_order", 'w')
file_prediction_obj = open("predictionAndLabels_family.txt", 'r')
file_family_index_obj = open("all.vect.family.index", 'r')
file_order_index_obj = open("all.vect.order.index", 'r')

order_id = {}# real tax-id -> 0,1,2,3
for line in file_family_index_obj:
    line = line.rstrip()
    fields = line.split()
    order_id[fields[1]] = fields[0]

family_id = {}# 0,1,2 -> real tax-id
for line in file_family_index_obj:
    line = line.rstrip()
    fields = line.split()
    family_id[fields[0]] = fields[1]

order_dict = {}


for line in file_rank_obj: # real tax-id family -> real order_id
    line = line.rstrip()
    tax_dict = {}
    taxonomy_list = line.split()[1].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:
        if "0" in tax_dict: # superkingdom exists
            order_dict[tax_dict['16']] = tax_dict['11']
file_order_family_obj = open("file_order_family.index", 'w')
for key in order_dict:
	file_order_family_obj.write(key+' '+order_dict[key]+'\n')


for line in file_prediction_obj:# all family id from 0...]: 

    line = line.rstrip()
    fields = line[1:-2].split(",")
    newline = "(" + order_id[order_dict[family_id[str(int(float(fields[0])))]]] + \
              "," + order_id[order_dict[family_id[str(int(float(fields[1])))]]] + ")"
    file_out_obj.write(newline+'n')


file_out_obj.close()


