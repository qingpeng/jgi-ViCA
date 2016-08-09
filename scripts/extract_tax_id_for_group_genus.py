file_rank_obj = open("all.vect.rank", 'r')
file_out_obj = open("genus_group.list", 'w')

tax_bact = set()
tax_virus = set()
tax_euk = set()
tax_arch = set()

for line in file_rank_obj:
    line = line.rstrip()
    tax_dict = {}
    taxonomy_list = line.split()[1].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:
        if "0" in tax_dict:
            if tax_dict["0"] == "2157":
                tax_arch.add(tax_dict["20"])
            elif tax_dict["0"] == "2":
                tax_bact.add(tax_dict["20"])
            elif tax_dict["0"] == "2759":
                tax_euk.add(tax_dict["20"])
            elif tax_dict["0"] == "10239":
                tax_virus.add(tax_dict["20"])
                
for item in tax_bact:
    file_out_obj.write(item+' bact\n')
for item in tax_virus:
    file_out_obj.write(item+' virus\n')
for item in tax_euk:
    file_out_obj.write(item+' euk\n')
for item in tax_arch:
    file_out_obj.write(item+' arch\n')
    
file_out_obj.close()


