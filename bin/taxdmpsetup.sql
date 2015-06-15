create table gi_taxid_prot(
        gi int primary key,
        taxid int);
        .separator '    '
.import trainingdata/gi_taxid_prot.dmp gi_taxid_prot

