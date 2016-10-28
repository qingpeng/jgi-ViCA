read_PRC <- function(file_name)
{
    table <- read.table(file_name, sep = ",")
    random_100 = table[sample(nrow(table),1000),]
    random_100$V1 = as.numeric(gsub('\\(','',random_100$V1))
    random_100$V2 = as.numeric(gsub('\\)','',random_100$V2))
    return(random_100)
}

PRC_withPfam_run2 = read_PRC("/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam_run2_Pipeline/logistic/PRC/all.txt")
PRC_withPfam_run1 = read_PRC("/global/projectb/scratch/qpzhang/Full_Training/Pfam/PRC_all_logistic_log_scaling/all.txt")
PRC_withoutPfam = read_PRC("/global/projectb/scratch/qpzhang/Full_Training/Pfam/PRC_no_pfam_logistic/all.txt")
PRC_PfamOnly = read_PRC("/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam_run2_Pipeline/Pfam_only_logistic/PRC/all.txt")
PRC_PfamVfam = read_PRC("/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam_run/logistics_report/PRC/all.txt")

ggplot() + 
       geom_point(aes(V1, V2, colour = 'with Pfam,auPRC:0.9761'), data = PRC_withPfam_run2,size=0.8) + 
       geom_point(aes(V1, V2, colour = 'No Pfam'), data = PRC_withoutPfam,size=0.8) + 
       geom_point(aes(V1, V2, colour = 'Pfam Only, auPRC:0.9691'), data = PRC_PfamOnly,size=0.8) +
       geom_point(aes(V1, V2, colour = 'with Pfam + Vfam, auPRC:0.9775'), data = PRC_PfamVfam,size=0.8) +
       scale_colour_discrete("Features") +
        xlab("Recall") + ylab("Precision") 

