library("Biostrings")
library(stringr)
library(dplyr)

### добавление Uniprot ID к найденным белкам
# Pisum_sativum_v1a_prot - фаста файл белковых последовательностей 
pea_db <- readBStringSet("C:/Users/admin/Pisum_sativum_v1a_prot.fasta")
pea_names = names(pea_db)
sequence = paste(pea_db)
pea_db <- data.frame(pea_names, sequence)
pea_db$sequence <- as.character(pea_db$sequence)
pea_db$pea_names <- tolower(pea_db$pea_names)

#Pisum_sativum_genome_names.csv - названия белков по номенклатуре Uniprot
pea_prot_names <- read.csv("C:/Users/admin/Pisum_sativum_genome_names.csv", header = T)
pea_prot_names <- pea_prot_names[, 2:3]
colnames(pea_prot_names)[1] <- "pea_names" 
colnames(pea_prot_names)[1] <- "first_accession"

pea_db_test <- inner_join(pea_db, pea_prot_names, by = "pea_names")
pea_db_test <- inner_join(pea_db, pea_prot_names, by = "first_accession")
pea_db_test <- pea_db_test[, c(3,2,1)]
writeFasta(data = pea_db_test, filename = "Pea_proteome.fasta")

#######
#Pea_coat_comp.csv - таблица с данными по диффэкспресии
pea_table <- read.csv("C:/Users/admin/Pea_coat_comp.csv", header = T, sep = ";", dec = ",")
colnames(pea_table)[12:21] <- c("wild_076", "wild_077", "wild_078", "wild_079", "wild_080", "cult_081", "cult_082", "cult_083", "cult_084", "cult_085")
pea_table$first_accession <- gsub(pea_table$Accession, pattern = ";.*", replacement = "")


pea_db_test$full_name <- gsub(pea_db_test$full_name, pattern = "psat.*uniprot: ", replacement = "")
pea_table$Description <- gsub(pea_table$Description, pattern = "uniprot: ", replacement = "")
colnames(pea_db_test)[3] <- "Description"

pea_seq <- inner_join(pea_db_test, pea_table, by = c("first_accession", "Description"))
pea_seq <- pea_seq[, c(3, 2)]
writeFasta(data = pea_seq, filename = "Pea_seq.fasta")
##################################
