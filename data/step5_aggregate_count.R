library(data.table)
library(dplyr)
library(stringr)

# Define the set of NCs that had non-zero values from Serratus
nc_have <- list.files("nc_pulls/",pattern = "^NC_")
nc_have_simple <- sapply(strsplit(nc_have, "[.]"), function(x) x[1])

# Redo logic to see what's missting
tbl <- fread("Table_human_viruses.txt", header = FALSE)
tbl$reactivation_candidate <- grepl("Herpesviridae|Polyomaviridae|Parvoviridae|Adenoviridae",tbl$V2)

tbl$serratus_hits <- sapply(1:dim(tbl)[1], function(x){
  ncs <- str_trim(strsplit(tbl[["V6"]][x], ",")[[1]])
  sapply(ncs, function(one_nc){
    ifelse(any(grepl(one_nc, nc_have)),
           (fread(paste0("nc_pulls/",nc_have[grep(one_nc, nc_have)])) %>% dim())[1],
           0)
  }) %>% sum()
})

tbl$serratus_high_confidence_hits <- sapply(1:dim(tbl)[1], function(x){
  ncs <- str_trim(strsplit(tbl[["V6"]][x], ",")[[1]])
  sapply(ncs, function(one_nc){
    ifelse(any(grepl(one_nc, nc_have)),
           (fread(paste0("nc_pulls/",nc_have[grep(one_nc, nc_have)])) %>%  filter(n_reads >= 100 & score > 50) %>% dim())[1],
           0)
  }) %>% sum()
})

write.table(tbl[,c(1,2,6,7,8)], file = "../output/Serratus_hits_all_viruses.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
tbl %>% filter(reactivation_candidate)
tbl
