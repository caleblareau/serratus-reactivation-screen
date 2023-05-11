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
           (fread(paste0("nc_pulls/",nc_have[grep(one_nc, nc_have)])) %>%  filter(n_reads >= 100 & score > 90) %>% dim())[1],
           0)
  }) %>% sum()
})
odf <- tbl[,c(1,2,6,7,9)]
odf %>% filter(reactivation_candidate)
colnames(odf) <- c("virus_name", "family", "genomes", "reactivation_candidate", "serratus_high_confidence_hits")

write.table(odf, file = "../output/Serratus_hits_all_viruses.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# get each entry per virus
ncs_reactive <- tbl %>% filter(reactivation_candidate) %>% pull(V6)
common_names <- tbl%>% filter(reactivation_candidate) %>% pull(V1)
lapply(1:length(ncs_reactive), function(i){
  one_nc <- ncs_reactive[i]
  fread(paste0("nc_pulls/",nc_have[grep(one_nc, nc_have)])) %>%  filter(n_reads >= 100 & score > 50) %>%
    mutate(viral_name = common_names[i])
}) %>% rbindlist() -> full_df
dim(full_df)
write.table(full_df, file = "../output/Serratus_libraries_all_hits_reactivation.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)


