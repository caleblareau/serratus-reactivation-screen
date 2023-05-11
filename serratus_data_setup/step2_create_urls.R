library(data.table)
library(dplyr)
library(stringr)

# Make the suite of URLs to pull data from
tbl <-  fread("Table_human_viruses.txt", header = FALSE)
tbl <- tbl[tbl$V6 != "Not available",]
nc_vec <-tbl[["V6"]]

# Loop over all possibilities
lapply(nc_vec, function(x){
  ncs <- str_trim(strsplit(x, ",")[[1]])
  lapply(ncs, function(one_nc){
    lapply(1:9, function(i){
      me <- paste0(one_nc, ".", i)
      data.frame(x = paste0("wget -O ",me, '.txt "https://api.serratus.io/matches/nucleotide?scoreMin=0&scoreMax=100&identityMin=75&identityMax=100&sequence=',me,'"'))
    }) %>% rbindlist() 
  }) %>% rbindlist() 
}) %>% rbindlist() -> all_wget

write.table(all_wget, file = "step3_run_wget_cmds.sh", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Redo logic to see what's missting
tbl <- fread("Table_human_viruses.txt", header = FALSE)
tbl <- tbl[tbl$V6 != "Not available",]
lapply(1:dim(tbl)[1], function(x){
  ncs <- str_trim(strsplit(tbl[["V6"]][x], ",")[[1]])
  lapply(ncs, function(one_nc){
    tbl$what <- one_nc
    tbl[x,]
  }) %>% rbindlist() 
}) %>% rbindlist() -> all_annotations

length(unique(all_annotations$V1))
dim(all_annotations)

nc_have <- list.files("nc_pulls/",pattern = "^NC_")
nc_have_simple <- sapply(strsplit(nc_have, "[.]"), function(x) x[1])

all_annotations %>% filter(!(what %in% nc_have_simple))


