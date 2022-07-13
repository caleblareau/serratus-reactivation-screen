library(data.table)
library(dplyr)

sra <- fread("../data/hu_SraRunInfo.csv.gz")
virus <- c("HHV6B")
pull_viral_info <- function(which_virus){
  virus <- fread(paste0("../data/manual/SerratusMatches-",which_virus,".csv.gz"))
  human_reactivation <- merge(virus, sra, by.x = "run_id", by.y = "Run")
  
  human_reactivation_filt <- human_reactivation %>%
    filter(n_reads >= 100 & score > 50) %>%
    arrange(desc(n_reads)) %>%
    mutate(virus_name = which_virus)
  human_reactivation_filt
}

lapply(virus, pull_viral_info) %>%
  rbindlist() %>% data.frame() -> all_viral_df

# Filter for GEO samples
all_viral_df_geo <- all_viral_df[grepl("GSM", all_viral_df$SampleName),]
all_viral_df[!grepl("GSM", all_viral_df$SampleName),] %>%
  filter(virus_name == "HHV6B")

library(GEOquery)

# Manually filter problematic ones
blacklist <- c("GSM4467089", "GSM4467088", "GSM4467090", "GSM4467091")
all_viral_df_geo <- all_viral_df_geo[!(all_viral_df_geo$SampleName %in% blacklist),]

all_viral_df_geo_get <- all_viral_df_geo[!(all_viral_df_geo$SampleName%in% gsub(".soft", "", list.files("../metadata/"))),]
dim(all_viral_df_geo_get)
titles_df <- lapply(all_viral_df_geo$SampleName, function(gsm){
  print(gsm)
  gds <- GEOquery::getGEO(gsm, destdir = "../metadata/")
  data.frame(gsm, 
             title = gds@header$title, 
             name1 = gds@header$source_name_ch1
  )
}) %>% rbindlist() %>% data.frame()


mdf <- merge(all_viral_df_geo, titles_df, by.x = "SampleName", by.y = "gsm")[,c("virus_name","run_id", "SampleName", "score", "n_reads", "title", "name1")] %>%
  arrange(desc(n_reads))

hhv6 <- mdf[grepl("T-cell|\ T\ cell|Tcell", mdf$name1) | grepl("T-cell|\ T\ cell|Tcell", mdf$title) ,] %>% distinct() %>%
  arrange(desc(score)) 

# read back in the T cells
others <- fread("../output/Tcell_hits_per_virus.tsv")[,-8] %>%
  filter(virus_name != "Human herpesvirus 6") %>% distinct
  

total_tcell_df <- rbind(hhv6, others) %>% distinct() %>% arrange(desc(score))
total_tcell_df
table(total_tcell_df$virus_name)
write.table(total_tcell_df, file = "../output/Tcell_Table-wHHV6B.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
