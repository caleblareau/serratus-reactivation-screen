library(dplyr)
library(data.table)
library(GEOquery)

# Import SRA meta data
sra <- fread("../meta/hu_SraRunInfo.csv.gz")
all_hits <- fread("../output/Serratus_libraries_all_hits_reactivation.tsv")
all_viral_df <- merge(sra, all_hits, by.y = "run_id", by.x = "Run")
all_viral_df_geo <- all_viral_df[grepl("^GSM", all_viral_df$SampleName),]

# Deeper dive with GEO meta data
all_viral_df_geo_get <- all_viral_df_geo[!(all_viral_df_geo$SampleName%in% gsub(".soft", "", list.files("../metadata/"))),]
dim(all_viral_df_geo_get)

# Black list a few samples that cause issues with the geo query
bl <- c("GSM4467089", "GSM4467088","GSM4467090", "GSM4467091")
all_viral_df_geo <- all_viral_df_geo[!c(all_viral_df_geo$SampleName %in% bl),]
titles_df <- lapply(all_viral_df_geo$SampleName, function(gsm){
  gds <- GEOquery::getGEO(gsm, destdir = "../metadata/")
  data.frame(gsm, 
             title = gds@header$title, 
             name1 = gds@header$source_name_ch1
  )
}) %>% rbindlist() %>% distinct() %>% data.frame()

# merge Serratus data with GEO meta data to find enriched patterns
mdf <- merge(all_viral_df_geo, titles_df, by.x = "SampleName", by.y = "gsm")[,c("viral_name","Run", "SampleName", "score", "n_reads", "title", "name1")] %>%
  arrange(desc(n_reads))

mdf$known_infection <- grepl("dpi|nfect|hpi|p.i.", mdf$title) | grepl("dpi|nfect|hpi|p.i.", mdf$name1)

# Look at the top hits per virus
mdf[ complete.cases(mdf)  ,] %>%
  distinct() %>%
  group_by(viral_name) %>% 
  filter(!known_infection) %>%
  top_n(5, wt = n_reads)  %>%
  arrange((viral_name)) %>%
  data.frame() %>%
  write.table(file = "../output/top_hits_per_virus.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


