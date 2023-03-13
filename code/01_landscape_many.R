library(data.table)
library(dplyr)
library(stringr)

# Import SRA meta data
sra <- fread("../data/hu_SraRunInfo.csv.gz")
nc_virus <- gsub(".txt", "", list.files("../data/nc_pulls/"))

pull_viral_info <- function(which_virus){
  virus <- fread(paste0("../data/nc_pulls/",which_virus,".txt"))
  human_reactivation <- merge(virus, sra, by.x = "run_id", by.y = "Run")
  
  human_reactivation_filt <- human_reactivation %>%
    filter(n_reads >= 100 & score > 50) %>%
    arrange(desc(n_reads)) %>%
    mutate(which_virus = which_virus)
  human_reactivation_filt
}

lapply(nc_virus, pull_viral_info) %>%
  rbindlist() %>% data.frame() -> all_viral_df
table(all_viral_df$which_virus)
all_viral_df$which_virus_short <- sapply(strsplit(all_viral_df$which_virus, "[.]"), function(x) x[1])

# Import virus annotation
tbl <- fread("../data/Table_human_viruses.txt", header = FALSE)

# Grep for those that can be latent / reactivated
tbl <- tbl[grepl(c("Herpesviridae|Polyomaviridae|Parvoviridae|Adenoviridae"),tbl$V2),]
tbl
tbl <- tbl[tbl$V6 != "Not available",]

#-----
# Chunk: annotate viral name
lapply(1:dim(tbl)[1], function(x){
  ncs <- str_trim(strsplit(tbl[["V6"]][x], ",")[[1]])
  lapply(ncs, function(one_nc){
    tbl$what <- one_nc
    tbl[x,]
  }) %>% rbindlist() 
}) %>% rbindlist() -> all_annotations

virus_vec <- all_annotations[["V1"]]; names(virus_vec) <- all_annotations[["what"]]
all_viral_df$virus_name <- virus_vec[as.character(all_viral_df$which_virus_short)]

all_viral_df[!is.na(all_viral_df$virus_name),c("virus_name", "BioSample")] %>%
  distinct() %>%
  group_by(virus_name) %>% summarize(count = n()) -> new_count_df

# Make basic plot
library(BuenColors)
cc <- as.character(new_count_df$virus_name )
new_count_df$virus_name <- factor(cc, levels = rev(cc))

write.table(new_count_df, file = "../output/source_data_bar_graph.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Make plot
ggplot(new_count_df, aes(x = virus_name, y = count)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "count of human RNA-seq BioSamples") +
  coord_flip()-> p1
p1
cowplot::ggsave2(p1, file = "../output/active_viruses_human_samples.pdf", width = 3, height = 1.5)

# Filter for GEO samples
all_viral_df_geo <- all_viral_df[grepl("GSM", all_viral_df$SampleName),]

library(GEOquery)

# Deeper dive with GEO meta data
length(unique(all_viral_df_geo$SampleName))
table(all_viral_df_geo$virus_name)

all_viral_df_geo_get <- all_viral_df_geo[!(all_viral_df_geo$SampleName%in% gsub(".soft", "", list.files("../metadata/"))),]
dim(all_viral_df_geo_get)

# Black list a few samples that cause issues
bl <- c("GSM4467089", "GSM4467088","GSM4467090", "GSM4467091")
all_viral_df_geo <- all_viral_df_geo[!c(all_viral_df_geo$SampleName %in% bl),]
titles_df <- lapply(all_viral_df_geo$SampleName, function(gsm){
  print(gsm)
  gds <- GEOquery::getGEO(gsm, destdir = "../metadata/")
  data.frame(gsm, 
             title = gds@header$title, 
             name1 = gds@header$source_name_ch1
  )
}) %>% rbindlist() %>% data.frame()

# merge Serratus data with GEO meta data to find enriched patterns
mdf <- merge(all_viral_df_geo, titles_df, by.x = "SampleName", by.y = "gsm")[,c("virus_name","run_id", "SampleName", "score", "n_reads", "title", "name1")] %>%
  arrange(desc(n_reads))
data.frame(sort(table(mdf$name1))) %>% arrange(desc(Freq)) %>% head()

mdf$known_infection <- grepl("dpi|nfect|hpi|p.i.", mdf$title) | grepl("dpi|nfect|hpi|p.i.", mdf$name1)

# Make T cell count
tcell_string <- "T-cell|\ T\ cell|Tcell|CD4T|CD8T|^T cell|^T-Lymph|T\ lymph"
mdf[complete.cases(mdf) & (grepl(tcell_string, mdf$name1) | grepl(tcell_string, mdf$title)),] %>% distinct()
mdf[grepl(tcell_string, mdf$name1) | grepl(tcell_string, mdf$title) ,] %>% distinct() 

# See all HHV6 hits
mdf %>% filter(virus_name == "Human herpesvirus 6")

# Look at the top hits per virus
mdf[ complete.cases(mdf)  ,] %>%
  distinct() %>%
  group_by(virus_name) %>% 
  filter(!known_infection) %>% top_n(5, wt = n_reads)  %>%
  arrange((virus_name)) %>%
  data.frame() %>%
  write.table(file = "../output/top_hits_per_virus.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
# Export T cells
mdf[complete.cases(mdf) & (grepl("T-cell|\ T\ cell|Tcell|CD4T|CD8T|^T cell", mdf$name1) | grepl("T-cell|\ T\ cell|Tcell|CD4T|CD8T |^T cell", mdf$title)),] %>% distinct() %>%
  group_by(virus_name) %>%
  write.table(file = "../output/Tcell_hits_per_virus.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# export iPSC cells
mdf[(grepl("iPSC|luripotent", mdf$name1) | grepl("iPSC|luripotent", mdf$title)) & complete.cases(mdf),] %>% distinct() %>%
  write.table(file = "../output/iPSC_hits_per_virus.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# Other exporatory analyses-- didn't find much
mdf[grepl("lood|PBMC", mdf$name1) ,] %>% distinct() %>%
  group_by(virus_name) %>% summarize(count = n())
mdf[grepl("lood|PBMC", mdf$name1) ,] %>% distinct() 

mdf %>% filter(virus_name %in% c(  "Human cytomegalovirus"))

mdf %>% filter(virus_name %in% c("JC polyomavirus")) 

(mdf[!(grepl("GM12|lymphoblastoid|LCL|EBV|Epstein|Lymphob|Raji", mdf$title) | grepl("GM12|lymphoblastoid|LCL|EBV|Epstein|Lymphob|Raji|LC_control", mdf$name1)),] %>% 
    filter(which_virus %in% c("EBV")))[150:300,]

(mdf[!(grepl("dpi|nfect|hpi", mdf$title) | grepl("dpi|nfect|hpi", mdf$name1)),] %>% filter(which_virus %in% c("CMV")) )
