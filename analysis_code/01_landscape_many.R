library(data.table)
library(dplyr)
library(stringr)
library(BuenColors)

# Import SRA meta data
sra <- fread("../meta/hu_SraRunInfo.csv.gz")
all_hits <- fread("../output/Serratus_libraries_all_hits_reactivation.tsv")
all_viral_df <- merge(sra, all_hits, by.y = "run_id", by.x = "Run")

# Now assemble data
small_table_for_supplement <- merge( all_hits, sra[,c("Run", "BioSample")], by.x = "run_id", by.y = "Run", all.x = TRUE) %>%
  filter()
small_table_for_supplement <- small_table_for_supplement[!is.na(small_table_for_supplement$BioSample),]
small_table_for_supplement <- small_table_for_supplement %>% filter(score >=50)
dim(small_table_for_supplement)

# Now summarize and collapse 
small_table_for_supplement %>%
  group_by(BioSample, viral_name) %>% summarize(count = n(), nBioSample = 1) %>%
  group_by(viral_name) %>% summarize(countLibrary = sum(count), countBioSample = sum(nBioSample)) -> new_count_df
sum(new_count_df$count)

sum(new_count_df$countBioSample)
sum(new_count_df$countLibrary)

# Make basic plot
cc <- as.character(new_count_df$viral_name )
new_count_df$viral_name <- factor(cc, levels = rev(cc))
new_count_df

if(FALSE){
  write.table(new_count_df, file = "../output/source_data_bar_graph.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  write.table(small_table_for_supplement, file = "../output/all_libraries_reactivation_annotated.tsv", 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Make plot
ggplot(new_count_df, aes(x = viral_name, y = countBioSample)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "count of human RNA-seq BioSamples") +
  coord_flip()-> p1
p1
cowplot::ggsave2(p1, file = "../output/active_viruses_human_samples.pdf", width = 3, height = 1.5)

