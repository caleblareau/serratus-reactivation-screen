library(data.table)
library(dplyr)
library(stringr)

# Import SRA meta data
sra <- fread("../meta/hu_SraRunInfo.csv.gz")
all_hits <- fread("../output/Serratus_libraries_all_hits_reactivation.tsv")
all_viral_df <- merge(sra, all_hits, by.y = "run_id", by.x = "Run")

small_table_for_supplement <- merge( all_hits, sra[,c("Run", "BioSample")], by.x = "run_id", by.y = "Run", all.x = TRUE) %>%
  filter(score > 90)
small_table_for_supplement[!is.na(small_table_for_supplement$BioSample),]
dim(small_table_for_supplement)
small_table_for_supplement %>% filter(run_id %in% c("ERR3240275", "ERR3240276"))
small_table_for_supplement[,c("BioSample", "viral_name")] %>%
  distinct() %>%
  group_by(viral_name) %>% summarize(count = n()) -> new_count_df
sum(new_count_df$count)

new_count_df

# Make basic plot
library(BuenColors)
cc <- as.character(new_count_df$viral_name )
new_count_df$viral_name <- factor(cc, levels = rev(cc))
write.table(new_count_df, file = "../output/source_data_bar_graph.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Make plot
ggplot(new_count_df, aes(x = viral_name, y = count)) +
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "count of human RNA-seq BioSamples") +
  coord_flip()-> p1
p1
#cowplot::ggsave2(p1, file = "../output/active_viruses_human_samples.pdf", width = 3, height = 1.5)


# Filter for annotated iPSCs
ipsc_biosamples <- unique(fread("../meta/access12March2023/ipsc-biosample_result.txt", header = FALSE)[[1]])
ipsc_samples <- unique(fread("../meta/access12March2023/ipsc-sra-query.csv.gz", header = TRUE)[["Sample Accession"]])

ipsc_sras <- sra %>%
  filter(Sample %in% ipsc_samples | BioSample %in% ipsc_biosamples)

# Do the same for T cells
tcell_biosamples <- unique(fread("../meta/access12March2023/tcell-biosample-query.txt", header = FALSE)[[1]])
tcell_samples <- unique(fread("../meta/access12March2023/tcell-sra-query.csv.gz", header = TRUE)[["Sample Accession"]])

tcell_sras <- sra %>%
  filter(Sample %in% tcell_samples | BioSample %in% tcell_biosamples)

# Export all iPSC/ T cell SRAs
write.table(tcell_sras,file = "../output/all_tcell_SRAs.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ipsc_sras,file = "../output/all_iPSC_SRAs.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### other misc
all_hits %>% filter(run_id %in% tcell_sras$Run) %>%
  filter(score > 90)

fread("../serratus_data_setup/manual/SerratusMatches-HHV6B.csv.gz") %>%
  filter(score >= 50) %>%
  filter(run_id %in% tcell_sras$Run)

