library(data.table)
library(dplyr)

# Import HHV-6B and annotate accordingly
sra <- fread("../meta/hu_SraRunInfo.csv.gz")
hhv6b <- fread(paste0("../serratus_data_setup/manual/SerratusMatches-HHV6B.csv.gz"))  %>%
  filter(score >= 50 & n_reads >= 100) %>%
  mutate(viral_name = "Human herpesvirus 6B")
hhv6b_annotated <- data.frame(merge(hhv6b, sra, by.x = "run_id", by.y = "Run"))
sample_name_vec <- sra$SampleName; names(sample_name_vec) <- sra$Run

# Import celltypes to filter
tcell_df <- fread("../output/all_tcell_SRAs.tsv")
ipsc_df <- fread("../output/all_iPSC_SRAs.tsv")

# First plot the proportion of T cells -- 
all %>%
  mutate(is_t_cell= run_id %in% tcell_df$Run) %>% 
  group_by(viral_name) %>% summarize(count = n(), sum_t = sum(is_t_cell)) %>%
  mutate(perc_tcell = sum_t/count*100) %>% filter(sum_t > 0) -> perc_plot_df

cc <- unique(perc_plot_df$viral_name)
perc_plot_df$viral_name <- factor(cc, levels = rev(cc))
perc_plot_df %>%
  ggplot(aes(x = viral_name, y = perc_tcell)) + 
  geom_bar(fill = "lightgrey", color = "black", stat = "identity") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10)) +
  labs(x = "", y = "") +
  coord_flip() -> pO
pO
cowplot::ggsave2(pO, file = "../output/percent_Tcells.pdf", width = 2.2, height = 1.5)



# Combine with Viralzone queries
all <- fread(paste0("../output/all_libraries_reactivation_annotated.tsv"))
combined_df <- rbind(
  hhv6b_annotated[,colnames(all)],
  all
)
combined_df$SampleName <- sample_name_vec[combined_df$run_id]

# see what are in geo
table(grepl("GSM", combined_df$SampleName))

# annotate for utility
combined_df <- combined_df[grepl("GSM", combined_df$SampleName),]
bl <- c("GSM4467089", "GSM4467088","GSM4467090", "GSM4467091") # problematic samples to query
combined_df <- combined_df[!(combined_df$SampleName %in% bl),]

titles_df <- lapply(combined_df$SampleName, function(gsm){
  gds <- GEOquery::getGEO(gsm, destdir = "../metadata/")
  data.frame(gsm, 
             title = gds@header$title, 
             name1 = gds@header$source_name_ch1
  )
}) %>% rbindlist() %>% distinct() %>% data.frame()
combined_df$title <- combined_df[["title"]]
combined_df$name1 <- combined_df[["name1"]]




