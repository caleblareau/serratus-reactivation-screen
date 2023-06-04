library(data.table)
library(dplyr)
library(BuenColors)

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

