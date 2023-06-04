library(dplyr)
library(data.table)
library(BuenColors)

# Ipmort T cell runs
t_cell_runs <- fread("../output/all_tcell_SRAs.tsv") %>% pull(Run)

# Specially compare HHV-6
hhv6b <- fread(paste0("../serratus_data_setup/manual/SerratusMatches-HHV6B.csv.gz"))
hhv6a <- fread(paste0("../serratus_data_setup/nc_pulls/NC_001664.4.txt"))

a_tcells <- hhv6a %>%
  filter(run_id %in% t_cell_runs) %>%
  filter(score >= 50 & n_reads >= 100) %>% pull(run_id)

b_tcells <- hhv6b %>%
  filter(run_id %in% t_cell_runs) %>%
  filter(score >= 50 & n_reads >= 100) %>% pull(run_id)
a_tcells %in% b_tcells

merged_hhv6 <- merge(hhv6b %>%
        filter(run_id %in% t_cell_runs) %>%
        filter(score >= 50 & n_reads >= 100),
      hhv6a, all.x = TRUE, by = "run_id") %>%
  arrange(desc(n_reads.x)) %>%
  mutate(rank = 1:n()) 

# Now plot
subset_merge_hhv6 <- merged_hhv6[,c("rank", "n_reads.x", "n_reads.y")]; colnames(subset_merge_hhv6) <- c("rank", "HHV6B", "HHV6A")
subset_merge_hhv6 %>% reshape2::melt(id.vars = "rank") -> melted_plot_df
p1 <- melted_plot_df %>%
  ggplot(aes(x = rank, y = value, color = variable)) + 
  geom_line(data = melted_plot_df, aes(x = rank, y = value, group = as.character(rank)), color = "black")+
  geom_point() + scale_y_log10() +
  scale_x_continuous(breaks = c(1, 5, 10, 15)) +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") + 
  scale_color_manual(values = c("firebrick", "dodgerblue3")) +
  labs(x = "Rank-ordered libraries", y = "# reads assigned to HHV-6 genome")
cowplot::ggsave2(p1, file = "../output/hhv6_strain_comparison.pdf", width = 1.8, height=1.8)
