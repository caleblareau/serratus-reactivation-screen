library(data.table)
library(dplyr)
library(stringr)

# Import SRA meta data
sra <- fread("../meta/hu_SraRunInfo.csv.gz")

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

