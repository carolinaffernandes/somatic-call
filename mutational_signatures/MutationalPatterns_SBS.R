# ---- Load libraries ----
library(MutationalPatterns)
library(BSgenome)
library(GenomicRanges)
library(gridExtra)
library(NMF)
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
cat("Conferindo o genoma de referência:", ref_genome, "\n")

setwd(work_dir)

group_file_name <- fread("groups.txt", header = FALSE, sep = "\t")
colnames(group_file_name) <- c("groups", "sample")

samples_per_group <- group_file_name %>%
  group_by(groups) %>%
  summarise(count = n_distinct(sample))
print(samples_per_group)

qtde_grupos <- n_distinct(group_file_name$groups)
cat("Quantidade de grupos:", qtde_grupos, "\n")


sample_names <- group_file_name$sample


sample_id <- sub(".hg38_multianno.txt", "", sample_names)


groups <- group_file_name$groups

myfiles <- lapply(sample_names, function(x) {
  fread(x, header = TRUE, sep = "\t", colClasses = c("Ref" = "character", "Alt" = "character"))
})


sample_gr <- list()

for (i in seq_along(sample_names)) {
  # Selecionar colunas necessárias (assumindo 1-5 = chr, start, end, ref, alt)
  file <- as.data.frame(myfiles[[i]][, 1:5])
  colnames(file) <- c("chr", "startvar", "endvar", "ref", "var")

  file$sample_names <- sample_names[i]
  file$sampleID <- sample_id[i]
  file$group <- groups[i]
  
  gr <- GRanges(
    seqnames = file$chr,
    ranges = IRanges(file$startvar, file$endvar),
    ref = file$ref,
    alt = file$var,
  )
  
  mcols(gr)$group <- file$group
  mcols(gr)$sampleID <- file$sampleID
  mcols(gr)$sample_names <- file$sample_names
  
  sample_gr[[i]] <- gr
  names(sample_gr)[i] <- sample_id[i]
}

cat("Objetos GRanges criados para todas as amostras.\n")

seqinfo_hg38 <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
seqinfo(sample_gr[[1]]) <- seqinfo_hg38[names(seqinfo(sample_gr[[1]]))]
for (i in seq_along(sample_gr)) {
  seqinfo(sample_gr[[i]]) <- seqinfo_hg38[names(seqinfo(sample_gr[[i]]))]
}

mut_mat <- mut_matrix(vcf_list = sample_grList, ref_genome = BSgenome.Hsapiens.UCSC.hg38)

type_occurrences <- mut_type_occurrences(sample_grList, ref_genome)
type_occurrences
ps0 <- plot_spectrum(type_occurrences, by = groups, CT = TRUE, legend = TRUE, error_bars = 'none')
ggsave("type_occurrences.png", ps0, width=30, height=20, units='cm', dpi=300)

plot96 <- plot_96_profile(mut_mat, condensed = FALSE, ymax = 0.1)
ggsave("plot_96_profile.png", plot96, width=30, height=20, units='cm', dpi=300)


genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_hg38

mut_mat_strand <- mut_matrix_stranded(
  vcf_list = sample_grList,
  ref_genome = BSgenome.Hsapiens.UCSC.hg38, genes_hg38, mode = "transcription")

head(mut_mat_strand)


plot192 <- plot_192_profile(mut_mat_strand, condensed = FALSE, ymax = 0.1)
ggsave("plot_192_profile.png", plot192, width=30, height=20, units='cm', dpi=300)

plot_96TR<- grid.arrange(plot96, plot192)
ggsave("plot_96e192.png", plot_96TR, width=20, height=20, units='cm', dpi=300)

strand_counts <- strand_occurrences(mut_mat_strand, by = groups)
strand_counts

strand_bias <- strand_bias_test(strand_counts)
strand_bias

write.table(strand_bias, file = "strand_bias.txt", row.names = FALSE)

strand_bias_notstrict <- strand_bias_test(strand_counts,
                                          p_cutoffs = 0.05,
                                          fdr_cutoffs = 0.2
)

strand_bias_notstrict

write.table(strand_bias_notstrict, file = "strand_bias.txt", row.names = FALSE)

ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand(strand_counts, mode = "absolute")
ps3 <- plot_strand_bias(strand_bias_notstrict, sig_type = "fdr")

plot_strandbias <- grid.arrange(ps1, ps2, ps3)
ggsave("plot_strandbias.png", plot_strandbias, width=20, height=25, units='cm', dpi=300)

cosmic_signatures <- get_known_signatures(muttype = "snv", source="COSMIC", genome = "GRCh38")

fit_res <- fit_to_signatures(mut_mat, cosmic_signatures)
plot_contribution(fit_res$contribution, coord_flip = FALSE)


contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               cosmic_signatures,
                                               n_boots = 1000,
                                               method = "strict"
)
plot_bootstrapped1 <- plot_bootstrapped_contribution(contri_boots)
plot_bootstrapped2 <- plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")
plot_bootstrapped1
plot_bootstrapped2

ggsave("bootstraped_absolute.png", plot_bootstrapped1, width = 30, height = 25, units = 'cm', dpi = 300)
ggsave("bootstraped_relative.png", plot_bootstrapped2, width = 30, height = 25, units = 'cm', dpi = 300)

plot_bootstrapped3 <- plot_bootstrapped_contribution(contri_boots, mode = "absolute", plot_type = "barplot")
plot_bootstrapped3

selected_sig1 <- which(rowSums(fit_res$contribution) > 0)
plot_contribution(fit_res$contribution[selected_sig1,],
                  coord_flip = FALSE,
                  mode = "relative")
selected_sig1

strict_refit = fit_to_signatures_strict(mut_mat, cosmic_signatures, max_delta = 0.0016)
fit_res_strict <- strict_refit$fit_res
strict_refit

selected_sig2 <- which(rowSums(fit_res_strict$contribution) > 0)
prelative <- plot_contribution(fit_res_strict$contribution[selected_sig2,],
                  coord_flip = FALSE,
                  mode = "relative")
selected_sig2
pabsolute <- plot_contribution(fit_res_strict$contribution[selected_sig2,],
                  coord_flip = FALSE,
                  mode = "absolute")

ggsave("SBS-2.png", dpi = 600, width=15, height=17, unit = 'cm', plot = prelative)
ggsave("SBS-2-absolute.png", dpi = 600, width=15, height=17, unit = 'cm', plot = pabsolute)

library(patchwork)

combined_plot <- prelative | pabsolute  # empilha horizontal; / se vertical
combined_plot
ggsave("SBS-2-combined.png", dpi = 600, width = 30, height = 17, unit = "cm",
       plot = combined_plot)
