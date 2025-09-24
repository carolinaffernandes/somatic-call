library(data.table)
library(gridExtra)
library(NMF)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(MutationalPatterns)
library(GenomicRanges)
library(VariantAnnotation)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
seqinfo_hg38 <- seqinfo(get(ref_genome))
library(ref_genome, character.only = TRUE)

vcf_dir <- setwd("/home/carolis/2025/Teste/vcf_dir")
group_file_name <- fread(file.path(vcf_dir, "groups.txt"), header = FALSE, sep = "\t")
colnames(group_file_name) <- c("groups", "sample")

samples_per_group <- group_file_name %>%
  group_by(groups) %>%
  summarise(count = n_distinct(sample))
print(samples_per_group)

sample_names <- group_file_name$sample
sample_id <- sub(".hg38_multianno.txt", "", sample_names)
groups <- group_file_name$groups

myfiles <- lapply(sample_names, function(x) {
  fread(x, header = TRUE, sep = "\t", colClasses = c("Ref" = "character", "Alt" = "character"))
})

##### PER SAMPLE -------------------------------------------------------------------------------

sample_gr <- list()

for (i in seq_along(sample_names)) {
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
    seqinfo = seqinfo_hg38
  )
  
  mcols(gr)$group <- file$group
  mcols(gr)$sampleID <- file$sampleID
  mcols(gr)$sample_names <- file$sample_names
  
  sample_gr[[i]] <- gr
  names(sample_gr)[i] <- sample_id[i]
}

cat("Objetos GRanges criados para todas as amostras.\n")
sample_grList <- GRangesList(sample_gr)
sample_grList <- GRangesList(sample_gr)
names(sample_grList) <- sample_id


##### PER GROUP  -------------------------------------------------------------------------------

group_gr <- list()
unique_groups <- unique(groups)

for (g in unique_groups) {
  idx <- which(groups == g)
  file <- do.call(rbind, lapply(idx, function(i) {
    df <- as.data.frame(myfiles[[i]][, 1:5])
    colnames(df) <- c("chr", "startvar", "endvar", "ref", "alt")
    df$sampleID <- sample_id[i]
    df$group <- groups[i]
    df
  }))
  
  gr <- GRanges(
    seqnames = file$chr,
    ranges = IRanges(file$startvar, file$endvar),
    ref = file$ref,
    alt = file$alt,
    seqinfo = seqinfo_hg38
  )
  
  mcols(gr)$sampleID <- file$sampleID
  mcols(gr)$group <- file$group
  
  group_gr[[g]] <- gr
}

groups_grList <- GRangesList(group_gr)
names(groups_grList) <- as.vector(unique(groups))

##CHECK TIME ---------------------------------------------------------------------------------

resumo_amostras <- tibble(
  sampleID = names(sample_grList),
  grupo = sapply(sample_grList, function(gr) unique(mcols(gr)$group)),
  registros = lengths(sample_grList)
)
resumo_amostras
resumo_grupos <- tibble(
  grupo = names(groups_grList),
  registros = lengths(groups_grList)  
)
resumo_grupos
soma_por_grupo <- resumo_amostras %>%
  group_by(grupo) %>%
  summarise(registros_calculados = sum(registros))
conferencia <- resumo_grupos %>%
  left_join(soma_por_grupo, by = "grupo") %>%
  mutate(
    igual = registros == registros_calculados,
    status = ifelse(igual, "OK", "Erro")
  )

cat("\nConferência final (GRangesList x soma de amostras):", conferencia$status)

##### PER GROUPS -------------------------------------------------------------------------------

## ----mut_type_occurrences---------------------------------------------------------------------
type_occurrences_groups <- mut_type_occurrences(sample_grList, ref_genome)
type_occurrences_groups

## ----plot_spectrum--------------------------------------------------------------------------
plot_type_occurrences_groups <- plot_spectrum(type_occurrences_groups, CT = TRUE, by = groups, 
                                              indv_points = TRUE)
plot_type_occurrences_groups


##### PER SAMPLES -------------------------------------------------------------------------------

## ----mut_type_occurrences---------------------------------------------------------------------
type_occurrences_samples <- mut_type_occurrences(sample_grList, ref_genome)
type_occurrences_samples

## ----plot_spectrum--------------------------------------------------------------------------
plot_type_occurrences_samples <- plot_spectrum(type_occurrences_samples, CT = TRUE, 
                                               by = sample_id)
plot_type_occurrences_samples

## ----mut_matrix-------------------------------------------------------------------------------
mut_mat <- mut_matrix(vcf_list = sample_grList, ref_genome = ref_genome)
head(mut_mat)

## ----plot_96_profile, fig.wide=TRUE-----------------------------------------------------------
plot_96_profile(mut_mat)

## ----mut_mat_extendend_context----------------------------------------------------------------
mut_mat_ext_context <- mut_matrix(sample_grList, ref_genome, extension = 2)
head(mut_mat_ext_context)

## ----plot_profile_heatmap, fig.wide=TRUE------------------------------------------------------
plot_profile_heatmap(mut_mat_ext_context, by= groups)
plot_river(mut_mat_ext_context)
