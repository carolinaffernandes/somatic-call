# carrega bibliotecas
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(data.table)
library(ggplot2)

setwd(work_dir)

group_file_name = read.csv("groups.txt", header = FALSE, sep = "\t")
colnames(group_file_name) = c("groups", "sample")

samples_per_group = as.data.frame(setDT(group_file_name)[, .(count = uniqueN(sample)), by = groups])
samples_per_group

qtde_grupos = length(samples_per_group$groups)
qtde_grupos

sample_names = group_file_name$sample
sample_names

groups = group_file_name$groups
groups

sample_id = sub(".hg38_multianno.txt", "", sample_names)
sample_id

myfiles = lapply(sample_names, read.table, header=TRUE, dec=",", sep = "\t", colClasses = c("Ref" = "character",
                                                                                            "Alt" = "character"))
listfile=list() #empty list for output 

for(i in 1:length(sample_names))
{
    # transforma cada arquivo em dataframe
    file=as.data.frame(myfiles[i]) 
    # filtra apenas as colunas que vou precisar
    file=file[,c(1:5)]
    # renomeia as colunas
    colnames(file) = c("chr", "startvar", "endvar", "ref", "var")
    
    # adiciona a coluna com o nome do arquivo (amostra)
    file$sample_names = sample_names[i]
    # adiciona a coluna com o ID da amostra (amostra)
    file$sample_id = sample_id[i]
    # adiciona coluna com o nome do grupo (amostra também)
    file$group = groups[i]
    
    
    listfile[[i]] = file # armazena numa lista
}


merged.dataframe=do.call("rbind", listfile)

merged_vr = VRanges(seqnames = merged.dataframe$chr,
                    ranges = IRanges(merged.dataframe$startvar,merged.dataframe$endvar),
                    ref = merged.dataframe$ref,
                    alt = merged.dataframe$var,
                    study = merged.dataframe$group,
                    sampleNames = merged.dataframe$sample_id)
head(merged_vr)

merged_motifs = mutationContext(merged_vr, BSgenome.Hsapiens.UCSC.hg38)
plotMutationSpectrum(merged_motifs, "study")


merged_mm = motifMatrix(merged_motifs, group = "study", normalize = TRUE)

rowSums(merged_mm, )
rowSums(merged_mm == 0): verifica quantas colunas da matriz tem valor == 0
# >> se for apenas 1 coluna, mant?m a linha
# >> se forem as 2 columas, remove a linha 
merged_mm = merged_mm[rowSums(merged_mm == 0) != ncol(merged_mm),]

gof_nmf = assessNumberSignatures(merged_mm, nSigs = 2, nReplicates = 30)

plotNumberSignatures(gof_nmf)

write.table(gof_nmf, "RSS_ExplainedVariance_numberNovelSignatures.txt", quote = FALSE)

# RSS: Position-dependent motif characterization using non-negative matrix factorization
ggsave("RSS_ExplainedVariance_numberNovelSignatures.svg", dpi = 600, width=10, height=15, units = 'cm')
ggsave("RSS_ExplainedVariance_numberNovelSignatures.png", dpi = 600, width=10, height=15, units = 'cm')


#gof_pca = assessNumberSignatures(merged_mm, n_sigs, pcaDecomposition)
#plotNumberSignatures(gof_nmf)

nSigs = n_sigs_nmf
sigs_nmf = identifySignatures(merged_mm, nSigs, nmfDecomposition)
sigs_nmf@signatures
sigs_nmf@samples
sigs_nmf@fitted


rm(pHeatMap)
pHeatMap = plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
pHeatMap = pHeatMap + scale_fill_gradient(low="#F9EEFC", high="#2F015A")
pHeatMap

ggsave("SomaticSig_SignatureMap_heatmap.svg", dpi = 600, width=10, height=20, units = 'cm')
ggsave("SomaticSig_SignatureMap_heatmap.png", dpi = 600, width=10, height=20, units = 'cm')

rm(p)

p = plotSignatures(sigs_nmf)
p = p + xlab("")
p = p + theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major=element_blank(),
              
              axis.ticks = element_line(color = 'black', linewidth = 0.2),
              axis.title.y = element_text(size = 11, color = 'black'),
              axis.text.y = element_text(size = 8, color = 'black'),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 7, color = 'black'),
              
              # altera tamanho dos labels (facet_grid)
              strip.text.y = element_text(size = 10, face = 'bold'), 
              strip.background.y = element_rect(color = "black"),
              strip.text.x = element_text(size = 10)) 


p = p + scale_fill_manual(values = c("#33cccc",  #C>A
                                     "#cc0000",  #C>G
                                     "#00b300",  #C>T
                                     "#ff1a75",  #T>A
                                     "#ffcc00",  #T>C
                                     "#400080")) #T>G

p

ggsave("SomaticSig_ReconstructedSpectrum.svg", dpi = 600, width=20, height=7, units = 'cm')
ggsave("SomaticSig_ReconstructedSpectrum.png", dpi = 600, width=20, height=7, units = 'cm')

rm(p)

p = plotObservedSpectrum(sigs_nmf)


p = p + xlab("")

p = p + theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.grid.major=element_blank(),
              
              axis.ticks = element_line(color = 'black', linewidth = 0.2),
              axis.title.y = element_text(size = 12, color = 'black'),
              axis.text.y = element_text(size = 8, color = 'black'),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 7, color = 'black'),

              # altera layout dos labels (facet_grid) 
              strip.text.y = element_text(size = 9, face = 'bold'),
              strip.background.y = element_rect(color = "black"),
              strip.text.x = element_text(size = 12)) 

p = p + scale_fill_manual(values = c("#F2F44DFF", 
                                     "#F3A002FF"))

p

ggsave("SomaticSig_ObservedSpectrum.svg", dpi = 600, width=30, height=12, units = 'cm')
ggsave("SomaticSig_ObservedSpectrum.png", dpi = 600, width=30, height=12, units = 'cm')


sigs_nmf_bkp = sigs_nmf
sigs_nmf@observed = sigs_nmf@observed[!grepl('CT', rownames(sigs_nmf@observed)),]

rm(p)
p = plotObservedSpectrum(sigs_nmf) 
p = p + theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.grid.major=element_blank(),
              
              axis.ticks = element_line(color = 'black', linewidth = 0.2),
              axis.title.y = element_text(size = 12, color = 'black'),
              axis.text.y = element_text(size = 8, color = 'black'),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 7, color = 'black'),
              
              # altera layout dos labels (facet_grid) 
              strip.text.y = element_text(size = 9, face = 'bold'),
              strip.background.y = element_rect(color = "black"),
              strip.text.x = element_text(size = 12))  

p = p + scale_fill_manual(values = c("#F2F44DFF", 
                                     "#F3A002FF")) 

ggsave("SomaticSig_ObservedSpectrum_SEM_CT.svg", dpi = 600, width=20, height=16, units = 'cm')
ggsave("SomaticSig_ObservedSpectrum_SEM_CT.png", dpi = 600, width=20, height=16, units = 'cm')



sigs_nmf = sigs_nmf_bkp
sigs_nmf@observed = sigs_nmf@observed[grepl('CT', rownames(sigs_nmf@observed)),]


rm(p)
p = plotObservedSpectrum(sigs_nmf) 
p = p + theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.grid.major=element_blank(),
              
              axis.ticks = element_line(color = 'black', linewidth = 0.2),
              axis.title.y = element_text(size = 12, color = 'black'),
              axis.text.y = element_text(size = 8, color = 'black'),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 7, color = 'black'),
              
              # altera layout dos labels (facet_grid) 
              strip.text.y = element_text(size = 9, face = 'bold'),
              strip.background.y = element_rect(color = "black"),
              strip.text.x = element_text(size = 12))   

p = p + scale_fill_manual(values = c("#F2F44DFF", 
                                     "#F3A002FF")) 
p


ggsave("SomaticSig_ObservedSpectrum_SO_COM_CT.svg", dpi = 600, width=5.6, height=16, units = 'cm')
ggsave("SomaticSig_ObservedSpectrum_SO_COM_CT.png", dpi = 600, width=5.6, height=16, units = 'cm')



rm(p)

# S4 class: subset with @
sigs_nmf@samples

write.table(sigs_nmf@samples, "SomaticSig_deNovoSampleMap.txt", quote = FALSE)

sigs_nmf_dt = as.data.table(sigs_nmf@samples,keep.rownames = TRUE)
sigs_nmf_dt = reshape2::melt(sigs_nmf_dt)

setnames(sigs_nmf_dt, old = 'rn', new = 'sample')
setnames(sigs_nmf_dt, old = 'variable', new = 'signature')

sigs_nmf_dt$signature = factor(sigs_nmf_dt$signature)
sigs_nmf_dt$sample = factor(sigs_nmf_dt$sample, levels = c('CO', 'XPV'))

p = plotSampleMap(sigs_nmf)
rm(p)
p = ggplot(sigs_nmf_dt)
p = p + geom_tile(aes(y = sample,
                      x = signature, 
                      fill = value),
                  color = "black")

p = p + scale_fill_gradient(low = "#ffffff",
                            high = "#172869")
p = p + xlab("") + ylab("")

p = p + theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major=element_blank(),
              
              axis.line = element_line(colour = "black", linewidth = 1),
              axis.text = element_text(size=10, colour = "black"),
              axis.text.x = element_text(angle = 0),
              axis.text.y = element_text(margin = margin(0, 6, 0, 0), color = "#425269"),
              axis.ticks = element_blank(),
              
              legend.title = element_blank(),
              legend.text = element_text(size = 7))

p = p + guides(fill = guide_colourbar(barwidth = 0.5,
                                      barheight = 5))

p

ggsave("SomaticSig_SampleMap.svg", dpi = 600, width=7, height=6, units = 'cm')
ggsave("SomaticSig_SampleMap.png", dpi = 600, width=7, height=6, units = 'cm')
