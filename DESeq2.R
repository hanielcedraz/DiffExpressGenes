suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("edgeR"))


folder = "folder_conteining_genes_counts" 
annof <- "annotation_file_from_biomart.txt"

samples <- read.table("samples_all.txt",header=T, as.is=T)
countf <- file.path(folder, paste0(samples$SAMPLE_ID, '_ReadsPerGene.counts'))

counts = readDGE(countf)$counts
dim(counts)
colnames(counts)
# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features:
noint = rownames(counts) %in% c("N_multimapping","N_noFeature","N_ambiguous")
cpms = cpm(counts)

mean(colSums(counts[!noint,])/colSums(counts))


# In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates, 
keep = rowSums(cpms >1) >=3 & !noint
counts = counts[keep,]
class(counts)


seqy <- DESeqDataSetFromMatrix(counts, samples, design = ~Sexo + Treatment + Sexo*Treatment)

model <- model.matrix(~0 + Sexo + Treatment + Sexo*Treatment + Treatment:Mae, samples)
colnames(model)
all.zero <- apply(model, 2, function(x) all(x==0))
all.zero

idx <- which(all.zero)
m1 <- model[,-idx]
unname(m1)

counts <- DESeq(seqy, full = model)

coef(counts)

res <- results(counts, contrast = )

class(res)