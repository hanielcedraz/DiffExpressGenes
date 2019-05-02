suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))


#setwd('/Volumes/HANIEL_RED/Susana/')
folder = "folder_containing_gene_counts_reads" 
annof <- "annnoatation_from_biomart.txt"

samples <- read.table("samples_all.txt",header=T, as.is=T)
countf <- file.path(folder, paste0(samples$SAMPLE_ID, '_ReadsPerGene.counts'))


counts = readDGE(countf)$counts
class(counts)
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

dge <- DGEList(counts=counts)

samples <- data.frame(SAMPLE_ID = samples$SAMPLE_ID, Mae = as.factor(samples$Mae), Tissue = as.factor(samples$Tissue), Sexo = as.factor(samples$Sexo), Treatment = as.factor(samples$Treatment))


design <- model.matrix(~0+ model) #Set the used model
colnames(design) <- c() #Set the names of col design




keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=3)



dim (design)
dim(counts)
dupC <- duplicateCorrelation(object = logCPM, design = design, block = samples$Mae)
dupC$consensus.correlation

#fit <- lmFit(counts, design, correlation = dupC$consensus)



fit <- lmFit(logCPM, design)

#make contrasts
cm <- makeContrasts(
    
    levels = design
)

nonEstimable(design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
names(fit2)



tab <- topTable(fit2, coef = 'contrast')
