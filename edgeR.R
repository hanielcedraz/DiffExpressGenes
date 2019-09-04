#
# 1. Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count:

Install_Multiples_Packages <- function(packages) {
    pack <- packages[!(packages %in% installed.packages()[,'Package'])];
    if (length(pack)) {
        install.packages(pack, repos = 'https://cran.rstudio.com/')
    }

    for (package_i in packages) {
        suppressPackageStartupMessages(library(package_i, character.only = TRUE, quietly = TRUE))
    }
    write(paste("Packages instaled and load successfuly"), stderr())
}

Install_Multiples_Packages(c('edgeR', 'gplots'))

#suppressPackageStartupMessages(library("DESeq2"))
#suppressPackageStartupMessages(library("edgeR"))


folder = "folder_conteining_genes_counts" #Change it to your folder
annof <- "annotation_file_from_biomart.txt" #change it to your biomart file

samples <- read.table("samples.txt",header=T, as.is=T)
countf <- file.path(folder, paste0(samples$SAMPLE_ID, '_ReadsPerGene.counts'))

counts = readDGE(countf)$counts
dim(counts)

# 2. Filter weakly expressed and noninformative (e.g., non-aligned) features:
noint = rownames(counts) %in% c("N_multimapping","N_noFeature","N_ambiguous")
cpms = cpm(counts)

mean(colSums(counts[!noint,])/colSums(counts))


# In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates, 
keep = rowSums(cpms >1) >=3 & !noint
counts = counts[keep,]
head(counts)
dim(counts)


sample_id = apply(samples[,c("Animal_ID", "Treatment")], 
              1, function(x) paste(na.exclude(x), collapse = ".")) 
condition = apply(samples[,c("Tissue", "Treatment")], 
              1, function(x) paste(na.exclude(x), collapse = ".")) 

# 4. Create a DGEList object (edgeRs container for RNA-seq count data):
d = DGEList(counts=counts, group=condition)

# 5. Estimate normalization factors using:
d = calcNormFactors(d)


# 6. Inspect the relationships between samples using a multidimensional scaling (MDS) plot, as shown in Figure 4:
pdf("MDS-any_information.pdf")
plotMDS(d, labels=sample_id,
        col = rainbow(length(levels(factor(condition))))[factor(condition)], cex = 0.6, main = "MDS")
dev.off()



# 7. Estimate tagwise dispersion (simple design) using:
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

# 8. plot the mean-variance relationship:
pdf("mean.variance-any_information.pdf")
plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE, 
            main = "MeanVar")
plotBCV(d, main = "BCV")
dev.off()


# edgeR Completely Randomized Design
# 9. Create a design matrix to specify the factors that are expected to affect expression levels:
design = model.matrix( ~ -1 + condition)
colnames(design) =c("Level_1","Level_2")
rownames(design) = samples$Animal_ID
design

# 10. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
d2 = estimateGLMTrendedDisp(d, design)
d2 = estimateGLMTagwiseDisp(d2, design)


# 11. Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)
f

# 12. Perform a likelihood ratio test, specifying the difference of interest
### Level_1 - Level_2
de = glmLRT(f, contrast = c(1, -1))


degeFC = de$table
head(degeFC)



# 13. Use the topTags function to present a tabular summary of the differential expression statistics
tt = topTags(de, n = nrow(d))
head(tt$table)
table(tt$table$FDR < 0.05)


# 14. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:

nc = cpm(d, normalized.lib.sizes = TRUE)
rn = rownames(tt$table)
head(nc[rn, order(samples$Treatment)], 5)

write.table(data.frame((nc[rn, order(samples$Treatment)])), 
            sep = "\t", file = "cpms_any_information.txt")


# 15. Plot the M (log-fold change) versus A (log-average expression)
deg = rn[tt$table$FDR < .05]
pdf("smear-any_information.pdf")
plotSmear(de, de.tags = deg, main = "Smear")
abline(h = c(-2, 0, 2), col = c("dodgerblue", "yellow", "dodgerblue"), lty = 2)
dev.off()


# 16. Save the result table as a CSV file:
### attach annotation

anno <- read.table(annof, sep = "\t", header = T, 
                   comment.char = "", quote = "", as.is = T)
write.csv(data.frame(tt$table, anno[match(rownames(tt$table), 
                    anno$Gene.stable.ID),]),
              file="toptags_any_information.csv")


#--------------------------------HEATMAP--------------------------


logcounts = cpm(d, log=TRUE)
head(logcounts)

var_genes = apply(logcounts, 1, var)
head(var_genes)


# Selecting the most 500 variable genes:
highly_variable_lcpm = logcounts[names(sort(var_genes, decreasing=TRUE))[1:500],]
dim(highly_variable_lcpm)



# Plotting heatmap:
png("heatmap_.png",    # create PNG for the heat map        
    width = 5*600,        # 5 x 300 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(highly_variable_lcpm,
          #main = "DEG 21 days", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=redgreen(75), # use on color palette defined earlier
          scale = "row",
          #breaks=col_breaks,    # enable color transition at specified limits
          Rowv =TRUE,
          labCol = samples$SAMPLE_ID,
          Colv=,
          distfun = dist,
          hclustfun = hclust,
          dendrogram = "both",
)
dev.off()  


