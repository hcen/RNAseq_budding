#install.packages("tidyverse")
library(tidyverse)

#install.packages("readxl")
#library(readxl)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")

#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)


#BiocManager::install("clusterProfiler")
#BiocManager::install("ReactomePA")
library("clusterProfiler")
library("ReactomePA")

#install.packages("pheatmap")
library(pheatmap)
#install.packages("UpSetR")
library(UpSetR)

#BiocManager::install("biomaRt")
library("biomaRt")

#install.packages("devtools") # most recent version of complexheatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) 
#install.packages("gridtext")
library(gridtext)
library(scales)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ##Set working directory to where this file is.
getwd()


#dir.create("input")
#dir.create("output")
#dir.create("figures")

# prepare raw counts and meta data ===================

# read each sample

body1 = read.table("input/Body-1-JZ-082322.counts.genes.coverage.txt") %>% 
  rename(body1=V2)
head(body1)

body2 = read.table("input/Body-2-JZ-081922.counts.genes.coverage.txt") %>% 
  rename(body2=V2)
head(body2)

body3 = read.table("input/Body-3-JZ-082322.counts.genes.coverage.txt") %>% 
  rename(body3=V2)
head(body3)

body4 = read.table("input/Body-4-JZ-082322.counts.genes.coverage.txt") %>% 
  rename(body4=V2)
head(body4)

bud1 = read.table("input/Bud-1-JZ-081922.counts.genes.coverage.txt") %>% 
  rename(bud1=V2)
head(bud1)

bud2 = read.table("input/Bud-2-JZ-081922.counts.genes.coverage.txt") %>% 
  rename(bud2=V2)
head(bud2)

bud3 = read.table("input/Bud-3-JZ-081922.counts.genes.coverage.txt") %>% 
  rename(bud3=V2)
head(bud3)

bud4 = read.table("input/Bud-4-JZ-081922.counts.genes.coverage.txt") %>% 
  rename(bud4=V2)
head(bud1)

# merge all samples into a raw.count data frame

raw.counts <- body1 %>% full_join(body2,by="V1") %>%
  full_join(body3,by="V1") %>%
  full_join(body4,by="V1") %>%
  full_join(bud1,by="V1") %>%
  full_join(bud2,by="V1") %>%
  full_join(bud3,by="V1") %>%
  full_join(bud4,by="V1") %>%
  column_to_rownames(var="V1") # Use gene names as column names
View(raw.counts)

# save raw counts for GEO upload -------
raw.geo <- raw.counts %>%
  mutate(genes=rownames(raw.counts)) %>%
  relocate(genes, .before = body1)
rownames(raw.geo) <- NULL
View(raw.geo)  

write.table(raw.geo,"output/raw_counts.txt", row.names = F, sep = "\t")
raw.geo <- read.table("output/raw_counts.txt", header = T)
raw.geo <- read.delim("output/raw_counts.txt")
# end ---------


raw.counts <- raw.counts[rowSums(raw.counts)!=0,] # Remove genes with all 0
dim(raw.counts) # 22454 genes detected
summary(colSums(raw.counts)) # read depth: min 1443222 median 1654439 mean 1685228 max 1991751 

raw.counts <- round(raw.counts,0) # Counts are not integers. Rounded up because DESeq2 requires integers.
raw.counts <- raw.counts[rowSums(raw.counts)!=0,] # Counts<0.5 became 0. Remove genes with all 0 again
dim(raw.counts) # 19237 genes left


# prepare meta data

sample <- colnames(raw.counts)
group <- c(rep("Body",4),rep("Bud",4))
meta.data <- data.frame(sample, group)
row.names(meta.data) <- meta.data$sample
all(colnames(raw.counts)==rownames(meta.data))
meta.data$group <- factor(meta.data$group, levels=c("Body",
                                                    "Bud"))
View(meta.data)

# DESeq2 and PCA =============

count.data.set <- DESeqDataSetFromMatrix(countData=raw.counts, 
                                         colData=meta.data, design= ~ group) 

# Filter out low count
nrow(count.data.set) # 19237 genes

keep <- rowSums(counts(count.data.set)>=5) >= ncol(meta.data)*0.25  # Keep genes with >5 counts in at least 25% samples. 
count.filter <- count.data.set[keep,]
nrow(count.filter) # 14883 genes left.

# create DESeq object
count.data.set.object <- DESeq(count.filter)

vsd <- vst(count.data.set.object)
View(vsd)
norm.data = assay(vsd)
head(norm.data)
write.table(norm.data, sep="\t",file="output/Norm_data.txt", row.names=TRUE,col.names=NA,quote=FALSE)

# cluster
sampleDists <- dist(t(norm.data),  method = "euclidean") # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters=hclust(sampleDists)
plot(clusters)

# PCA
# modify the PCA function to change format ---------------------------
getMethod("plotPCA","DESeqTransform")
plotPCA.format <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", 
                                fill = "group")) + geom_point(size = 4, shape=21,alpha=0.6) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed() #+  geom_label_repel((aes(label=sample)))
  }
  .local(object, ...)
}


#------------------------------------------------------------------------

# use color blind friendly palette
cbPalette <- c("#E69F00", #lightorange
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#F0E442", #yellow
               "#999999" #grey
)

#install.packages("ggrepel")
library(ggrepel) # https://ggrepel.slowkow.com/articles/examples.html

plotPCA.format(vsd, intgroup=c("group"))+ 
  geom_text_repel(aes(label=colData(vsd)$sample),size=3,
                  color="grey50",
                  box.padding   = 0.4,
                  point.padding = 0,
                  #force=1,
                  #force_pull=10,
                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
                  segment.color = 'grey50')+
  scale_fill_manual(values=cbPalette) +
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank()
  )+
  theme(aspect.ratio=1/1)
# FYI, this is the default function plotPCA()
#plotPCA(vsd, intgroup=c("group"))

ggsave(filename="figures/PCA2.png",width=9,height=7,units="cm",dpi=400)

# sample matrix

library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
View(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         legend = FALSE,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
?pheatmap
pDist<-pheatmap(sampleDistMatrix,
                legend = FALSE,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                cellwidth = 10, cellheight = 10,
                treeheight_row=10,
                treeheight_col=10,
                col = colors)

pDist

png(file = "figures/dist_matrix.png",
    width = 6, 
    height = 6, 
    units = "cm", res = 400)

pDist
# plot upset above
dev.off()

save_pheatmap_png <- function(x, filename, width=800, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
View(pDist)
pDist$gtable
save_pheatmap_png(pDist, "figures/pDist'.png")
dev.off()
#

# plot the VST normalized counts on heatmap =============


View(norm.data)
m<- norm.data

#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)
m.z <- t(scale(t(m))) #%>% as.data.frame()
View(m.z)
colnames(m.z)
ceiling(max(abs(m.z)))
#m.z[m.z>5] <- NA
#m.z <- t(scale(t(m.z)))
#m.z[is.na(m.z)] <- 5

number_of_body <- length(grep("body", colnames(m.z)))
number_of_bud <- length(grep("bud", colnames(m.z)))

end_index_body <- grep("body", colnames(m.z))[number_of_body]
end_index_bud <- grep("bud", colnames(m.z))[number_of_bud]


###

#heatmap_list_total <- 
heatmap.b <- Heatmap(m.z, #matrix_Z_score_total,
                         name = "Z score",
                         show_row_names = FALSE,
                         show_column_names = FALSE,
                         show_row_dend = TRUE,
                         # row_labels = gt_render(m.anno$protein.symbol),
                         row_names_gp = gpar(fontsize = 8),
                         column_names_side = "top",
                         column_dend_side = "bottom",
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "ward.D2",
                     row_dend_side = "left",
                     row_dend_width = unit(5, "mm"),
                         layer_fun = function(j, i, x, y, width, height, fill) {
                           mat = restore_matrix(j, i, x, y)
                           ind = unique(c(mat[, c(end_index_body#, 
                                                  #end_index_bud
                                                  )]))
                           grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                     y = y[ind], 
                                     width = unit(0.03, "inches"), 
                                     height = unit(1/nrow(m.z), "npc"),
                                     gp = gpar(col = "white")
                           )
                         },
                         col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
                         top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                              , height = unit(12, "mm")
                         )),
                         
                         column_order = 1:ncol(m.z),
                         height = 
                           
                           
                           unit(80, "mm"), 
                         
                         width = ncol(m.z)*unit(6, "mm"),
                         border_gp = gpar(col = "black"),
                         show_heatmap_legend = TRUE,
                         heatmap_legend_param = list(
                           title = "Z-score",
                           title_position = "topleft",
                           legend_height = unit(4, "cm")))

draw(heatmap.b)
#

png(file = "figures/heatmap_all.png",
    width = 3, 
    height = 4, 
    units = "in", res = 600)


draw(heatmap.b)
dev.off()
#
seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

####Condition labels

#Condition label 1


grid.rect(x = (loc2$x - loc1$x)*(end_index_body)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_body)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[1], 0.5),
                    col = alpha(cbPalette[1], 0.5)
          )
)
grid.text(expression("Body"), 
          x = (loc2$x - loc1$x)*(end_index_body)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2

grid.rect(x = (loc2$x - loc1$x)*(end_index_body + 
                                   end_index_bud)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_bud)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[2], 0.5),
                    col = alpha(cbPalette[2], 0.5)
          )
)
grid.text(expression("Bud"), 
          x = (loc2$x - loc1$x)*(end_index_body + 
                                   end_index_bud)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))





###Top label gaps

#Vertical lines
grid.rect(x = end_index_body/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.03, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))


dev.off()


# DEG ============

#resultsNames(count.data.set.object)


#res <- results(count.data.set.object, contrast=c("group",exp,ctrl),alpha=0.05)
res <- results(count.data.set.object, contrast=c("group","Bud","Body"),alpha=0.05)


# save result summary
summary(res) 
out <- capture.output(summary(res))
cat("Bud vs Body", out, file="output/results_summary.txt", sep="\n", append=TRUE)

# save results - all genes
res <-  na.omit(res) # omit NA
res  <-  res[order(res$padj),] 
res <- as.data.frame(res)
res$symbol <- rownames(res)
res$entrez.s <- mapIds(org.Hs.eg.db, keys=res$symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first")
res$entrez.a <- mapIds(org.Hs.eg.db, keys=res$symbol, column="ENTREZID", keytype="ALIAS", multiVals="first")
res<- res %>% mutate(entrez=coalesce(entrez.s,entrez.a))

res$ensembl.s <- mapIds(org.Hs.eg.db, keys=res$symbol, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
res$ensembl.a <- mapIds(org.Hs.eg.db, keys=res$symbol, column="ENSEMBL", keytype="ALIAS", multiVals="first")
columns(org.Hs.eg.db)
res <- res %>% mutate(ensembl=coalesce(ensembl.s,ensembl.a)) 

res <- res %>% dplyr::select(-c(entrez.s,entrez.a,ensembl.s,ensembl.a)) 

View(res)
write.table(res, sep="\t",file=paste0("output/Results.txt"), row.names=TRUE,col.names=NA,quote=FALSE)

# save results - significant genes (DEG)
res.sig <- res[res$padj <= 0.05,]
dim(res.sig) # 8553 DE

write.table(res.sig, sep="\t",file=paste0("output/Results_DE.txt"), row.names=TRUE,col.names=NA,quote=FALSE)

res.up <- res.sig[res.sig$log2FoldChange>0,]
res.down <- res.sig[res.sig$log2FoldChange<0,]

dim(res.up) #4325 up
dim(res.down) # 4228 down




# Selected gene list ===============

selected <- c("PDX1", "NKX6-1", "NEUROG3", "NEUROD1", "CHGA", 
              "CHGB", "ISL1", "MNX1", "NKX2-2", "RFX6", "PAX4", 
              'PAX6', 'ARX', 'HHEX', 'FOXA2', 'FEV', 'INS', 'GCG', 
              'GHRL', 'GCK', 'G6PC2', 'SYP', 'SIX2', 'MAFA', 'MAFB', 
              'UCN3', 'ENTPD3', 'SLC30A8', 'KCNJ11', 'ABCC8', 
              'ROBO2', 'CDH1', 'SOX9', 'KRT18', 'KRT19', 'CFTR', 
              'SPP1', 'MUC1', 'TPM1', 'MCM2', 'MCM3', 'ONECUT1', 
              'AMY2B', 'NKX2-5', 'PBX1', 'PBX2', 'PBX3', 'VIM', 
              'CD90 (THY1)', 'COL1A1', 'YAP1', 'NOTCH1', 
              'NOTCH2', 'NOTCH3', 'HES1', 'CTNNB1', 'WNT3A', 
              'PCNA', 'CDK2', 'SLIT1', 'SLIT2', 'SLIT3', 'EPHB1', 
              'EPHB2', 'EPHB3', 'EPHB4')


selected[-which(selected %in% res$symbol)]
# "MAFA"        "CFTR"        "NKX2-5"      "CD90 (THY1)"

raw.counts[c(selected[-which(selected %in% res$symbol)],"THY1"),]
# body1 body2 body3 body4 bud1 bud2 bud3 bud4
# MAFA       1     1     1     1    2    1    2    1
# CFTR       1     1     1     2    1    1    2    0
# NKX2-5     2     1     2     2    1    0    0    0
# NA        NA    NA    NA    NA   NA   NA   NA   NA
# THY1      33    31    38    36    7    6    7    2

selected.new <- c(selected[which(selected %in% res$symbol)],"THY1")

selected.new[-which(selected.new %in% res.sig$symbol)]
# "KRT19" "AMY2B" "PBX2"  "SLIT3" are not DE


# plot the selecte gene (VST normalized counts) on heatmap ------

# plot selected genes only
View(norm.data)
m<- norm.data[selected.new,]

# plot selected genes on the heatmap of all DE genes
m <- norm.data[res.sig$symbol,]
at <- match(selected.new,rownames(m))
labels <- selected.new

#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)
m.z <- t(scale(t(m))) #%>% as.data.frame()
View(m.z)
colnames(m.z)
ceiling(max(abs(m.z)))

number_of_body <- length(grep("body", colnames(m.z)))
number_of_bud <- length(grep("bud", colnames(m.z)))

end_index_body <- grep("body", colnames(m.z))[number_of_body]
end_index_bud <- grep("bud", colnames(m.z))[number_of_bud]

## selected genes
heatmap.b <- Heatmap(m.z, #matrix_Z_score_total,
                     name = "Z score",
                     
                     show_row_names = TRUE,
                     
                     show_column_names = FALSE,
                     show_row_dend = TRUE,
                     row_labels = gt_render(row.names(m)),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_side = "top",
                     column_dend_side = "bottom",
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "ward.D2",
                     row_dend_side = "left",
                     row_dend_width = unit(5, "mm"),
                     layer_fun = function(j, i, x, y, width, height, fill) {
                       mat = restore_matrix(j, i, x, y)
                       ind = unique(c(mat[, c(end_index_body#, 
                                              #end_index_bud
                       )]))
                       grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                 y = y[ind], 
                                 width = unit(0.03, "inches"), 
                                 height = unit(1/nrow(m.z), "npc"),
                                 gp = gpar(col = "white")
                       )
                     },
                     col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
                     top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                          , height = unit(12, "mm")
                     )),
                     
                     column_order = 1:ncol(m.z),
                     height = 
                       
                       unit(180, "mm"), 
                     
                     width = ncol(m.z)*unit(6, "mm"),
                     border_gp = gpar(col = "black"),
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param = list(
                       title = "Z-score",
                       title_position = "topleft",
                       legend_height = unit(4, "cm"))) 
draw(heatmap.b)

### all DE genes + label selected genes
heatmap.b <- Heatmap(m.z, #matrix_Z_score_total,
                     name = "Z score",
                     
                     #show_row_names = TRUE,
                     show_row_names = FALSE,
                     
                     show_column_names = FALSE,
                     show_row_dend = TRUE,
                     row_labels = gt_render(row.names(m)),
                     row_names_gp = gpar(fontsize = 8),
                     column_names_side = "top",
                     column_dend_side = "bottom",
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "ward.D2",
                     row_dend_side = "left",
                     row_dend_width = unit(5, "mm"),
                     layer_fun = function(j, i, x, y, width, height, fill) {
                       mat = restore_matrix(j, i, x, y)
                       ind = unique(c(mat[, c(end_index_body#, 
                                              #end_index_bud
                       )]))
                       grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                 y = y[ind], 
                                 width = unit(0.03, "inches"), 
                                 height = unit(1/nrow(m.z), "npc"),
                                 gp = gpar(col = "white")
                       )
                     },
                     col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
                     top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                          , height = unit(12, "mm")
                     )),
                     
                     column_order = 1:ncol(m.z),
                     height = 
                       
                       unit(180, "mm"), 
                     
                     width = ncol(m.z)*unit(6, "mm"),
                     border_gp = gpar(col = "black"),
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param = list(
                       title = "Z-score",
                       title_position = "topleft",
                       legend_height = unit(4, "cm"))) +
  rowAnnotation(label = anno_mark(at = at, labels = labels,
                                  labels_gp = gpar(col = "black", fontsize = 8))
  )

draw(heatmap.b)


## plot selected genes
png(file = "figures/heatmap_selected.png",
    width = 4, 
    height = 8, 
    units = "in", res = 600)

## plot all DE + selected
png(file = "figures/heatmap_allDE_selected.png",
    width = 4, 
    height = 8, 
    units = "in", res = 600)


draw(heatmap.b)
dev.off()
#
seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

####Condition labels

#Condition label 1


grid.rect(x = (loc2$x - loc1$x)*(end_index_body)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_body)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[1], 0.5),
                    col = alpha(cbPalette[1], 0.5)
          )
)
grid.text(expression("Body"), 
          x = (loc2$x - loc1$x)*(end_index_body)/2/ncol(m.z),
          y = 0.35,
          just = c("center", "top"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2

grid.rect(x = (loc2$x - loc1$x)*(end_index_body + 
                                   end_index_bud)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_bud)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[2], 0.5),
                    col = alpha(cbPalette[2], 0.5)
          )
)
grid.text(expression("Bud"), 
          x = (loc2$x - loc1$x)*(end_index_body + 
                                   end_index_bud)/2/ncol(m.z),
          y = 0.35,
          just = c("center", "top"),
          gp = gpar(fontsize = 11))





###Top label gaps

#Vertical lines
grid.rect(x = end_index_body/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.03, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))


dev.off()


# ORA for correlated genes ==============

geneList <- res$log2FoldChange
names(geneList) <- as.character(res$entrez)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

all_genes <- as.character(res$entrez)

x <- res.sig
x <- res.up
x <- res.down

ora<- function(x){ # Somehow when the code is wrapped inside the function, it can't convert ID to symbol. So run the code inside line-by-line instead...
  
  #KEGG ORA
  KEGG <- enrichKEGG(gene         = x$entrez,
                     organism     = 'hsa',
                     universe      = all_genes,
                     pvalueCutoff=1, pAdjustMethod="BH", 
                     qvalueCutoff=1)
  KEGG <<- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  KEGG.df <<- as.data.frame(KEGG)
  head(KEGG.df)
  #KEGG.df2 <- KEGG.df
  #KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\ ",replacement = "_") 
  #KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\,",replacement = ".")
  #write.xlsx2(KEGG.df, file="output/correlated_genes_ORA.xlsx", sheetName = "KEGG",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  #write.csv(KEGG.df2, sep="\t",file="output/correlated_genes_KEGG.csv", row.names=TRUE,col.names=NA,quote=FALSE)
  
  #Reactome ORA
  react <- enrichPathway(gene         = x$entrez,
                         organism     = 'human',
                         universe      = all_genes,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff=1, pAdjustMethod="BH", 
                         qvalueCutoff=1)
  react <<- setReadable(react, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  react.df <<- as.data.frame(react)
  head(react.df)
  #react.df2 <- react.df
  #react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\ ",replacement = "_") 
  #react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\,",replacement = ".")
  #write.xlsx2(react.df, file="output/correlated_genes_ORA.xlsx", sheetName = "Reactome",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gobp <<- enrichGO(gene        = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gobp.df <<- as.data.frame(gobp)
  head(gobp.df)
  #write.xlsx2(gobp.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_BP",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gomf <<- enrichGO(gene       = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gomf.df <<- as.data.frame(gomf)
  head(gomf.df)
  #write.xlsx2(gomf.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_MF",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gocc <<- enrichGO(gene      = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gocc.df <<- as.data.frame(gocc)
  head(gocc.df)
  #write.xlsx2(gocc.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_CC",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
}

# save ORA all correlated genes
write.csv(KEGG.df, file = "output/ORA_kegg.csv")
write.csv(react.df, file = "output/ORA_reactome.csv")
write.csv(gobp.df, file = "output/ORA_gobp.csv")
write.csv(gomf.df, file = "output/ORA_gomf.csv")
write.csv(gocc.df, file = "output/ORA_gocc.csv")

#save ORA upregulated genes
kegg.up <- KEGG.df

write.csv(KEGG.df, file = "output/ORA_up_kegg.csv")
write.csv(react.df, file = "output/ORA_up_reactome.csv")
write.csv(gobp.df, file = "output/ORA_up_gobp.csv")
write.csv(gomf.df, file = "output/ORA_up_gomf.csv")
write.csv(gocc.df, file = "output/ORA_up_gocc.csv")


#save ORA downregulated genes
kegg.down <- KEGG.df

write.csv(KEGG.df, file = "output/ORA_down_kegg.csv")
write.csv(react.df, file = "output/ORA_down_reactome.csv")
write.csv(gobp.df, file = "output/ORA_down_gobp.csv")
write.csv(gomf.df, file = "output/ORA_down_gomf.csv")
write.csv(gocc.df, file = "output/ORA_down_gocc.csv")


# plot pathway erichment results ============
df.up <-  kegg.up[kegg.up$p.adjust<0.005,] %>% mutate(Log10adj.P=log10(p.adjust)*(-1)) #rank.down
df.up$GeneRatio_num <- sapply(df.up$GeneRatio, function(x) eval(parse(text=x)))
df.up$title <- "Up in Bud"

df.down <-  kegg.down[kegg.down$p.adjust<0.005,] %>% mutate(Log10adj.P=log10(p.adjust)*(-1)) #rank.down
df.down$GeneRatio_num <- sapply(df.down$GeneRatio, function(x) eval(parse(text=x)))
df.down$title <- "Up in Body"

p.up <- ggplot(df.up, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(limits=c(0,5), 
                        #midpoint = 2.1, 
                        low = "white", high = "red", space = "Lab" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank()) +
  facet_grid(. ~ title) +
  theme(strip.background =element_rect(fill=alpha(cbPalette[2], 0.5) ))+
  theme(strip.text = element_text(colour = 'black'))

p.up

p.down <- ggplot(df.down, 
               aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(limits=c(2,20), 
                        #midpoint = 2.1, 
                        low = "white", high = "red", space = "Lab" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x  = element_blank()) +
  facet_grid(. ~ title) +
  theme(strip.background =element_rect(fill=alpha(cbPalette[1],0.5) ))+
  theme(strip.text = element_text(colour = 'black'))
p.down

#install.packages("cowplot")
library(cowplot)

plot_grid(p.up, p.down, ncol = 1, align = "v")

ggsave(filename="figures/pathways.png",width=20,height=25,units="cm",dpi=400)

# end ============
2.48*12+2.48*12+2.48*4+2.48*4+26.98*2
