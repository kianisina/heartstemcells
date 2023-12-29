#BiocManager::install("DESeq2")
#read data
countdata <- read.table("H_paired_end_readcounts.txt", header=TRUE, row.names=1)
#adjust data
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("X.prj.heartstemcells.SK22.Mapping.H.sortedByCoordinate.", "", colnames(countdata))
countdata
countdata <- as.matrix(countdata)
head(countdata)

strain <- factor(c(rep("H",15)))
condition <- factor(c(rep(c("control", "Ser", "Ser_p38", "Ser_TGFb", "NFkB"),3)))

coldata <- data.frame(row.names = colnames(countdata), strain, condition)
coldata$group <- factor(paste0(coldata$strain,"_",coldata$condition))
coldata

library(DESeq2)
#DESeq2
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group)
dds <- DESeq(dds)


coldata
# save degs
get_degs <- function(condition_1, condition_2){
  if(condition_1 != condition_2) {
    res <- results(dds, contrast = c("group", condition_1, condition_2), alpha = 0.05)
    res <- res[order(res$padj),]
  
    condition_names <- rownames(coldata[(coldata$group == condition_1) | (coldata$group == condition_2),])
    res_counts <- merge(as.data.frame(res),as.data.frame(counts(dds[,condition_names],normalized=T)), by="row.names",sort=F)
    names(res_counts)[1] <- "Gene"
    res_counts <- subset(res_counts, padj <= 0.05)
    write.csv(res_counts, file = paste0("DEGs...",condition_1," vs ",condition_2, "sig.csv"),row.names = F)
    
    #write.csv(res_counts, file = paste0("DEGs...",condition_1," vs ",condition_2, ".csv"),row.names = F)
  }
}
get_degs("H_Ser", "H_control")
for (co in c("Ser", "Ser_p38", "Ser_TGFb", "Ser_NFkB")) {
  get_degs("H_control", paste0("H_",co))
  
}
for (co in c("control", "Ser_p38", "Ser_TGFb", "Ser_NFkB")) {
  get_degs(paste0("H_",co),"H_Ser" )
  
}


#plot

png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
# heatmap with all conditions
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()




#Vulcano
getVulcano <- function(nameOfDEGs) {
  data_table <- read.delim("DEGs...H_Ser_NFkB vs H_Ser(p<0.05).csv", sep = ",",header = T)
  nameOfDEGs <- substring(nameOfDEGs, 8, nchar(nameOfDEGs)-4)
  
  head(data_table)
  # remove all NA's
  data_table[data_table ==0] <- NA
  data_table <- data_table[,colSums(is.na(data_table))<nrow(data_table)]
  data_table <- na.omit(data_table)
  #we create plot
  library(ggplot2)
  library(ggrepel)
  
  
  #names( rowData(data_table) )
  p <- ggplot(data= data_table, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color = cut(-log10(padj),c(-Inf,-log10(0.05),Inf))),size = 1, show.legend = F) +
    scale_color_manual(values = c("grey", "red") )
    #geom_text_repel(data = head(data_table, 90), aes(label=Gene), size = 3, max.overlaps = 14))
  #ggsave(paste0("vulcano_",nameOfDEGs,".png"))
  p
}
BiocManager::install('EnhancedVolcano')
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(data_table,
                lab = data_table$Gene,
                x = 'log2FoldChange',
                pCutoff = 0.005,
                
                y = 'padj')

ggplot(data= data_table, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = cut(-log10(padj),c(-Inf,-log10(0.05),Inf))),size = 1, show.legend = F) +
  scale_color_manual(values = c("grey", "red","blue") )+
  #geom_hline(yintercept = 1.3, col = 'red')+
  geom_text_repel(data = head(data_table, 90), aes(label=Gene), size = 3, max.overlaps = 14)
  
library(ggplot2)
library(ggrepel)  # Make sure ggrepel is installed

ggplot(data = data_table, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(log2FoldChange > 0, "red", "blue")), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("red", "blue")) +
  geom_hline(yintercept = 1.3, col = 'red') 
  #geom_text_repel(data = head(data_table, 90), aes(label = Gene), size = 3, max.overlaps = Inf)
library(ggplot2)
library(ggrepel)
if(data_table$log2FoldChange<0 & data_table$padj < 0.05) {
  print(data_table$Gene)
}
selected_geness <- data_table$Gene[data_table$log2FoldChange < -4 & data_table$padj < 0.05]
selected_genessup <- data_table$Gene[data_table$log2FoldChange > 0 & data_table$padj < 0.05]
selected_geness
selected_genessup
length(selected_geness)
length(selected_genessup)
ifelse(data_table$log2FoldChange<0 & data_table$padj < 0.05,data_table$Gene,"fa")
ggplot(data = data_table, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(log2FoldChange<0 & padj < 0.05 , "Downregulated", ifelse(log2FoldChange>0 & padj<0.05,"Upregulated","Not Significant"))), size = 1, show.legend = TRUE) +
  scale_color_manual(name = "Regulation",values = c("blue","grey","red")) +

  #geom_hline(yintercept = -log10(0.05), col = 'grey') 
  geom_text_repel(data = head(data_table, 90), aes(label = Gene), size = 3, max.overlaps = 14)

head(data_table)

getVulcano("DEGs...H_Ser_TGFb vs H_Ser(p<0.05).csv")

for (co in c("control", "Ser_p38", "Ser_TGFb", "Ser_NFkB")) {
  getVulcano(paste0("DEGs...H_",co," vs H_Ser(p<0.05).csv"))
  
}
for (co in c("Ser", "Ser_p38", "Ser_TGFb", "Ser_NFkB")) {
  getVulcano(paste0("DEGs...H_control vs H_", co, ".csv"))
  
}





#pca


#read the csv data
data_table <- rld
head(data_table)
assay(data_table)

# remove all NA's
data_table[data_table ==0] <- NA
data_table <- data_table[,colSums(is.na(data_table))<nrow(data_table)]
data_table <- na.omit(data_table)
data_table

#remove target_id column
data_table <-  data_table[-c(1)]
head(data_table)

#transpose the data
data_table <- t(assay(data_table))
data_table <- t(data_table)

#perform principal component analysis
pca <- prcomp(data_table)
prcomp()

s <- summary(pca)
s
pca$x

scores <- as.data.frame(pca$x)
head(scores)
#create plot
library(ggplot2)

pca12 <- ggplot(data = scores, aes(x=PC1,y=PC2)) +
  geom_point() +
  labs(x= paste0("PC1 == ", s$importance[2,1]*100,"%"),
       y= paste0("PC2 == ", s$importance[2,2]*100,"%"))

pca12
row.names(data_table)
lst <- c(1:100)
seq_len(100)
assay(data_table)
sortTheData= order(rowVars(assay(data_table)), decreasing = TRUE)
sortTheData
temp <- assay(data_table)[sortTheData, ]
row.names(temp)
temp <- temp[c('WEE1', 'CUTA', 'AAAS', 'RFX2', 'COX5B'),]
#WEE1, CUTA, AAAS, RFX2 und COX5B
length(sortTheData)
pca = prcomp(t(temp))
pca
rowData(data_table)
names(rowData(data_table))
names(colData(data_table))

data_table
colnames(data_table)
#we color the points 
library(ggrepel)
scores
pca12 <- ggplot(data = scores, aes(x=PC1,y=PC2, color = condition )) +
  geom_point() +
  geom_text_repel(aes(label = rownames(scores))) +
  labs(x= paste0("PC1 == ", s$importance[2,1]*100,"%"),
     y= paste0("PC2 == ", s$importance[2,2]*100,"%"))

pca12
ggsave("pcaaa.png")






#GO-term
packageurl <- "https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
remove.packages("clusterProfiler")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")


nBiocManager::install("clusterProfiler")
BiocManager::install("clusterProfiler")
BiocManager::install(version = '3.9')
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(pathview)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)






#Go term and KEGGs





getPathview <- function(cond1) {
  res <- results(dds, contrast = c("group", "H_NFkB", "H_Ser"))
  res <- res[order(res$padj),]
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  res$entrez = mapIds(org.Hs.eg.db,
                      keys=row.names(res), 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  res$geneSymbol =   mapIds(org.Hs.eg.db,
                            keys=row.names(res), 
                            column="GENENAME",
                            keytype="SYMBOL",
                            multiVals="first")
  res$ensGene =   mapIds(org.Hs.eg.db,
                         keys=row.names(res), 
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
  library(gage)
  library(gageData)
  
  data(kegg.sets.hs)
  data(sigmet.idx.hs)
  
  kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
  foldchanges = res$log2FoldChange
  names(foldchanges) = res$entrez
  
  keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
  attributes(keggres)
  head(keggres$greater)
  head(keggres$less)
  lapply(keggres, head)
  #UP
  upRegulatedPathways = data.frame(id=rownames(keggres$greater), keggres$greater)
  #write.csv(upRegulatedPathways, file = paste0("upRegulatedPathways...", cond1, " vs H_control", ".csv"))
  upRegulatedPathways$stat.mean 
  
  rownames(upRegulatedPathways) <- gsub("\\d", "", rownames(upRegulatedPathways))
  rownames(upRegulatedPathways) <- gsub("hsa ", "", rownames(upRegulatedPathways))
  upRegulatedPathways <- subset(upRegulatedPathways, p.val <= 0.05)
  
  ggplot(data = upRegulatedPathways, aes(x= reorder(rownames(upRegulatedPathways), stat.mean) , y = stat.mean) )+ geom_bar(stat = "identity",width = 0.9 ,fill="darkblue")+
    labs(y = "fold enrichment", x = "KEGG Terms") + theme(axis.text.y=element_text(size=7), axis.text.x=element_text(size=13), axis.title = element_text(size=13)) + coord_flip()   
  
  pax6_up_plot <- ggplot(upRegulatedPathways) +
    geom_col(aes(x = reorder(rownames(upRegulatedPathways), stat.mean), 
                 y = stat.mean), 
             fill = "lightblue", 
             color = "black",
             width=0.6, 
             position = position_dodge(width=0.4)) + 
    coord_flip() +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "fold enrichment", 
                       limits = c(0, 3),
                       expand = c(0, 0)) +
    
    ggtitle("KEGG Terms") +
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'black', 
                                   size = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 10,
                                      face = "italic"),
          axis.text.x = element_text(size = 6),
          plot.title = element_text(size = 9,
                                    vjust = -3)) +
    geom_text(aes(x = reorder(rownames(upRegulatedPathways), stat.mean),
                  y = stat.mean /2,
                  label = reorder(rownames(upRegulatedPathways), stat.mean)),
              size = 8 *.36)+
    geom_text(aes(x = reorder(rownames(upRegulatedPathways), stat.mean),
                  y = stat.mean,
                  label = paste0("", round(p.val, digits = 6))),
              size = 3.5,
              hjust = -0.1,
              vjust = 0.5)
  
  ggsave(paste0("upRegulatedPathways...", cond1, " - H_Ser.png"))
  keggrespathways <- rownames(keggres$greater)[1:5]
  # Extract the IDs part of each string
  keggresids = substr(keggrespathways, start=1, stop=8)
  pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
  #Down
  
  downRegulatedPathways = data.frame(id=rownames(keggres$less), keggres$less)
 
  rownames(downRegulatedPathways) <- gsub("\\d", "", rownames(downRegulatedPathways))
  rownames(downRegulatedPathways) <- gsub("hsa ", "", rownames(downRegulatedPathways))
  downRegulatedPathways <- subset(downRegulatedPathways, p.val <= 0.05)
  #ggplot(data = downRegulatedPathways, aes(x= reorder(rownames(downRegulatedPathways), stat.mean) , y = stat.mean) )+ geom_bar(stat = "identity",width = 0.9 ,fill="darkblue")+
    #labs(y = "fold enrichment", x = "KEGG Terms") + theme(axis.text.y=element_text(size=7), axis.text.x=element_text(size=13), axis.title = element_text(size=13)) + coord_flip()   
  head(downRegulatedPathways,21)
  downRegulatedPathways <- head(downRegulatedPathways, 21)
  ggplot(downRegulatedPathways) +
    geom_col(aes(x = reorder(rownames(downRegulatedPathways), stat.mean), 
                 y = stat.mean), 
             fill = "lightblue", 
             color = "black",
             width=0.8, 
             position = position_dodge(width=0.4)) + 
    coord_flip() +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "fold enrichment", 
                       limits = c(-4, 0),
                       expand = c(0, 0)) +
    
    ggtitle("KEGG Terms") +
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'black', 
                                   size = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 10,
                                      face = "italic"),
          axis.text.x = element_text(size = 13),
          plot.title = element_text(size = 9,
                                    vjust = -3)) +
    geom_text(aes(x = reorder(rownames(downRegulatedPathways), stat.mean),
                  y = stat.mean / 2,
                  label = reorder(rownames(downRegulatedPathways), stat.mean)),
              size = 8 * 0.36) +
    geom_text(aes(x = reorder(rownames(downRegulatedPathways), stat.mean),
                  y = stat.mean,
                  label = sprintf("%.6f", p.val)),
              size = 3.5,
              hjust = 1.1,
              vjust = 0.5,
              y = -3.5)
              
  pax6_up_plot
  
  ggsave(paste0("downRegulatedPathways... - H_Ser.png"),height = 500)
  #write.csv(downRegulatedPathways, file = paste0("downRegulatedPathways...", cond1, " vs H_control", ".csv"))
  
  #keggrespathways <- rownames(keggres$less)[1:5]
  # Extract the IDs part of each string
  #keggresids = substr(keggrespathways, start=1, stop=8)
  pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
  
}

getPathview("H_Ser")

library(pathview)
library(gage)
library(gageData)
getGos <- function(cond1,term) {
  res <- results(dds, contrast = c("group", "H_NFkB","H_Ser" ))
  res <- res[order(res$padj),]
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  res$entrez = mapIds(org.Hs.eg.db,
                      keys=row.names(res), 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  res$geneSymbol =   mapIds(org.Hs.eg.db,
                            keys=row.names(res), 
                            column="GENENAME",
                            keytype="SYMBOL",
                            multiVals="first")
  res$ensGene =   mapIds(org.Hs.eg.db,
                         keys=row.names(res), 
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
  foldchanges = res$log2FoldChange
  names(foldchanges) = res$entrez
  data(go.sets.hs)
  data(go.subs.hs)
  
  if(term == "BP") {
    gobpsets = go.sets.hs[go.subs.hs$BP]
    gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
    lapply(gobpres, head)
    #greaterGo
    greaterGO <- as.data.frame(gobpres$less)
    greaterGO <- subset(greaterGO,p.val <= 0.05)
    rownames(greaterGO) <- gsub("GO:", "", rownames(greaterGO))
    rownames(greaterGO) <- gsub("^\\d+\\s+", "", rownames(greaterGO))
    greaterGO <- as.data.frame(greaterGO)
    greaterGO <- greaterGO[order(greaterGO$p.val),]
    num_rows <- nrow(greaterGO)
    greaterGO <- head(greaterGO, 20)
    greaterGO <- tail(greaterGO, 20)
    #greaterGO <- greaterGO - head(greaterGO, 9)
    #ggplot(data = greaterGO, aes(x= reorder(rownames(greaterGO), stat.mean) , y = stat.mean) )+ geom_bar(stat = "identity",width = 0.9 ,fill="darkblue")+
      #labs(y = "fold enrichment", x = "GO Terms") + theme(axis.text.y=element_text(size=5), axis.text.x=element_text(size=10), axis.title = element_text(size=13)) + coord_flip() + ylim(0, 3)
    
    ggplot(greaterGO) +
      geom_col(aes(x = reorder(rownames(greaterGO), stat.mean), 
                   y = stat.mean), 
               fill = "lightblue", 
               color = "black",
               width=0.8, 
               position = position_dodge(width=0.4)) + 
      coord_flip() +
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "fold enrichment", 
                         limits = c(-3, 0),
                         expand = c(0, 0)) +
      
      ggtitle("GO Terms") +
      theme(panel.border = element_blank(),
            axis.line = element_line(color = 'black', 
                                     size = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10,
                                        face = "italic"),
            axis.text.x = element_text(size = 13),
            plot.title = element_text(size = 9,
                                      vjust = -3)) +
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                    y = stat.mean /2,
                    label = reorder(rownames(greaterGO), stat.mean)),
                size = 8 *.36,
                hjust = 0.5
                )+
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                    y = stat.mean,
                    label = paste0("", round(p.val, digits = 6))),
                size = 3.5,
                hjust = -0.1,
                vjust = 0.5,
                y=-2.8)
    
    
    ggsave(paste0("greaterGOs_BP...", cond1, " - H_Ser-Scale0to-3.3(updatedVersion1).png"))
    
    #write.csv(greaterGO, file = paste0("greaterGOs_BP...", cond1, " vs H_control", ".csv"))  
  }else if(term == "CC") {
    gobpsets = go.sets.hs[go.subs.hs$CC]
    gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
    lapply(gobpres, head)
    #greaterGo
    greaterGO <- as.data.frame(gobpres$less)
    greaterGO <- subset(greaterGO,p.val <= 0.05)
    rownames(greaterGO) <- gsub("GO:", "", rownames(greaterGO))
    rownames(greaterGO) <- gsub("\\d", "", rownames(greaterGO))
    greaterGO <- as.data.frame(greaterGO)
    greaterGO <- greaterGO[order(greaterGO$p.val),]
    #ggplot(data = greaterGO, aes(x= reorder(rownames(greaterGO), stat.mean) , y = stat.mean) )+ geom_bar(stat = "identity",width = 0.9 ,fill="darkblue")+
      #labs(y = "fold enrichment", x = "GO Terms") + theme(axis.text.y=element_text(size=5), axis.text.x=element_text(size=10), axis.title = element_text(size=13)) + coord_flip() + ylim(0, 3)
    
    pax6_up_plot <- ggplot(greaterGO) +
      geom_col(aes(x = reorder(rownames(greaterGO), stat.mean), 
                   y = stat.mean), 
               fill = "lightblue", 
               color = "black",
               width=0.8, 
               position = position_dodge(width=0.4)) + 
      coord_flip() +
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "fold enrichment", 
                         limits = c(-3, 0),
                         expand = c(0, 0)) +
      
      ggtitle("GO Terms") +
      theme(panel.border = element_blank(),
            axis.line = element_line(color = 'black', 
                                     size = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10,
                                        face = "italic"),
            axis.text.x = element_text(size = 13),
            plot.title = element_text(size = 9,
                                      vjust = -3)) +
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                    y = stat.mean /2,
                    label = reorder(rownames(greaterGO), stat.mean)),
                size = 8 *.36,
                hjust= 0.5)+
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                    y = stat.mean,
                    label = paste0("", round(p.val, digits = 6))),
                size = 3.5,
                hjust = -0.1,
                vjust = 0.5,
                y=-2.8)
    
    
    ggsave(paste0("greaterGOs_CC...", cond1, " - H_Ser-Scale0to-3.3(updatedVersion).png"))
    #write.csv(greaterGO, file = paste0("greaterGOs_CC...", cond1, " vs H_control", ".csv"))
  }else if(term == "MF") {
    gobpsets = go.sets.hs[go.subs.hs$MF]
    gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
    lapply(gobpres, head)
    #greaterGo
    greaterGO <- as.data.frame(gobpres$less)
    greaterGO <- subset(greaterGO,p.val <= 0.05)
    greaterGO <- tail(greaterGO, 21)
    
    rownames(greaterGO) <- gsub("GO:", "", rownames(greaterGO))
    rownames(greaterGO) <- gsub("\\d", "", rownames(greaterGO))
    greaterGO <- as.data.frame(greaterGO)
    greaterGO <- greaterGO[order(greaterGO$p.val),]
    
    #ggplot(data = greaterGO, aes(x= reorder(rownames(greaterGO), stat.mean) , y = stat.mean) )+ geom_bar(stat = "identity",width = 0.9 ,fill="lightblue")+
      #labs(y = "fold enrichment", x = "GO Terms") + theme(axis.text.y=element_text(size=5), axis.text.x=element_text(size=10), axis.title = element_text(size=13)) + coord_flip() + ylim(0, 3)+
      #geom_text(aes(label = reorder(rownames(greaterGO), stat.mean)),size= 9/ .pt,position = position_dodge(width = 0.9), hjust = 1,vjust = 0.5, colour = "black")
    
    BiocManager::install("devEMF")
    library(devEMF)
    require(devEMF)
    emf(file = "trial.emf")
    emf(file = "trial.emf", emfPlus = FALSE)
    
     ggplot(greaterGO) +
      geom_col(aes(x = reorder(rownames(greaterGO), stat.mean), 
                   y = stat.mean), 
               fill = "lightblue", 
               color = "black",
               width=0.8, 
               position = position_dodge(width=0.4)) + 
      coord_flip() +
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "fold enrichment", 
                         limits = c(-3, 0),
                         expand = c(0, 0)) +
      
      ggtitle("GO Terms") +
      theme(panel.border = element_blank(),
            axis.line = element_line(color = 'black', 
                                     size = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10,
                                        face = "italic"),
            axis.text.x = element_text(size = 13),
            plot.title = element_text(size = 9,
                                      vjust = -3)) +
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                    y = stat.mean /2,
                    label = reorder(rownames(greaterGO), stat.mean)),
                size = 8 *.36,
                hjust = 0.5)+
      geom_text(aes(x = reorder(rownames(greaterGO), stat.mean),
                  y = stat.mean,
                  label = paste0("", round(p.val, digits = 6))),
              size = 3.5,
              hjust = -0.1,
              vjust = 0.5,
              y=-2.8)
    
    
    dev.off()
    
    ggsave(paste0("greaterGOs_MF...", cond1, " - H_Ser(UpdatedVersion).png"))
    #write.csv(greaterGO, file = paste0("greaterGOs_MF...", cond1, " vs H_control", ".csv"))
  }
  
}
library(ggplot2)
getGos("H_control", "MF")
for (co in c("H_control", "H_Ser_p38", "H_Ser_TGFb", "H_NFkB")) {
  getGos(co, "BP")
  
}




#HEATMAPS

rld <- rlogTransformation(dds)
head(assay(rld))
topVarGenes <- head(order(-rowVars(assay(rld))),50)

hist(assay(rld))
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[rld$condition])
colors <- colorRampPalette( rev(brewer.pal(10, "RdBu")) )(255)
color <- colorRampPalette(rev(brewer.pal(n=7,name= "RdYlBu")))(100)
# Sample distance heatmap

sampleDists <- assay(rld)[ topVarGenes, ]
sampleDists <- sampleDists - rowMeans(sampleDists)
sampleDists <- sampleDists[order(rownames(sampleDists) %in% rownames(sampleDists)[topVarGenes]), ]

library(gplots)
install.packages("pheatmap")
library(pheatmap)
summary(sampleDists)
str(sampleDists)
my_hclust_gene <- hclust(dist(sampleDists), method = "complete")
library(dendextend)
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
png("heatmapSerControlTgfb top50.png", w=1000, h=1000, pointsize=20)

pheatmap(sampleDists, 
         #cluster_rows = FALSE,  # Do not cluster rows
         cluster_cols = FALSE,  # Do not cluster columns
         color = color,        # Set heatmap colors
         fontsize_row = 8,       # Set font size for row labels
         fontsize_col = 8,
         angle_col = 45,        # Rotate column labels by 45 degrees
         #cellwidth = 15,       # Set the width of each cell
         cellheight = 8, # Set the height of each cell
         #annotation_row = my_gene_col,
         #cutree_rows = 3,
         #cutree_cols = 1,
         
         margin = c(8, 8))
dev.off()

png("heatmapSerControlTgfb top50.png", w=1000, h=1000, pointsize=20)
heatmap.2(sampleDists, trace="none",
          col=colors,
          ColSideColors=mycols,
          margin=c(10,10),
          scale = "row")
dev.off()

getHeatmap <- function(cond) {
  strain <- factor(c(rep("H",6)))
  condition <- factor(c(rep(c(cond, "Ser"),3)))
  
  coldata <- data.frame(row.names = colnames(countdata), strain, condition)
  coldata$group <- factor(paste0(coldata$strain,"_",coldata$condition))
  coldata
  
  library(DESeq2)
  #DESeq2
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group)
  dds <- DESeq(dds)
  rld <- rlogTransformation(dds)

  topVarGenes <- head(order(-rowVars(assay(rld))),200)
  
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[rld$condition])
  colors <- colorRampPalette( rev(brewer.pal(10, "RdBu")) )(255)
  # Sample distance heatmap
  
  sampleDists <- assay(rld)[ topVarGenes, ]
  sampleDists <- sampleDists - rowMeans(sampleDists)
  library(gplots)
  png(paste("heatmap...H_", cond," vs H_Ser top50.png"), w=1000, h=1000, pointsize=20)
  heatmap.2(sampleDists, trace="none",
            col=colors,
            ColSideColors=mycols,
            margin=c(10, 10),
            scale = "row")
  dev.off()
  
}
getHeatmap("Ser_NFkB")
strain <- factor(c(rep("H",9)))
condition <- factor(c(rep("control", 3), rep("Ser", 3), rep("NFkB", 3)))
condition <- factor(c(rep(c("control", "Ser", "Ser_NFkB"),1)))
co <- subset(countdata, select = c(H3_control, H5_control, H12_control))
se <- subset(countdata, select = c(H3_Ser, H5_Ser, H12_Ser))
tg <- subset(countdata, select = c(H3_Ser_TGFb, H5_Ser_TGFb, H12_Ser_TGFb))
nf <- subset(countdata, select = c(H3_Ser_NFkB, H5_Ser_NFkB, H12_Ser_NFkB))
control <- rowMeans(subset(countdata, select = c(H3_control, H5_control, H12_control)))
ser <- rowMeans(subset(countdata, select = c(H3_Ser, H5_Ser, H12_Ser)))
p38 <- rowMeans(subset(countdata, select = c(H3_Ser_p38, H5_Ser_p38, H12_Ser_p38)))
tgfb <- rowMeans(subset(countdata, select = c(H3_Ser_TGFb, H5_Ser_TGFb, H12_Ser_TGFb)))
nfkb <- rowMeans(subset(countdata, select = c(H3_Ser_NFkB, H5_Ser_NFkB, H12_Ser_NFkB)))
control <- round(control)
tgfb <- round(tgfb)
ser <- round(ser)
as.matrix(control,ser,tgfb)
summary(ser)
control <- lapply(control,as.numeric)
ser <- lapply(ser,as.numeric)
tgfb <- lapply(tgfb,as.numeric)
data.frame(control,ser,tgfb)
head(countdata)
heatmapWithThreeCondition <- data.frame(co, se, nf)
conditionsInFrame <- data.frame(control, ser, tgfb)
coldata <- data.frame(row.names = colnames(c("H_control", "H_Ser", "H_Ser_NFkB")), strain, condition)
coldata$group <- factor(paste0(coldata$strain,"_",coldata$condition))
coldata
rownames(conditionsInFrame)
colnames(conditionsInFrame) <- c("H_control", "H_Ser", "H_Ser_TGFb")
conditionsInFrame <- lapply(conditionsInFrame,as.numeric)
head(conditionsInFrame)
head(data.frame(conditionsInFrame))
conditionsInFrame <- data.frame(conditionsInFrame)
conditionsInFrame$tgfb <- lapply(conditionsInFrame$tgfb,as.numeric)
conditionsInFrame$nfkb <- lapply(conditionsInFrame$nfkb,as.numeric)
conditionsInFrame$p38 <- lapply(conditionsInFrame$p38,as.numeric)
conditionsInFrame$ser <- lapply(conditionsInFrame$ser,as.numeric)
conditionsInFrame$control <- lapply(conditionsInFrame$control,as.numeric)
round(conditionsInFrame)
feri1 <- as.matrix(conditionsInFrame)
feri2 <- data.matrix(feri1)

# Convert to matrix
library(DESeq2)
#DESeq2
ncol(conditionsInFrame) == nrow(coldata)
str(conditionsInFrame)
sum(!is.finite(conditionsInFrame))
sum(is.na(conditionsInFrame))
# Ensure 'feri' is a data frame
conditionsInFrame <- as.data.frame(conditionsInFrame)

# Check the data types and dimensions again
str(feri1)
conditionsInFrame[] <- lapply(conditionsInFrame, as.numeric)

countdata
dds <- DESeqDataSetFromMatrix(countData=heatmapWithThreeCondition , colData=coldata, design = ~group)
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
results(dds)
#vignette("DESeq2")
condition <- factor(c(rep(c("Ser", "control", "TGFb"),1)))
topVarGenes <- head(order(-rowVars(assay(rld))),50)

library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[condition])
colors <- colorRampPalette( rev(brewer.pal(10, "RdBu")) )(255)
# Sample distance heatmap
head(sampleDists)
sampleDists <- assay(rld)[ topVarGenes, ]
sampleDists <- sampleDists - rowMeans(sampleDists)
head(rowData(rld))
row.names(rowData(rld))
colnames(rowData(rld))
colnames(rld)
colVars(rld)
#write.csv(assay(rld), file = paste0("sampleassay...h", ".csv"))
H_control <- rowMeans(subset(sampleDists, select = c(H3_control, H5_control, H12_control)))
H_Ser <- rowMeans(subset(sampleDists, select = c(H3_Ser, H5_Ser, H12_Ser)))
H_Ser_p38 <- rowMeans(subset(sampleDists, select = c(H3_Ser_p38, H5_Ser_p38, H12_Ser_p38)))
H_Ser_TGFb <- rowMeans(subset(sampleDists, select = c(H3_Ser_TGFb, H5_Ser_TGFb, H12_Ser_TGFb)))
H_Ser_NFkB <- rowMeans(subset(sampleDists, select = c(H3_Ser_NFkB, H5_Ser_NFkB, H12_Ser_NFkB)))
sampleDists <- cbind(H_control,H_Ser,H_Ser_NFkB)
write.csv(sampleDists, file = paste0("ListOfGenesInHeatmap....csv"))
library(gplots)
png(paste("heatmap...ser-control top200c.png"), w=1000, h=800, pointsize=20)
heatmap.2(sampleDists, trace="none",
          col=colors,
          ColSideColors=mycols,
          margin=c(10, 10),
          scale = "row")
dev.off()


#PCA for genes


#WEE1, CUTA, AAAS, RFX2 und COX5B
