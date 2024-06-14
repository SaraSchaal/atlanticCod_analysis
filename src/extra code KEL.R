TO DO: Revisit this graph when we have info for the LG group that each sample has.

### Identify clusters in the PCA in the PCA

TO DO: Just do a dendogram
```{r}

# Note: kmeans sometimes does not find the correct # clusters, or it switches the order. Make sure to # set a seed and check the results!
set.seed(8245982)
clust <- kmeans(svd7$u, 6, nstart=20)
str(clust)

table(clust$cluster)

samp_full$LGall_cluster <- as.character(unlist(clust$cluster))

# How pops map to clusters
table(samp_full$LGall_cluster, samp_full$Pop)

head(samp_full)



samp_full$G_LG_PC1 <- svd7$u[,1]
samp_full$G_LG_PC2 <- svd7$u[,2]

ggplot(samp_full) + geom_point(aes(x=G_LG_PC1, y=G_LG_PC2, col = LGall_cluster)) + scale_colour_discrete()

pdf("figures/-clusters.pdf", width=5, height=5)
ggplot(samp_full) + geom_point(aes(x=LG7_PC1, y=LG7_PC2, col = LG7cluster)) + scale_colour_discrete()
dev.off()
```

TO DO: Revisit this plot once we know which samples have which haplotypes

## First heatmap

Based on 5000 loci - showing haplotypes

Use to get an order for plotting
```{r}

set.seed(2345690)
G_LG_sub <- G_LG[,sort(sample(1:ncol(G_LG), 500))] 
head(colnames(G_LG_sub))

#a <- heatmap(G_LG_sub,
#        scale="none",
#        #Rowv = NA,
#        Colv = NA,
#        labRow = sort(samp_full$LG7cluster),
#        cexRow = 0.3,
#        cexCol=0.1)
#head(a$rowInd)

#samp_full$G_LG_heatmapOrder <- a$rowInd
```

## Second heatmap

```{r, eval=FALSE, echo=FALSE}
head(samp_full)



set.seed(2345690)
G_LG_sub1 <- G_LG[,sort(sample(1:ncol(G_LG), 500))] 
head(colnames(G_LG_sub))

my_heatmap <- pheatmap(G_LG_sub1, silent = TRUE)
#do not run str() on this object
G_LG_heatmapOrder <- (my_heatmap$tree_row$order)

heatOrder <- order(samp_full$G_LG_heatmapOrder,
                   #samp_full$Region,
                   #samp_full$Ecotype,
)


G_LG_sub <- G_LG_sub1[heatOrder,]

samp_ann <- data.frame(
  Region = samp_full$Region[heatOrder],
  Pop = samp_full$Pop[heatOrder],
  Ecotype = samp_full$Ecotype[heatOrder]
)
rownames(samp_ann) <- rownames(G_LG_sub) # for pheatmap the annotation has to have the same rownames as the matrix
head(samp_ann)
tail(samp_ann)

head(muts)

chrom_GLG_sub <- unlist(lapply(colnames(G_LG_sub), function(x){strsplit(x,split="__")[[1]][1]}))
head(chrom_GLG_sub)

muts_ann <- data.frame(
  chrom=chrom_GLG_sub
)
rownames(muts_ann) <- colnames(G_LG_sub)


dim(G_LG_sub)

#plot_ord <- 

plot_gap_chrom <- c(
  max(grep(breakpoints$chrom[1],colnames(G_LG_sub))),
  max(grep(breakpoints$chrom[2],colnames(G_LG_sub))),
  max(grep(breakpoints$chrom[3],colnames(G_LG_sub))),
  max(grep(breakpoints$chrom[4],colnames(G_LG_sub)))
)  


annotation_colors = list(
  Region = c(Iceland = "#5977ff", 
             GoM = "#f74747"),
  Pop = c(Pop1 = adjustcolor("darkred",1), 
          Pop2 = adjustcolor("darkred",0.75),
          Pop3 = adjustcolor("darkred",0.5),
          Pop4 = adjustcolor("darkred",0.25),
          Pop5 = adjustcolor("darkred",0.1),
          Pop6 = adjustcolor("darkblue",1),
          Pop7 = adjustcolor("cornflowerblue",0.5),
          Pop8 = adjustcolor("darkblue",0.85),
          Pop9 = adjustcolor("cornflowerblue",0.7)
  ),
  chrom = c( "NC_044048.1" = "black", 
             "NC_044049.1" = "grey50",
             "NC_044054.1" = "grey75",
             "NC_044059.1" = "grey95")
)

# with row clustering
pheatmap(G_LG_sub , 
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         scale = "none",
         #breaks=c(0,1,2),
         show_colnames=FALSE,
         show_rownames=FALSE,
         #cutree_rows = 3,
         annotation_row = samp_ann,
         annotation_col = muts_ann,
         gaps_col = plot_gap_chrom,
         annotation_colors = annotation_colors
)
```