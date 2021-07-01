### Functions
## prepare counts - for TEtranscripts
prepcounts <- function(pathfile,patternsam){ 
  pathfile <- pathfile 
  setwd(pathfile)
  # all counts
  filenames = list.files(pathfile,pattern = "*.cntTable")
  datalist = lapply(filenames, function(x){read.csv(file=x,sep = "\t")})
  data <- Reduce(function(x,y) {merge(x,y)}, datalist)
  rownames(data) <- data$gene.TE
  setwd(wkdirpath)
  # modify sample names
  colnames(data) <- gsub(patternsam,"",x=colnames(data))
  colnames(data) <- gsub(suffixbam,"",x=colnames(data))
  # reorder <- sort(as.numeric(colnames(data)),decreasing = F)
  # reorder
  # data <- data[c(as.character(reorder))]
  # pre-filtering
  min_read <- 1
  data <- data[apply(data,1,function(x){max(x)}) > min_read,]
  data$gene.TE <- NULL
  data <- data[,order(colnames(data))]
  return(data)
}
## signatures/pathways
organisePATHS <- function(selectPATH, outPATH, file, pathways){
  if( is.null(pathways) ){
    # prepare pathway data
    paths = msigdbr(species = "Homo sapiens" , category = "H")  # no need to download the file / "HALLMARK"
    paths = paths %>% 
      split(x = .$gene_symbol, f = .$gs_name)
    pathways <-  paths
  }
  # read file
  stem <- read.table(file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  # get index for selected signature
  index <- which(colnames(stem) == selectPATH)
  # identify NAs
  stem$name <- ifelse( stem[[index]] != "", stem[[index]], NA)
  # remove NAs
  name <- subset(stem$name,!is.na(stem$name))
  # get gene names
  name <- name[2:length( name )]
  # include new signature to pathway data
  pathways$default <- name
  index2 <- which(names(pathways)=="default")
  names(pathways)[index2] <-  outPATH
  pathways$default <- NULL
  return(pathways)
}
## volcano plot with annotation of most significant genes
volcano_plot <- function(data,log2,threshold,tag){
  resdata<-as.data.frame(data)
  resdata$gene<-rownames(data)
  resdata<-resdata[complete.cases(resdata),]
  resdata$Legend <- ifelse(resdata$padj < 0.05 & resdata$log2FoldChange > log2, "Up",#,paste0("upregulated",log2),
                          ifelse(resdata$padj < 0.05 & resdata$log2FoldChange < (log2 * -1), "Down", #, paste0("downregulated",log2),
                                 ifelse(resdata$padj > 0.05 & abs(resdata$log2FoldChange)>log2, paste0("No Sig FC > ",log2),
                                        "No Sig"
                                 )))
  # top significant Genes
  top_genes <- rbind(resdata[resdata$Legend=="Down",],
                     resdata[resdata$Legend=="Up",])
  # order by padj values, lower on the top
  top_genes <- top_genes[order(top_genes$padj,decreasing = F),]
  # extract same number of upregulated and downregulated gene values
  top_genes <- rbind(head(top_genes[top_genes$log2FoldChange > 0,],5),
                     head(top_genes[top_genes$log2FoldChange < 0,],5))
  # prepare gene annotation
  resdata$label <- ifelse(resdata$gene %in% top_genes$gene,resdata$gene,"")
  #resdata$label <- na.omit(resdata$label)
  toplog2 <- c(subset(resdata$gene,resdata$log2FoldChange == max(resdata$log2FoldChange)),
               subset(resdata$gene,resdata$log2FoldChange == min(resdata$log2FoldChange)))
  toplog2pos <- c(which(resdata$log2FoldChange == max(resdata$log2FoldChange)),
               which(resdata$log2FoldChange == min(resdata$log2FoldChange)))
  resdata$label[toplog2pos] <- c(toplog2)
  cols <- c("No Sig" = "darkgrey",
            "Down" = "#25a2d0",
            "Up" = "#d91229",
            "No Sig FC > log2" = "#009e73"
            )
  # indicate selected log2FC in graph
  names(cols)[4] <- paste0("No Sig FC > ",log2)
  ## Make a basic ggplot2 object
  vol <- ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj), colour = Legend))
  # inserting manual colors as per color pallette  with term "scale_colour_manual(values = cols)" below
  
  print(vol +   
    #ggtitle(label = paste0("Volcano Plot: ", selection)) +
    geom_point(size = 1.5, alpha = .99, na.rm = T)  +
    geom_text_repel(data=(resdata ), aes(x=log2FoldChange, y=-log10(padj),label=label),
                    colour="black")+#,point.padding = unit(1.3, "lines"))+
    scale_color_manual(name="Differential Expression",values=cols)+
    theme_bw(base_size = 10) + 
    theme(legend.position = "right") + 
    xlab(expression(log[2](expression))) + # X-Axis label
    ylab(expression(-log[10]("padj"))) + # Y-Axis label
    geom_hline(yintercept = -log10(0.05), # horizontal dashed cut-off line (padj)
               colour="#990000", linetype="dashed") + 
    geom_vline(xintercept = (log2), colour="#990000")+ # upreguated genes cut-off line
    geom_vline(xintercept = -(log2), colour="#990000")+ # downreguated genes cut-off line
    scale_y_continuous(trans = "log1p",breaks=c(1,10,100,1000))+
    scale_x_continuous(limits = c(quantile(resdata$log2FoldChange)[1][[1]]-1, quantile(resdata$log2FoldChange)[5][[1]]+1))+
    labs(subtitle = tag)
  )
}
## colours - heatmap
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)] #https://slowkow.com/notes/pheatmap-tutorial/
}
## Heatmap: for variation
heatmap_var <- function(resobj,tag,col){ # delete tag here and in calling functions
  ids <- rownames( resobj )
  matVar_ids <- assay(vst(dds))[ ids , ]
  topVar <- head(order(rowVars( matVar_ids), decreasing = T), 50)
  topVar <- rownames(matVar_ids[ topVar ,])
  mattopVar <- assay(vst(dds, blind = F))[ topVar , ]
  mattopVar <- mattopVar - rowMeans( mattopVar )
  # colours
  # mat_breaks <- quantile_breaks(mattopVar, n = 11) 
  return(pheatmap( mattopVar , scale="row", clustering_distance_rows="correlation", 
                  clustering_distance_cols="correlation",
                  fontsize = 7, fontsize_row = 7, height=20, 
                  annotation = annot,
                  color = col,#dichromat(type = "tritan",viridis(length(mat_breaks) - 1)),#col,
                  angle_col = "315",
                  drop_levels = T,
                  annotation_colors = mycolors
                  ))
  #breaks = mat_breaks))
}
# # Heatmap: for most significant genes/tes
# heatmap_varSig <- function(resobj,filters,tag){
#   if(filters == "log2FC"){
#     namepdf <- paste0(contrast[2],"vs",contrast[3],"_heatVar_topSigFC_",tag,".pdf")
#     top_DEGref_up <- resobj[order(resobj$log2FoldChange,decreasing=TRUE),]
#     top_DEGref_down <- resobj[order(resobj$log2FoldChange,decreasing=FALSE),]
#     top_DEGref <- rbind(head(top_DEGref_up,25), # upregulated
#                         head(top_DEGref_down,25)) # downregulated
#   }else{
#     namepdf <- paste0(contrast[2],"vs",contrast[3],"_heatVar_topSig_",tag,".pdf")
#     top_DEGref <- resobj[order(resobj$padj,decreasing=FALSE),]
#     top_DEGref <- rbind(head(top_DEGref[ top_DEGref[, 'log2FoldChange'] > 0 , ],25), # upregulated
#                         head(top_DEGref[ top_DEGref[, 'log2FoldChange'] < 0 , ],25)) # downregulated
#   }
#   idref <- rownames( top_DEGref )
#   matSig <- assay(vst(dds, blind = F))[ idref, ]
#   matSig <- matSig - rowMeans( matSig )
#   # colours
#   mat_breaks <- quantile_breaks(matSig, n = 11) 
#   # plot
#   pdf(namepdf)
#   print(pheatmap( matSig , scale="row", clustering_distance_rows="correlation", 
#                    clustering_distance_cols="correlation",
#                    fontsize = 7, fontsize_row = 7, height=20, 
#                    main = paste(contrast[2],"vs",contrast[3],"by",filters, sep = " "),
#                    annotation = annot,
#                    color = col,#dichromat(type = "tritan",viridis(length(mat_breaks) - 1)),#col,
#                   drop_levels = T,angle_col = "315"))#,annotation_names_row = T))#dichromat(type = "tritan",rainbow(length(mat_breaks) - 1)),
#                   #breaks = mat_breaks))
#   dev.off()
# }

## Heatmap: for samples pairwise correlation
heatmap_paircor <- function(resobj, filters, tag){
  # colours
  heat.colors <- brewer.pal(6, "Blues")
  if(!is.null(filters)){
    namepdf <- paste0(contrast[2],"vs",contrast[3],"_heatPairCorr_topSig_",tag,".pdf")
    top_DEG <- resobj[order(resobj$padj,decreasing=FALSE),]
    top_DEG <- rbind(head(top_DEG[ top_DEG[, 'log2FoldChange'] > 0 , ],25), # upregulated
                          head(top_DEG[ top_DEG[, 'log2FoldChange'] < 0 , ],25)) # downregulated
    ids <- rownames( top_DEG )
    titlemain <- "Top 50"
  }else{
    namepdf <- paste0(contrast[2],"vs",contrast[3],"_heatPairCorr_",tag,".pdf")
    ids <- rownames( resobj )
    titlemain <- "All Sig"
  }
  mat_ids <- assay(vst(dds))[ ids, ] 
  cor_ids <- cor( mat_ids , method = "spearman") # pairwise correlation values for samples
  # plot
  return(pheatmap( cor_ids , color = heat.colors, border_color=NA, 
                   fontsize = 7, fontsize_row = 7, height=20,
                   annotation = annot, main = paste(titlemain,tag),
                   angle_col = "315",
                   drop_levels = T,
                   annotation_colors = mycolors))
}
## Phylo
phylo <- function(resultobj,tag){
  ids <- rownames( resultobj )
  mat_ids <- assay(vst(dds))[ ids, ]
  d <- cor( mat_ids , method = "spearman")
  hc <- hclust(dist(1 - d))
  pdf(paste0(contrast[2],"vs",contrast[3],"_dendrogram_",tag,".pdf"))
   print(ape::plot.phylo(as.phylo(hc), type = "p", edge.col = 4, edge.width = 1, cex = 0.7, 
             show.node.label = TRUE, no.margin = TRUE,font = 1))
  dev.off()
}
## Enrichment analysis
enrich <- function(resobj){
  resobj$symbol <- rownames(resobj)
  if( squire!="" ){
  resobj$symbol <- str_remove(res$symbol, ",.")
  }
  res_statFrame <- data.frame(resobj)
  # prepare data to generate a list of entrezID and stat values
  res_enrich <- res_statFrame %>% 	
  dplyr::select(symbol, stat) %>% 
  distinct() %>% 
  group_by(symbol) %>% 			  # group entrez-stats values for each gene
  dplyr::summarize(stat=mean(stat))	# decreasing ties in the preranked stats
  # input for Enrichment analysis
  geneList <- data.frame(res_enrich)
  # creates list from data frame
  ranks <- deframe(geneList)    
  # sorting by stats values (it can be sorted by log2FC)
  ranks_sort <- sort(ranks, decreasing = T)
  
  ## runs GSEA
  fgseaRes <- fgsea(pathways=paths, stats=ranks_sort, nperm = 10000)
  # plots non-collapsed pathways NES vs padj 
  pdf(paste0(contrast[2],"vs",contrast[3],".pdf"))
  limit <- c(0,1)#*c(-1, 1)
  print(ggplot(fgseaRes, aes(reorder(pathway, NES), NES, fill = padj)) +     
        geom_col() +
        scale_fill_distiller(type = "div", limit = limit)+
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title=paste0(contrast[2]," vs ",contrast[3])) + # write title
        theme_minimal() )
  dev.off()
}
## PCA
pcaplot <- function(trans,sampleselect,listcond){
  if(ncol(trans)>16){
  return(plotPCA(trans, intgroup = c(listcond)))
  }else{
  # another
  pobj <- pca(assay(trans), metadata = sampleselect, removeVar = 0.1)
  return(biplot(pobj,
               colby = complcondit[2], #"condition",
               colLegendTitle = paste('PCA:',project),
               # encircle config
               encircle = TRUE,
               encircleFill = TRUE,
               hline = 0, vline = c(-25, 0, 25),
               legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0))
        }
}
## Save csv files
savecsv <- function(resobj,tag,LFC,selection){ 
  # order by log2FC
  top_resobj <- resobj[order( resobj$log2FoldChange, decreasing = TRUE),]
  # extract results for both TEs and Genes
  top_resobj <- top_resobj[rownames( top_resobj ),]
  # write csv files
  if(LFC == ""){
  write.csv(top_resobj, file = paste0("res_",selection,"_",tag,".csv"))
  }else{write.csv(top_resobj, file = paste0("res_",selection,"_LFC_",tag,".csv"))}
}
## Prepare input for DESeq2
prepare_DEA <- function(){
  # counts
  data <- read.table( countmatrix , header = T,row.names = 1,sep = countSEP)
  # metadata
  sampleInfo <- read.csv(file.path(wkdirpath, metafile ), header = T, row.names = 1) # modify by read.table() and sep = metaSEP
#  sampleInfo <- read.table(file.path(wkdirpath, metafile ), header = T, row.names = 1, sep = metaSEP)
  ## DEA - check comparison
  if(SELECTION_COND_2!=""){
    indexselect <- colnames(sampleInfo[which(colnames(sampleInfo)==SELECTION_COND_2)])
    group <- list(levels(sampleInfo[[indexselect]]), project) # for DEA: selection
    group <- unlist(group)
  }else{
    group <- project
  }
  return(list(data=data,sampleInfo=sampleInfo,indexselect=indexselect,group=group))
}
## Run DESeq2
ddsobj <- function(){
    sampleInfo <- inputdea$sampleInfo
    indexselect <- inputdea$indexselect
    # set wprking directory, ensure main path
    setwd(wkdirpath)
    # each group-type
    selection <- inputdea$group[[j]] 
    # checking
    print(selection) 
    # create new directory
    dir.create(selection)    
    ## Select samples with if statement
    if( selection != project ){ # for ALL
      samples <- rownames(sampleInfo[which(sampleInfo[[indexselect]] == selection),])
    }else{
      samples <- rownames( sampleInfo )
    }
    ## Set new working directory
    # where to save output
    setwd(file.path(wkdirpath,selection))
    ## Filter samples
    # extract selected samples
    data_sub <- data[,samples]
    sampleInfo_sub <- sampleInfo[samples,]
    # annotations
    annot <- sampleInfo_sub[complcondit]
    # fix colours of main condition in Heatmaps
    names(mycolors$condition) <- unique(annot$condition)
    ## DEA
    # run DESeq2, counting differences by collection media too
    if( length(contrasts_formula) >1 ){
      #designform <- unlist(paste(complcontit_pca[2],complcontit_pca[1],sep = "+"))
      dds <- DESeqDataSetFromMatrix(countData = data_sub , colData = sampleInfo_sub,
                                    selectdesign)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = data_sub , colData = sampleInfo_sub,
                                    design = ~condition)
    }
    dim(dds)
    # filter counts
    keep <- rowSums(counts(dds)) >= 10
    # remove genes with <10 counts
    dds <- dds[keep, ] 
    # run differential expression analysis with DESeq2
    dds <- DESeq(dds,parallel = TRUE)    
    ## Transform data
    # used by Phylo and Sample-pairwise analyses
    vds <- vst(dds, blind = T ) # for PCA
    mat <- assay(vst(dds, blind = F)) # for heatmap   
    ## PCA
    pca <-  pcaplot(trans = vds,sampleselect = sampleInfo_sub,listcond = complcontit) #complcontit_pca
	## Output
    return(list(dds=dds,vds=vds,mat=mat,pca=pca,selection=selection,annot=annot,mycolors=mycolors))
}
## Graphs DEA
rundea <- function(pca,selection){
  # print the current contrast formula, it is just a feedback for the user
  print(paste(selection,contrast[2],"vs",contrast[3]))
  cont <- paste0(contrast[2],"vs",contrast[3])
  condit <- paste0(contrast[2],"/",contrast[3])
  # just checking if the pdf file already exists in the specified working directory
  # if it does not exist, create it and plot the graphs
  if(!file.exists(paste0(contrast[2],"vs",contrast[3],".pdf"))){
    print("results dds")
    # allow filtering for low counts, outlier genes, etc, using independentFiltering=T
    res <- results(dds,independentFiltering=T, alpha = 0.05, contrast = contrast, parallel = T)
    res <- na.omit(res)
    enrich(res)
   
    ## preparing data to generate volcano plot
    # shrinkage of Log2FC, correcting log2FC for low gene counts
    resLFC <- lfcShrink(dds, res = res, type = "ashr", contrast = contrast, parallel = T) 
    # transposable elements
    TE_only <- resLFC[grepl(":",rownames(resLFC)),]
    # genes results
    RefSeq_only <- resLFC[!rownames(resLFC) %in% rownames(TE_only),] #RefSeq_only <- subset(RefSeq_only, baseMean > 10)
    # lines
    LINE_only <- TE_only[grepl(":LINE",rownames(TE_only)),]
    
    ## VOLCANO PLOT
    # for TEs/LINES/GENES
    tvol <- volcano_plot(data = TE_only, log2 = cutlog2FC_tes, threshold = nname_tes, tag = tag_tes)
    lvol <- volcano_plot(data = LINE_only, log2 = cutlog2FC_ls, threshold = nname_ls, tag = tag_ls)
    gvol <- volcano_plot(data = RefSeq_only, log2 = cutlog2FC_g, threshold = nname_g, tag = tag_g)
    # extract legend
    legend <- get_legend(tvol)
    tvol <- tvol+ theme(legend.position="none")
    lvol <- lvol+ theme(legend.position="none")
    gvol <- gvol+ theme(legend.position="none")
    # empty graph
    blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
      cowplot::theme_nothing()
    # plot
    print("Volcano plot & QC")
    png(paste0(cont, "_Volcano.png"),width = 3100,height = 2000, res=250)
    grid.arrange(tvol,lvol,legend,gvol,pca,blankPlot,ncol=3, nrow =2,widths=c(2.3, 2.3 ,0.6))
    dev.off()
    
    ## PLOT MA: looking for log2FC and mean of counts representation
    # useful to fix a logs2FC cutoff
    # for TEs/LINES/GENES
    png(paste0(cont, "_MA.png"),width=1400, height=1200, res=250)
    par(mfrow = c(1,3))
    tma=DESeq2::plotMA(TE_only,main=tag_tes)
    lma=DESeq2::plotMA(LINE_only,main=tag_ls)
    gma=DESeq2::plotMA(RefSeq_only,main=tag_g)
    dev.off()
    
    ## DEg/tes
    DEGref <- RefSeq_only[which(RefSeq_only$padj < 0.05),]
    DEGtes <- TE_only[which(TE_only$padj < 0.05),] 
    DEGlines <- LINE_only[which(LINE_only$padj < 0.05),] 
    
    ### Heatmaps 
    print("Heatmaps")
    # performs hierarchical clustering on the vsd transformed expression matrix 
    # subsetted by the DEG/TEs identified in the above differential expression analysis
    ## - Pearson correlation-based distance measure and complete linkage for cluster joining
    gv=heatmap_var(DEGref,tag_g,colheat) # most variable genes
    tv=heatmap_var(DEGtes,tag_tes,colheat) # most variable tes
    
    ## - Spearman pairwise correlation
    gp=heatmap_paircor(DEGref,filters = NULL,tag_g) # all significant genes
    gps=heatmap_paircor(DEGref,filters = "yes",tag_g) # top 50 DE significant genes
    # plots
    png(paste0(cont, "_heatPairCorr",tag_g,".png"),width=3300, height=1700, res=300)
    gridExtra::grid.arrange(grobs=list(gp[[4]], gps[[4]]), 
                            ncol= 2 )
    dev.off()
    ## TEs
    tp=heatmap_paircor(DEGtes,filters = NULL,tag_tes)  # all significant tes
    tps=heatmap_paircor(DEGtes,filters = "yes",tag_tes) # top 50 DE significant tes
    # plots
    png(paste0(cont, "_heatPairCorr",tag_tes,".png"),width=3300, height=1700, res=300)
    gridExtra::grid.arrange(grobs=list(tp[[4]], tps[[4]]), 
                            ncol= 2)
    dev.off()
    
    ## same for LINEs
    if(length(rownames(DEGlines)) >= 10 ) {
      lv=heatmap_var(DEGlines,tag_ls,colheat) # most variable lines
      lp=heatmap_paircor(DEGlines, filters = NULL,tag_ls) # all significant lines
      lps=heatmap_paircor(DEGlines,filters = "yes",tag_ls) # top 50 DE significant lines
      # plots
      png(paste0(cont, "_heatVar_L.png"),width = 1900, height = 1800, res = 300)
      gridExtra::grid.arrange(grobs=list(lv[[4]]), 
                              ncol= 1)#+guides(fill=guide_legend(nrow=2,byrow=TRUE))## display plot
      dev.off()
    }
    if(length(rownames(DEGtes)) >= 10 ){
      # plots
      #png(paste0(cont, "_heatVar_GT.png"),width=3300, height=1900, res=300)
      png(paste0(cont, "_heatVar_GT.png"),width = 3300 ,height = 1900 ,res = 300)
      gridExtra::grid.arrange(grobs=list(tv[[4]], gv[[4]]), 
                              ncol= 2)#+guides(fill=guide_legend(nrow=2,byrow=TRUE))## display plot
      dev.off()
    }
    
    ## Dendrogram - Spearman correlation
    print("Dendrograms")
    phylo(RefSeq_only,tag_g)
    phylo(TE_only,tag_tes)
    phylo(LINE_only,tag_ls)
    print(paste("End analysis for", selection,cont))
    
    ## Save output for each comparison
    savecsv(DEGref,tag_g,"",selection)
    savecsv(DEGref,tag_g,"LFC",selection)
    savecsv(DEGtes,tag_tes,"",selection)
    savecsv(DEGtes,tag_tes,"LFC",selection)
  }
}
metadata
