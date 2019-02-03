### inputs 
data.file = './input/data1.xlsx'
human.file = './input/human.gene.rds'
### outputs 

out.dir <- './results/'


### global options 
options(scipen=999)
### shared variable 
library("RColorBrewer")
getPalette = colorRampPalette(brewer.pal(9, "Set1")) # expand color pallete

### load necessary library 
options(java.parameters = "-Xmx11g" )
library("XLConnect")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library(corrplot)
library("ggplot2")
library(gplots)
library("biomaRt")
library("edgeR")
library("pheatmap")
library(ggrepel)
library(reshape)
library(data.table)
library(plyr)
library("sva")
library("limma")
library(heatmap3)
library(gridExtra)
library(Rtsne)
library("sleuth")
library(VennDiagram)
library(knitr)
library(pander)
library(Vennerable)
library(corrplot)
library("psych")
library(survminer)
library ( gridExtra)
library(forcats)
library ( factoextra)
library ( FactoMineR)
library(vegan)
library(kableExtra)
library(survival)
library( stringr )
library ( dendextend)
library(d3heatmap)
library("PerformanceAnalytics")
library("openxlsx")
library (cluster)
require(pathview)
library ( emmeans )
library(lsmeans)
library("openxlsx")
library( stringr )
library( preprocessCore )
library ( ggbeeswarm )




############## functions. 


get.counts <- function (count_dir,ens,sum.this){
  
  data <- data.frame()
  what.choose <- c()
  sampleFiles <- grep('.tab',list.files(count_dir),value=T) # find all the files
  sampleFiles <- sort(sampleFiles)
  # loop through each file and extract only antisense 
  d <- 1;  # this is intitialize the "data" data.frame
  for (file in sampleFiles){
    rnaInfo <- read.table(paste0(count_dir,file), header=FALSE,sep="\t",stringsAsFactors = FALSE)
    rnaInfo <- rnaInfo[5:nrow(rnaInfo),]  # get rid of the first few rows
    rnaInfo$V1 <-gsub('\\.\\d*',"", rnaInfo$V1) # clean up name of gene_id 
    
    #rnaInfo <- rnaInfo[,c(1,4)] # delete extra column. 
    rnaInfo2 <- merge(rnaInfo,ens, by.x="V1", by.y="ensembl_gene_id",all.x=FALSE )
    
    ## take care of transcript counts 
    rnaInfo <- rnaInfo2[!duplicated(rnaInfo2$V1), ]  # ok we need to get rid of duplicates. since some genes can have multiple transcripts. 
    row.names(rnaInfo) <- make.names(rnaInfo$external_gene_name, unique=TRUE) 
    
    id <- gsub("\\ReadsPerGene.out.tab","", file) # filename is ugly but contains the id
    
    ### create intial dataframe with ensembl_id and symbol 
    
    if (d==1){
      data <- data.frame(rnaInfo[ , sum.this ],stringsAsFactors = FALSE)
      colnames(data) <- c(id)
      row.names(data)<-row.names(rnaInfo)
      d <- 0
    }
    
    data[id]<-rnaInfo[, sum.this] # add to column (sample) to data
    print (id)
    what.choose <- c( what.choose, sum ( rnaInfo[, sum.this])/ sum ( rnaInfo$V2) )
    
  }
  
  
  
  
  return( list ( data=data, ratio=what.choose) )  
  
}

# Coefficient of Variation
cov <- function(x){
  x <- sd(x)/ mean ( abs ( x)  ) * 100
  return (x)
}

# data.frame, title, key, name of the column for group (eg control vs drug), name of the col corresponding to dataframe column
# returns a cor-cor boxplot as well as mds
box.outlier <- function (data.boxplot,title.boxplot,key.boxplot, group.id = "targets", tname="name" ){
  
  col <- getPalette ( length (  key [,group.id]  )) 
  
  cor.cpm <- cor( data.boxplot, use="complete.obs", method="pearson")
  cor.cpm <- melt(cor.cpm)
   
  # figour out what samples belongst to what group
  cor.cpm$state <- factor ( 
    as.character ( 
    sapply( as.character ( cor.cpm$X1 ), function(x) 
      key.boxplot [key.boxplot[[tname]] == x, group.id] 
      
    )
  ) )
  

  # order by state - this is your group
  cor.cpm <- cor.cpm[order ( cor.cpm$state), ]
  cor.cpm$X1 <- factor ( cor.cpm$X1, levels = unique (cor.cpm$X1))

  g1<- ggplot(cor.cpm, aes(y=value, x=X1)) +
    geom_violin()+ 
    geom_jitter(shape=19, position=position_jitter(0.07), aes( colour=state), size=3 ) +
    theme_bw() +
    ggtitle(title.boxplot) +
    ylab(" ") +
    xlab("") +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(angle = 90, size=11.5),
          axis.title.x = element_text(size=22),
          
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = .5) +
    scale_colour_manual(values = col )
  
  
  ### plot mds 
  
  mds.semi <-plotMDS( data.boxplot )
  mds.semi <- data.frame(x=mds.semi$x, y=mds.semi$y)
  
  mds.gene <- ggplot( mds.semi , aes(x=x, y=y, shape= key[,group.id] ))  +
    geom_point(size=4) +
    
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Dimension 1", y="Dimension 2") + ggtitle("MDS Supervised")   
  
  mds.gene
  
  
  
  
  
  
  return (list ( g1=g1,cor=cor.cpm, mds=mds.gene) )
  
}

## generates mds plots 
# g1 is for shapes 
# g2 is for color 
mds.pre <- function ( data.mds, title, key.mds, g1, g2, size=5){
  
  mds.semi <-plotMDS( data.mds )
  mds.semi <- data.frame(x=mds.semi$x, y=mds.semi$y)
  
  mds.gene <- ggplot( mds.semi , aes(x=x, y=y, shape= key[,g1], col=key[,g2] ))  +
    geom_point(size=5) +
    
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Dimension 1", y="Dimension 2") + ggtitle(title)   
  
  return ( mds.gene )
  
  
  
}

## draw venn 
### create venn diagrams
## square are circles 
draw.venn <- function ( L1,L2,l1, l2, dow=TRUE, type1="circles", l1.c="#a6cee3", lm.c = "#fdbf6f", l2.c = "#b2df8a"  ){
  l <- list()
  l[[l1]]<- L1
  l[[l2]] <- L2
  g1 <- 
    compute.Venn(
      Venn ( l ) , 
      doWeights = dow, type=type1
    )
  
  gp <- VennThemes(g1)
  
  
  gp[["Face"]][["10"]]$fill <-  l1.c
  gp[["Face"]][["11"]]$fill <-  lm.c # middle
  gp[["Face"]][["01"]]$fill <-  l2.c
  gp$Set[["Set1"]]$col <- "grey"
  gp$Set[["Set2"]]$col <- "grey"
  gp$SetText[["Set1"]]$col <- "black"
  gp$SetText[["Set2"]]$col <- "black"
  grid.newpage()
  plot(g1, gp = gp)
  
  hm <- recordPlot()
  plot.new() ## clean up device
  
  return (hm)
  
  
}

### create histogram of post analysis 
# this function takes a df with 3 columns. 
# 1. p.value 2. logFC and fdr; it then returns a historgram distrubtion for these


hist.limma <- function (result.df){
  
  h1 <- ggplot(data=result.df, aes(result.df$logFC)) +
    geom_histogram(
      binwidth = .2,
      col="grey", 
      fill="#6e9be5", 
      alpha = 1
    ) + 
    labs(title="logFC distrubtion") +
    labs(x="logFC", y="Count")
  
  
  h2 <- ggplot(data=result.df, aes(result.df$P.Value)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#efd37f", 
      alpha = 1
    ) + 
    labs(title="p.value distrubtion") +
    labs(x="p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  h3 <- ggplot(data=result.df, aes(result.df$adj.P.Val)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#bae8a2", 
      alpha = 1
    ) + 
    labs(title="FDR distrubtion") +
    labs(x="corrected p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  ### 
  qfc <- quantile (result.df$logFC, probs= seq(0, 1, by = .05)   )
  qp <- quantile (result.df$P.Value, probs= seq(0, 1, by = .05)   )
  qf <- quantile (result.df$adj.P.Val, probs= seq(0, 1, by = .05)   )
  
  q <- data.frame (lofc=qfc, p.value=qp, fdr=qf)
  
  
  return (list(logfc=h1, p.value = h2, fdr=h3, quantile=q))
}


##Overrepresentation of pathways and Gene Ontology
# r.sub are all the columns name containing the non data fields 
# g1 and g2 column names for are the color and shape from the key
# exp.group is the name of the column that contains the two group
# exp.this and normal.this is what the name of each group is called.
#sample.id are the name of the field containing the ids. 
# GENE_SYMBOL is the field containg the name of your gene. 

#result.gene = df.mmd.hem

#g1 = "group" 
#g2="group" 
#new.key = key.c.mmd.hem
#r.sub = rsub
#exp.group="group"; exp.this="MMD-hem"; normal.this= "Control"; sample.id="ID"; GENE_SYMBOL = "gene"; fdr=.05;  #fold_thres = 0; title1 = "mmd.hem"; top10 =NA



plot.post <- function (result.gene, g1,g2, new.key, r.sub,exp.group="group",exp.this="exp", normal.this="control", sample.id="sample.id",p.val=.05, fdr = .05, fold_thres = 1.5, top10 = NA, GENE_SYMBOL= "GENE_SYMBOL", title1="") {  
  
  new.key$group = gsub ( "-",".", new.key$group)
  exp.this = gsub ( "-",".", exp.this  )
  #### 
  normal.color <- "#9BB3DF"
  exp.color <- "#DCB335"
  new.key$colorCodes <- normal.color
  new.key[new.key[[exp.group]] == exp.this, ]$colorCodes <- exp.color
  ###
  
  
  d.filter <- result.gene[result.gene$P.Value < p.val & result.gene$adj.P.Val < fdr & (result.gene$logFC > fold_thres | result.gene$logFC < -fold_thres )  ,]
  
  mds.semi <-plotMDS( d.filter[,!colnames(d.filter) %in% r.sub  ]  )
  
  mds.semi <- data.frame(x=mds.semi$x, y=mds.semi$y)
  
  mds.gene <- ggplot( mds.semi , aes(x=x, y=y, shape= new.key[,g1], col=new.key[,g2] ))  +
    geom_point(size=8) +
    
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Dimension 1", y="Dimension 2") + ggtitle(paste0(title1, " MDS Supervised") )   
  
  mds.gene
  
  ### heatmap 
  
  my_palette <- colorRampPalette(c("blue", "#f2f2f2", "#fc8600"))(n = 75)
  par(mar=c(0,0,0,0))
  # bottom, left, top and right 

  hm= d3heatmap(scale ( as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ) ), labRow=d.filter$gene, 
  ColSideColors = this.color, distfun = function(x) dist(x,method = 'euclidean'), 
  hclustfun= function(x) hclust(x,method = 'ward.D2')  ) 
  
  
  my_palette <- colorRampPalette(c("blue", "#f2f2f2", "#fc8600"))(n = 75)
  par(mar=c(10,0,0,0))
  # bottom, left, top and right 
  heatmap3(as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ), 
           Rowv=F, dendrogram="column", density.info="none", 
           trace="none", key=TRUE, keysize=1.5, labRow=NA, cexCol=1.0, scale="row", 
           col = my_palette, showRowDendro=TRUE,   margins=c(8,2),
           main= paste0(title1, " Heatmap" ), ColSideColor = new.key$colorCodes, ColSideLabs = "")
  legend("right",legend=c(normal.this,exp.this),fill=c(normal.color,exp.color) )
  
  hm2 <- recordPlot()
  plot.new() ## clean up device
  
  
  
  
  
  #### Scatterplot here 
  
  # scatter plot now
  ## prepping scatter plot   
  
  result.gene$class <- 'no-change'
  result.gene[ result.gene$P.Value < p.val & (result.gene$logFC > fold_thres) & result.gene$adj.P.Val < fdr, ]$class <- 'up'
  if ( nrow ( result.gene[ result.gene$P.Value < p.val & (result.gene$logFC < -fold_thres) & result.gene$adj.P.Val < fdr, ] ) > 0 ){
  result.gene[ result.gene$P.Value < p.val & (result.gene$logFC < -fold_thres) & result.gene$adj.P.Val < fdr, ]$class  <- 'down'
  }
  
  
  ##############
  
  
  overexp_col <- '#fc8600'
  underexp_col <- 'blue'
  nochange_col <- '#F7F6F5'
  
  
  c <- normal.this
  e <- exp.this
  
  # name of average 
  ave.control <- paste0(c,".ave")
  ave.exp <- paste0(e,".ave")
  
  # get all the colname from key sample.id needs to match that of data
  id.control <- new.key[new.key[[exp.group]] %in% normal.this, sample.id]
  id.exp <- new.key[new.key[[exp.group]] %in% exp.this, sample.id]
  
  result.gene[[ave.exp]] <- rowMeans(result.gene[, id.exp  ] )
  result.gene[[ave.control]] <- rowMeans(result.gene[, id.control ] )
  
  result.gene$AVE <- rowMeans( result.gene[, c(  id.control , id.exp )] )
  
  
  
  class <- "class"
  result.gene <- result.gene[order(result.gene$adj.P.Val, result.gene$P.Value),]
  
  if ( is.na (top10) ){
    top10 <- head (result.gene[ !result.gene$class == "no-change" & !result.gene[[GENE_SYMBOL]] == "",  GENE_SYMBOL ],10 )
    top10 <- unique (top10)
  }
  
  label.this <- result.gene[ result.gene[[GENE_SYMBOL]] %in% top10, ]
  
  scatterplot1 <- ggplot(result.gene, aes_string(ave.control ,ave.exp, col="class"), inherit.aes = FALSE )  + 
    scale_color_manual(values=c("up"=overexp_col, "no-change"=nochange_col, "down"=underexp_col), breaks=c("up", "down", "no-change"), labels=c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change")) +
    
    
    geom_point() +
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x=paste0(c ," average (log2 cpm)" ), y= paste0(e ," average (log2 cpm)" ) ) + 
    ggtitle(paste0(title1, "Scatter.plot ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( label.this[[GENE_SYMBOL]])
      ), 
      
      point.padding = unit(.55, "lines"),
      box.padding = unit(2.25, "lines"),
      nudge_y = 0.5,
      fontface=2, 
      color="black"
    )
  
  ### 
  
  
  
  
  ma_p_h <- ggplot(result.gene , aes(x=AVE,y=logFC, col=factor(class) ) )+ 
    scale_color_manual(values=c("up"=overexp_col, "no-change"=nochange_col, "down"=underexp_col), breaks=c("up", "down", "no-change"), labels=c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change")) +
    scale_size_manual(guide="none", values=c("N"=1, "Y"=5)) +
    scale_alpha_manual(guide="none", values=c("N"=0.3, "Y"=1)) +
    #scale_shape_manual(breaks='Y',label='KRAS',values=c("Y"=17,"N"=16))+
    geom_point() +
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Ave Expression", y="FC-log2") +  
    ggtitle(paste0(title1, " MA ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( factor ( label.this[[GENE_SYMBOL]])  )
      ), 
      
      point.padding = unit(2.55, "lines"),
      box.padding = unit(2.25, "lines"),
      nudge_y = 0.5,
      fontface=1, 
      color="black"
    ) 
  
  ### volcano 
  
  
  volcano.g <- ggplot(result.gene , aes(logFC, -log10(P.Value))) +
    geom_point(aes(col=class)) +
    scale_color_manual(values=c(underexp_col, 'black',overexp_col) , name="sig", breaks = c('up','down','no-change'), labels = c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change") ) + 
    ggtitle(paste0(title1, " Volcano.plot ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( label.this[[GENE_SYMBOL]])
      ), 
      fontface=2, 
      color="red",
      
      point.padding = unit(.25, "lines"),
      box.padding = unit(.25, "lines"),
      nudge_y = 0.1
    ) +theme_bw(base_size = 20)  +  theme(legend.text=element_text(size=15)) + guides(color = guide_legend(override.aes = list(size=10)))
  
  
  pm <- list ( 
    mds = mds.gene, 
    hm = hm, hmstat = hm2, 
    scatter = scatterplot1, 
    ma = ma_p_h, 
    volcano = volcano.g, 
    data = result.gene
  )
  
  return (pm)
  
}


### 
## go ontology and kegg
# 1 g.name is the field for the gene name 
# 2 m.name is the filed for the second dataframe containing the entrezgene conversion 
# 3 e.name is the field for entrez name in the annotation frame 
# 2 dataframe. One with class up or down, the field needs to be named as such 
# 3 a dataframe to convert gene name to entrez 
# 4 species. "Hs" (human), "Mm" (mouse), 

get.go <- function (g,mmid,g.name="GENE_SYMBOL", m.name="external_gene_name", e.name="entrezgene", species="Hs",bg, fdr=.05){
  # hg38 must be defined 
  # data frame must contain class to define up and down
  
  g <- merge(g,mmid, by.x=g.name, by.y=m.name )
  
  tm.this <- function (x){
    x$Term <- strtrim(x$Term,60)
    return (x)
  }
  
  up <- unique ( g[g$class == "up",e.name] )
  down <-unique (  g[g$class == "down",e.name] )
  
  # bg is now force to have input so you need to provide a background! 
      bg = merge(data.frame(bg=bg),mmid, by.x="bg", by.y=m.name )
      bg <- unique ( mmid[ , e.name] )

  
  
  go <- goana(list(Up=up, Down=down),
              species=species,
              universe=bg, FDR=fdr)
  
  

  
    # note that I commented out fdr, reason is because fdr already screened above
  
  
  bp <- tm.this ( topGO(go, 'BP', number=Inf) )
  #bp$P.Upfdr = p.adjust(bp$P.Up , n = nrow( bp) ) 
  #bp$P.Downfdr = p.adjust(bp$P.Down , n = nrow( bp) )
  bp = bp[bp$P.Down < .05 | bp$P.Up < .05, ]
  #bp = bp[bp$P.Downfdr < fdr | bp$P.Upfdr < fdr, ]
  bp$P.Up = ifelse ( bp$P.Up < .05, "*", "ns")
  bp$P.Down = ifelse ( bp$P.Down < .05, "*", "ns")
  
  
  
  cc <- tm.this ( topGO(go, 'CC', number=Inf) )
  #cc$P.Upfdr = p.adjust(cc$P.Up , n = nrow( cc) ) 
  #cc$P.Downfdr = p.adjust(cc$P.Down , n = nrow( cc) )
  cc = cc[cc$P.Down < .05 | cc$P.Up < .05, ]
  #cc = cc[cc$P.Downfdr < fdr | cc$P.Upfdr < fdr, ]
  cc$P.Up = ifelse ( cc$P.Up < .05, "*", "ns")
  cc$P.Down = ifelse ( cc$P.Down < .05, "*", "ns")
  
  
  mf <- tm.this ( topGO(go, 'MF', number=Inf) )
  #mf$P.Upfdr = p.adjust(mf$P.Up , n = nrow( mf) ) 
  #mf$P.Downfdr = p.adjust(mf$P.Down , n = nrow( mf) )
  mf = mf[mf$P.Down < .05 | mf$P.Up < .05, ]
  #mf = mf[mf$P.Downfdr < fdr | mf$P.Upfdr < fdr, ]
  mf$P.Up = ifelse ( mf$P.Up < .05, "*", "ns")
  mf$P.Down = ifelse ( mf$P.Down < .05, "*", "ns")
  
  kegg <- kegga(list(Up=up, Down=down),
                species=species,
                universe=bg, FDR=fdr)
  
  kegg = kegg[(kegg$P.Up < .05 | kegg$P.Down < .05 ), ]
  keggUP = kegg[order(kegg$P.Up), ]
  keggDown = kegg[order(kegg$P.Down), ]
  kegg$p = ifelse ( kegg$P.Up < kegg$P.Down, kegg$P.Up, kegg$P.Down)
  # the last one makes sure that p value is not specific to up or down rather it capture both and you can subset later
  keggn = kegg[ order(kegg$p), ]
  keggn$P.Up = ifelse ( keggn$P.Up < .05, "*", "ns")
  keggn$P.Down = ifelse ( keggn$P.Down < .05, "*", "ns")
  
  list ( bp=bp, cc=cc, mf=mf, kegg=keggn #, 
         #k.up = k.up, k.down = k.down
  )
  
}   


### histogram
# this function takes a df with 3 columns. Does not matter how big the df is as long as it has these columns 
# 1. p.value 2. logfc and fdr; it then returns a historgram distrubtion for these
# you must rename the columns how limma does. 
# logFC, P.Value, adj.P.Val
do.histo <- function (result.df){
  
  q.fc <- cut( result.df$logFC, breaks=c(-10,-1.5, -1, 0, 1, 1.5, 10))
  q.fc  <-  data.frame ( table ( q.fc ) ) 
  
  
  q.fdr <-  cut( result.df$adj.P.Val , breaks=c(0.05, .1, .2, .3, .5, .6, 1))
  q.fdr <-  data.frame ( table ( q.fdr) ) 
  q.p <-  cut( result.df$P.Value, breaks=c(0,0.01, .05, 1))
  q.p <-  data.frame ( table ( q.p ) ) 
  
  h1 <- ggplot(data=result.df, aes(result.df$logFC)) +
    geom_histogram(
      binwidth = .2,
      col="grey", 
      fill="#6e9be5", 
      alpha = 1
    ) + 
    labs(title="logFC distrubtion") +
    labs(x="logFC", y="Count")
  
  
  h2 <- ggplot(data=result.df, aes(result.df$P.Value)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#efd37f", 
      alpha = 1
    ) + 
    labs(title="p.value distrubtion") +
    labs(x="p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  h3 <- ggplot(data=result.df, aes(result.df$adj.P.Val)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#bae8a2", 
      alpha = 1
    ) + 
    labs(title="FDR distrubtion") +
    labs(x="corrected p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  return (list(logfc=h1, p.value = h2, fdr=h3, q.fc=q.fc, q.p=q.p, q.fdr=q.fdr ))
}



# melted df 
# pthis ( variable to plot) eg Total.reads
# ds vector of list of samples 
# yl = title of xlab
# color_g color 
# violSize 
gQCviolin <- function(qc_df,pthis,ds,yl, color_g, violSize=15){
  
  
  g1<-  ggplot(qc_df[qc_df$variable==pthis,] , aes(y=value, x=variable)) +
    geom_violin()+ 
    #geom_point(size=1, aes( colour=state2) ) +
    geom_jitter(shape=1, position=position_jitter(0.07), color=color_g ) +
    theme_bw() +
    ggtitle("") +
    ylab("") +
    xlab(yl) +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(), 
          axis.text.y = element_text(size=12),
          #axis.text.x = element_text(size=15),
          axis.title.x = element_text(size=22),
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) +
    
    ##############
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = .5) +
    geom_text_repel(data=qc_df[ qc_df$sample %in% ds & qc_df$variable==pthis,], aes(label=sample ), 
                    arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(.1, 'lines'),color="black", size=violSize   )
  
  
  return ( g1)
  
  
  
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# density plot 
# takes datamatrix cpm log preferbably
denseplot <- function(x,title){
  grid.newpage()
  nsamples <- ncol(y.match)
  col <- getPalette( nsamples)
  par(mfrow=c(1,2))
  mx = median ( lcpm )
  
  
  
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,.25), las=2, 
       main="", xlab="")
  title(main= title, xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", key$tube, text.col=col, bty="n")
  
  hm <- recordPlot()
  plot.new() ## clean up device
  
  return ( hm)
  
}

# takes a qc.df with the following columns. 
# "sample" "Total.reads" "Mapped.reads" "Unmapped.reads"
# "Uniquely.mapped.reads" "Multiple.mapped.reads"
# "Unmapped.reads.Mapped.reads.ratio" "Uniquely.mapped.reads.Multiple.mapped.reads.ratio"
# "Q30.pct" "Sense.anti.sense.ratio"
# "Orientation" "Exonic.Intronic.ratio"
# "Library.size..Mil." "Uniquely.Mapped.NonDupeReads" "dataset"
# then the key 
# its going to color by group so the name of the column.id in the key for the group
qcbar <- function (qc.df, key, group.col="group", tube="tube" ){
  
  
  # ensure that the order is by group

  qc.df$group <- sapply ( qc.df$sample, function(x) key[key$tube == x,group.col ]  )
  
  qc.df <- qc.df[order ( qc.df$group ) , ]
  qc.df$sample <- factor (qc.df$sample, levels = unique (qc.df$sample))
  

  
  colqc <- getPalette ( length (  unique ( qc.df$group )  )) 
  
  Uniquely.Mapped.NonDupeReads <- ggplot(data=qc.df, aes(x=sample , y=Uniquely.Mapped.NonDupeReads, fill=group)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) +
    geom_hline(yintercept = 10000000, color = "red") +
    ggtitle("No Dup & unique mapped, PASS > 10mil") +
    theme(axis.text.x = element_text(angle = 90, size=10) )
  
  

  Mapped.read <- ggplot(data=qc.df, aes(x=sample , y=Mapped.reads, fill= group )) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) +
    geom_hline(yintercept = 35000000, color = "red") +
    ggtitle("Total mapped, Ideal is > 35 million") +
    theme(axis.text.x = element_text(angle = 90, size=10) )
  

  Uniquely.mapped.reads<- ggplot(data=qc.df, aes(x=sample , y=Uniquely.mapped.reads, fill= group )) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) +
    geom_hline(yintercept = 30000000, color = "red")  +
    ggtitle("Unique mapped, Ideal is > 30 million")+
    theme(axis.text.x = element_text(angle = 90, size=10) )
  
  
 
  library.size <- ggplot(data=qc.df, aes(x=sample , y=Library.size..Mil., fill= group )) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) +
    geom_hline(yintercept = 15, color = "red")  +
    ggtitle("Nonscaled Library Size, Ideal is > 15 million")+
    theme(axis.text.x = element_text(angle = 90, size=10) )
  
  unique2multiple <- ggplot(data=qc.df, aes(x=sample , y=Uniquely.mapped.reads.Multiple.mapped.reads.ratio, fill= group )) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) + ylab("mapped: unique to multiple") +
    theme(axis.text.x = element_text(angle = 90, size=10) )
  
  exon2introl <- ggplot(data=qc.df, aes(x=sample , y=Exonic.Intronic.ratio, fill= group )) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = colqc ) + ylab("exon to intron")+
    theme(axis.text.x = element_text(angle = 90, size=10) )
  
  
  
  return ( list (library.size =library.size , Uniquely.mapped.reads=Uniquely.mapped.reads,
                 Mapped.read =Mapped.read , Uniquely.Mapped.NonDupeReads=Uniquely.Mapped.NonDupeReads, 
                 unique2multiple=unique2multiple, exon2introl =exon2introl 
  )  ) 
  
}


### plot by gene 

# subset by gene temp <- v.datap[row.names(v.datap) == "Cyp2a5",]
# key must contain column "tube" matching colomn name of df
# key must contain a column "group" to define each group 
# key.c is redundant as group unless you specify if you do then it could add an additional annotation
# example although you are grouping by column "group" A & B you might want to call it something else
# yl describe what type of data this is 
# example plot.gene ( temp <- v.datap[row.names(v.datap) == "Cyp2a5",], key=key, ptext="",title="")
# add == 1 will add additional column in the beginning for annotations 
## example to add an additional column and annotations 
# plot.gene (  v.datap[row.names(v.datap) == "Cyp2a5",], key=key, ptext="",title="Gapdh looking good",add=1) + annotate(geom = 'text', label = "Based on Limma", x = -Inf, y = Inf, hjust = -.20, vjust = 3, fontface =2) + annotate(geom = 'text', label = "p (a vs b) = .05", x = -Inf, y = Inf, hjust = -.20, vjust = 4.5, fontface =2)+ annotate(geom = 'text', label = "p (c vs d) = .05", x = -Inf, y = Inf, hjust = -.20, vjust = 6, fontface =2)
# key.c  the column that will labe the groups
# pid is the column that you want to label the plot with
plot.gene <- function (temp, key.single, ptext='',title='', yl="log2", add=0, pid="tube", key.c="group",vps= 10, ytitle="log2"  ) {
  
  
  temp <- melt (as.matrix ( temp )  )
  colnames ( temp ) = c("gene","variable","value")
temp$group <- sapply(temp$variable, function(x) key.single[key.single[[pid]] == x, key.c ])
temp$pid <- sapply(temp$variable, function(x) key.single[key.single[[pid]] == x, pid ])
colsingl <- getPalette ( length (  unique ( temp$group )  )) 
# create color scheme
nu <- as.character(factor( temp$group ,labels=colsingl))


temp$group <- factor (temp$group, levels = unique (temp$group))
mean.temp <- temp %>% group_by(group ) %>%
    summarise(
        mean = mean(value)
    ) %>%
    data.frame()



temp$meangroup <- sapply(temp$group, function(x) mean.temp[mean.temp$group == x, ]$mean )
temp$dumm <- "mean"

yp <- max( temp$meangroup)




g <- ggplot(temp, aes(x=factor ( group) , y=value, fill= variable)) +
  geom_bar( stat = "identity", position = position_dodge(), colour='grey' ) +
  scale_fill_manual(values = nu  ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 0, size=15.5),
        axis.text.y = element_text(size=15.5)
  )  +
  scale_colour_manual(name="",values = "black") +
  xlab("Groups") + ylab(yl) + ggtitle(title)  +facet_wrap(gene~., ncol=2, scales = "free") 

g2 = ggplot(temp, aes(x=factor ( group) , y=value )) +
    geom_violin()+ 
    #geom_point(size=1, aes( colour=state2) ) +
    geom_jitter(shape=19, position=position_jitter(0.07), aes( colour=group) , size = vps, alpha= .5) +
    theme_bw() +
    ggtitle( title ) +
    ylab( ytitle) +
    xlab("")  +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(), # comment this out to display cancer.subtype
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=15),
          axis.title.x = element_text(size=22),
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + scale_color_manual(values = unique ( nu)) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = .5) +facet_wrap(gene~., ncol=2) 

bee = ggplot(temp, aes(x=factor ( group) , y=value )) +
    geom_quasirandom(aes(fill = group ), size = vps, shape = 21, alpha = .5) +
    #geom_violin()+ 
    #geom_point(size=1, aes( colour=state2) ) +
    #geom_jitter(shape=19, position=position_jitter(0.07), aes( colour=group) , size = vps, alpha= .5) +
    theme_bw() +
    ggtitle( title ) +
    ylab( ytitle) +
    xlab("")  +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(), # comment this out to display cancer.subtype
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=15),
          axis.title.x = element_text(size=22),
          axis.title.y     = element_text(size=22), 
          legend.text      =element_text(size=12)
    ) + scale_color_manual(values = unique ( nu)) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = .5) +facet_wrap(gene~., ncol=2) 



return ( list ( g=g, violin=g2, bee= bee ) )
}

###
# result.gsva: provide a gsva merged limma matrix
# this.color is the color in the order of the colomns for the sample. eg WT,WT,KO,KO (blue,blue,red,red), gthis to subset
# labfont from .1 to 1 so that the side row fonts can show up
#bottom, left, top and right margins

gsva.hm = function(  result.gsva, this.color , rsub, gthis="KEGG", pv= .05, fdr = .05, fold_thres=.2, 
                     labfont = 1, title="", mg = c(60,450,40,20)){
    result.gene = result.gsva[ grepl ( gthis, result.gsva$gene ) , ]
    
    d.filter <- result.gene[result.gene$P.Value < pv & result.gene$adj.P.Val < fdr & (result.gene$logFC > fold_thres|result.gene$logFC < -fold_thres )  ,]
    
    
    my_palette <- colorRampPalette(c("blue", "#f2f2f2", "#fc8600"))(n = 75)
    
    heatmap3(as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ), 
             Rowv=F, dendrogram="column", density.info="none", cexRow = labfont , 
             trace="none", key=TRUE, keysize=1.5, labRow=d.filter$gene, cexCol=1.0, scale="row", 
             col = my_palette, showRowDendro=TRUE,   margins=c(8,15),
             main= paste0("gsva ", title ), ColSideColor = this.color, ColSideLabs = "")
    
    hm <- recordPlot()
    plot.new() ## clean up device
    
    
    
    #d= d3heatmap(as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ), 
    #ColSideColors = this.color, distfun = function(x) dist(x,method = 'euclidean'), labRow=d.filter$gene, 
    #hclustfun= function(x) hclust(x,method = 'ward.D2')  ) 
    
    df = d.filter[,!colnames(d.filter) %in% rsub  ] 
    row.names ( df ) = d.filter$gene
    d2 = heatmaply( df  , colors = my_palette ,cexRow = labfont ,
                    col_side_colors = this.color, distfun = function(x) dist(x,method = 'euclidean'), 
              hclustfun= function(x) hclust(x,method = 'ward.D2') , margins = mg  ) 
    
    
    return ( d2 )
}


### gsea stuff 

do.gsea = function ( gfold, gname, pathways.in, db=NA, fdr = .05, pv= .05 ){
   
    
    g.this2 = gfold # grab logfc
    names ( g.this2 ) = gname # associate with name of gene 
    g.this2 = sort ( g.this2 ) # probably not necessary but I do it anyway
    
    
    
    # if its mouse then you need a db to convert to human 
    if ( !is.na(db) ){
     
     df = data.frame ( gene=names ( g.this2), logFC = as.numeric ( g.this2) )
     df = merge ( db, df, by.x="MGI.symbol", by.y="gene")
     #df = df[ !duplicated ( df$HGNC.symbol), ]
     g.this2 = df$logFC # grab logfc
     names ( g.this2 ) = df$HGNC.symbol # associate with name of gene 
     
     
    }
    
    gsea.result <- fgsea(pathways.in , g.this2, minSize=15, maxSize=50000, nperm=1000)
    gsea.result = gsea.result[ gsea.result$pval < pv, ] 
    
    gsea.result = gsea.result[ gsea.result$padj < fdr, ] 
    gsea.result = gsea.result[ order ( gsea.result$pval), ]
    
    
    
    kegg.gsea = gsea.result[ grepl("KEGG", gsea.result$pathway), ]
    hall.gsea = gsea.result[ grepl("HALLMARK_", gsea.result$pathway), ]
    react.gsea= gsea.result[ grepl("REACT", gsea.result$pathway), ]
    biocart.gsea= gsea.result[ grepl("BIOCARTA", gsea.result$pathway), ]
    
  
    #kegg.gsea.plot = plotGseaTable(pathways.in[ktemp], g.this2, gsea.result, 
     #                              gseaParam = 0.5) 
    
    return ( list ( 
        kegg = kegg.gsea, 
        hall = hall.gsea, 
        react = react.gsea,
        biocart = biocart.gsea,
        gthis = g.this2
    ))
}



###
### legacy stuff 

qc.all <- function(dataset){
  
  
  
  
  qc_df$UMEM <- as.numeric(qc_df$UMEM )
  
  
  total.sample <- nrow(qc_df)
  qc_df <- melt(qc_df, c("sample","dataset"))
  unique(qc_df$variable)
  qc_df$value <- as.numeric( qc_df$value  )
  
  ### create different column names 
  qc_df$sample2 <- paste0(qc_df$sample,' (',qc_df$value,')')
  qc_df$mil <- paste0(qc_df$sample,' (',round(qc_df$value/1000000,digits=2),')')
  
  # dataset is a vector 
  #dataset <- qc_df[qc_df$dataset %in% 'ASC23', ]$sample
  
  
  g.umem <- gQCgraph(qc_df,"UMEM",dataset,'UMEM','#d8ac11') + geom_hline(yintercept = 10000000, color = "red")
  mr <- gQCgraph(qc_df,"mapped.reads",dataset,'Mapped reads','#f47442')
  unique2mulitple <- gQCgraph(qc_df,"unique.multiple", dataset,'ratio of unique to multiple','#228e0c')
  anti.sense <- gQCgraph(qc_df,"anti.sense", dataset,'ratio of anti to sense','#245ab7')
  exon.intron <- gQCgraph(qc_df,"exon.intron", dataset,'ratio of exon to intron','#a01956')
  library.size <- gQCgraph(qc_df,"library.size", dataset,'Count Library Size','#d38217')
  
  g = list( g.umem=g.umem, mr=mr, unique2mulitple=unique2mulitple, anti.sense=anti.sense,exon.intron=exon.intron, library.size=library.size,qc_df=qc_df )
  
  return (g)
  
}

gQCgraph<- function(qc_df,pthis,ds,yl, color_g){
  
  
  
  mr <- ggplot(qc_df[qc_df$variable==pthis,], aes(x=variable, y=value, color=variable)) + 
    geom_boxplot(outlier.colour="grey", outlier.shape=20,
                 outlier.size=4) + 
    scale_colour_manual(values=c(color_g)) +
    theme(legend.position="none", 
          legend.title=element_blank(), 
          legend.key = element_blank(),
          # axis.text.x = element_text(angle = 0, vjust=1, hjust = 0.1, size=11.5)
          axis.text.x = element_text(angle = 0, size=11.5), 
          axis.text.y = element_text(size=7) 
    ) +
    labs(x = "", y="" ) +
    ggtitle( "" ) +
    geom_text_repel(data=qc_df[ qc_df$sample %in% ds & qc_df$variable==pthis,], aes(label=sample ), 
                    arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(.25, 'lines'),color="#8ba0e8"   )
  return (mr)
  
}

# add ability to decode sig GO terms to DGE genes

# main is the GO output 
# df is all the DGE genes 
# this function requires a mouse conversion df 

getGOgenes.mouse <- function( main, df, mouse ) {
    
    #main = GO.d3.gfp.ank$bp
    goterm = row.names ( main )
    
    df =  df$gene
    df =  mouse[mouse$external_gene_name %in% df, ]$entrezgene
    df = unique ( df )
    
    gogene = c()
    
    for ( go in goterm){
        
        x <- org.Mm.egGO2ALLEGS # org.Hs.egPATH2EG
        Rkeys(x) <- go
        EG <- mappedLkeys(x)
        
        igene = intersect ( EG, df)
        # convert it back 
        igene = unique ( mouse[mouse$entrezgene %in% igene, ]$external_gene_name )
        gogene = c(gogene, toString(igene))
        
    }
    
    main$gene = gogene # add genes to goterm df 
    
    return ( main )
}



### specific to this analysis. 

# function plots violin plots for each sample  ( column ) in data.frame
# you need to provide a key where there is a column ID matching to a column called group, this will color code the group

boxv <- function ( temp, key.temp, title="" ){
    
    temp = melt ( temp )
    colnames(temp) = c("ID", "value")
    temp$group = as.character ( sapply ( temp$ID,  function(x)  unique ( key.temp [ key.temp$ID  == x, ]$group)  ) )    
    g1 = ggplot(temp, aes(y=value, x=ID)) +
        geom_violin()+ 
        geom_jitter(shape=19, position=position_jitter(0.07), aes( colour=group), size=1, alpha = .5 ) +
        theme_bw() +
        ggtitle( title ) +
        ylab(" ") +
        xlab("") +
        theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(),
              
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(angle = 90, size=11.5),
              axis.title.x = element_text(size=22),
              
              axis.title.y     = element_text(size=22), 
              legend.text      =element_text(size=12)
        ) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                         geom = "crossbar", width = .5)
    return ( g1 )
}


