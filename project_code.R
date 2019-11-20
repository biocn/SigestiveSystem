setwd('E:/work_backup/201909_Gastrointestinal_Pan-Adenocarcinomas/source')

##########base info##############
m6a.genes=c('METTL3','METTL14','WTAP','KIAA1429','RBM15','RBM15B'
            ,'ZC3H13','RBMX','YTHDC1','YTHDC2','YTHDF3','YTHDF1','YTHDF2'
            ,'HNRNPC','HNRNPA2B1','IGF2BP1','IGF2BP2','IGF2BP3','FTO','ALKBH5')
cancer.code=c('STAD','READ','PAAD','LIHC','ESCA','COAD','CHOL')
library(hgu133plus2.db)
m6a.gene.id=select(hgu133plus2.db,keys=m6a.genes,keytype = 'SYMBOL',columns = c('ENTREZID','ENSEMBL'))
m6a.gene.id[4,2]=25962
m6a.gene.id[4,3]='ENSG00000164944'
m6a.gene.id=m6a.gene.id[-10,]
############m6a exp################
#for(abbr in cancer.code){
#  tpm=read.csv(paste0('exp/',abbr,'.tpm.txt'),sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
  #print(c(abbr,length(grep('-01$',colnames(tpm))),length(grep('-11$',colnames(tpm)))))
#}
library(corrplot)
library(psych)

for (abbr in cancer.code) {
  data=read.csv(paste0('exp/',abbr,'.tpm.txt'),sep = '\t',check.names = F,row.names = 1,stringsAsFactors = F)
  m6a_gene<-as.matrix(data)
  M=cor(t(data),method = "pearson")
  col3 <- colorRampPalette(c("#1f6650", "white", "#ea5e5e")) 
  res1<-cor.mtest(t(m6a_gene),conf.level = .95)
  pdf(file=paste0('Fig2A_',abbr,'_cor.pdf'),width = 6,height = 6)
  corrplot(M,p.mat = res1$p, insig = "label_sig",
           sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white", type = "upper",col=col3(100),tl.cex=0.8)
  dev.off()
}


all.fc=cbind()
all.ps=cbind()
vd1.sbs=rbind()
for (abbr in cancer.code) {
  data=read.csv(paste0('exp/',abbr,'.tpm.txt'),sep = '\t',check.names = F,row.names = 1,stringsAsFactors = F)
  normal=grep('-11$',colnames(data))
  tumor=grep('-01$',colnames(data))
  sp=c(rep('Tumor',length(tumor)),rep('Normal',length(normal)))
  dat.al=rbind()
  fcs=c()
  ps=c()
  for(i in 1:nrow(data)){
    fcs=c(fcs,mean(as.numeric(data[i, tumor]))/mean(as.numeric(data[i, normal])))
    ps=c(ps,t.test(as.numeric(data[i, tumor]),as.numeric(data[i, normal]))$p.value)
    dat=data.frame(Expression=as.numeric(data[i, c(tumor,normal)])
                   ,category=sp,type=rep(row.names(data)[i],length(sp)))
    dat.al=rbind(dat.al,dat)
  }
  all.fc=cbind(all.fc,fcs)
  all.ps=cbind(all.ps,ps)
  vd1.sbs=rbind(vd1.sbs,cbind(dat.al,fa=rep(abbr,nrow(dat.al))))
}
colnames(all.fc)=cancer.code
row.names(all.fc)=row.names(data)
colnames(all.ps)=cancer.code
row.names(all.ps)=row.names(data)

up=-log10(all.ps)
down=log2(all.fc)
UpColor <- colorRamp2(breaks = c(0, 5), colors = c("#FFFFFF","#eafbea"))
DnColor <- colorRamp2(breaks = c(-1,0, 4), colors = c("#0000FF",'#FFFFFF',"#ea5e5e"))

row_an <-  HeatmapAnnotation(type = m61.genes.type, 
                             show_annotation_name = F, 
                             col = list(type = c("Writers" = "#5AC9FA", "Readers" = "#FAC67A", "Erasers" = "#51B743")), 
                             show_legend = T,  
                             annotation_legend_param = list(title = "m6A group",nrow = 1), 
                             which = "row" 
)

DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                 unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
    
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                 unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    if(up[i, j]>=1.3){
      txt="****"
      if(up[i, j]>=1.3&up[i, j]<2){
        txt='*'
      }else if(up[i, j]>=2&up[i, j]<3){
        txt='**'
      }else if(up[i, j]>=3&up[i, j]<4){
        txt='***'
      }
      grid.text(label=txt,x=(x + 0.5*width),
                y=(y+ 0.5*height),just = c('right','top'))
    }
  }
}

p1 <- Heatmap(up, column_title = "Gene Expression across Tumor of digestive system"
              , rect_gp = gpar(type = "none")
              , show_heatmap_legend = F
              , cluster_rows = F
              , cluster_columns = T,
              left_annotation = row_an, 
              cell_fun = DiagFunc(up = up, down = down) 
) 
col_fun = colorRamp2(c(-1, 0, 4), c("#0000FF", "#FFFFFF", "#ea5e5e")) 
lgd <- Legend(title = "log2(fold change)", 
              col_fun = col_fun, 
              at = c(-1,0,1,2,4), 
              labels = c("-1","0","1","2",">4"),  
              direction = "horizontal" 
)
col_fun2 = colorRamp2(c(0, 5), c("#FFFFFF", "#eafbea")) 
lgd2 <- Legend(title = "-log10(p value)", 
               col_fun = col_fun2, 
               at = c(0,1,2,3,4,5),
               labels = c('0',"1","2","3","4",">5"),  
               direction = "horizontal" 
)
pdf('Fig2B.pdf',width = 8,height = 6)
draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom",heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()
###########case9###############
library(ggplot2)
library(tidyverse)
library(reshape2)
library(optparse)
plotMutiBar=function(dat,ist=F,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0){
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  df=data.frame(bg=paste0('R',1:nrow(dat)),dat)
  colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  pg=ggplot(melt(df), aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(df)-1)) {
      tmp=df[order(df[,1],decreasing = T),]
      tmp[,i]=cumsum(tmp[,i])
      tmp[,i+1]=cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      pg=pg+geom_segment(data=data.frame(tmp,STX=rep(i-1+lineW/2,nrow(tmp)),EDX=rep(i-lineW/2,nrow(tmp))), aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  if(showValue){
    pg=pg+geom_text(data=melt(df),aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  pg=pg+labs(x=xlb, y=ylb)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1))
  }
  return(pg)
}

data<-read.table('exp/cris9.txt',header =F ,stringsAsFactors = F)
m6a.cirs=data[match(m6a.genes,data[,1]),]
essential=apply(m6a.cirs[,2:4], 1, mean)
non.essential=100-essential
essential[is.na(essential)]=0
non.essential[is.na(non.essential)]=0
dat=data.frame(essential,non.essential)
row.names(dat)=m6a.genes
pdf('Fig2C.pdf',width = 8,height = 4)
plotMutiBar((t(dat)),showLine=F,showValue = F)
dev.off()

##################GSVA############
library('GSVA')
library(GSEABase)
c2KEGG <- getGmt("GSVASource/c2.cp.kegg.v6.2.symbols.gmt",
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=SymbolIdentifier())
for(code in cancer.code){
  sy.exp=read.csv(paste0('exp/',code,'.symbol.tpm.txt'),sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
  ssGSEA.kegg <- gsva(as.matrix(sy.exp), c2KEGG,method='ssgsea',
                      min.sz=10, max.sz=500, verbose=TRUE)
  save(ssGSEA.kegg,file = paste0(code,'.ssGSEA.kegg.RData'))
}

all.sig.m6a.kegg=rbind()
for(code in cancer.code){
  load(file = paste0(code,'.ssGSEA.kegg.RData'))
  row.names(ssGSEA.kegg)=gsub('KEGG_','',row.names(ssGSEA.kegg))
  exp=read.csv(file=paste0('exp/',code,'.tpm.txt'),sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
  m6a.kegg.cor=cor(t(exp),t(ssGSEA.kegg))
  all.sig.r=rbind()
  for(i in 1:nrow(m6a.kegg.cor)){
    t.inds=which(m6a.kegg.cor[i,]>0.6)
    c.kegg=colnames(m6a.kegg.cor)[t.inds]
    all.sig.r=rbind(all.sig.r,cbind(Gene=rep(row.names(m6a.kegg.cor)[i]
                                             ,length(c.kegg))
                                    ,KEGG=c.kegg,R=m6a.kegg.cor[i,t.inds]))
  }
  all.sig.m6a.kegg=rbind(all.sig.m6a.kegg,cbind(all.sig.r,rep(code,nrow(all.sig.r))))
}
lst.kegg=names(table(all.sig.m6a.kegg[,2]))
lst.kegg.rs=rbind()
g.c=c()
for(code in cancer.code){
  load(file = paste0(code,'.ssGSEA.kegg.RData'))
  row.names(ssGSEA.kegg)=gsub('KEGG_','',row.names(ssGSEA.kegg))
  exp=read.csv(file=paste0('exp/',code,'.tpm.txt'),sep = '\t',stringsAsFactors = F,check.names = F,row.names = 1)
  m6a.kegg.cor=cor(t(exp),t(ssGSEA.kegg))
  lst.kegg.rs=rbind(lst.kegg.rs,m6a.kegg.cor[,match(lst.kegg,colnames(m6a.kegg.cor))])
  ct=rep(code,20)
  names(ct)=row.names(m6a.kegg.cor)
  g.c=c(g.c,ct)
}
annotation_row = data.frame(
  GeneClass = g.c
)
row.names(lst.kegg.rs)=paste0(g.c,'_',row.names(lst.kegg.rs))
rownames(annotation_row) = row.names(lst.kegg.rs)
culst=pheatmap::pheatmap(t(lst.kegg.rs),cluster_cols = F
                         ,annotation_col = annotation_row,filename = 'Fig3B.pdf',width = 8,height = 6
                         ,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

matrix2vator=function(dat){
  vt=c()
  for(i in 1:nrow(dat)){
    vt=c(vt,as.numeric(dat[i,]))
  }
  return(vt)
}
table(cutree(culst$tree_row,3))
dim(lst.kegg.rs)
vd1.sbs=rbind()
for(i in seq(1,140,20)){
  cc=cancer.code[1+(i-1)/20]
  tmp=rbind()
  for( j in 1:3){
    t.inds=which(cutree(culst$tree_row,3)==j)
    pt=lst.kegg.rs[i:(i+19),t.inds]
    pt=matrix2vator(pt)
    tmp=rbind(tmp,cbind(Expression=pt,category=rep(paste0('C',j),length(pt))))
  }
  vd1.sbs=rbind(vd1.sbs,cbind(tmp,Type=rep(cc,nrow(tmp))))
}

colnames(vd1.sbs)=c('Expression','category','type')
vd1.sbs=as.data.frame(vd1.sbs)
vd1.sbs[,1]=as.numeric(as.character(vd1.sbs[,1]))

library(ggpubr)
p <- ggboxplot(vd1.sbs, x="type", y="Expression", color = "category", 
               palette = "jco"
               ,scales = "free_x"
               #, add = "jitter"
               #, facet.by = "fa"
               ,ncol=2, short.panel.labs = T)#按dose进
p=p+stat_compare_means(aes(group=category), label = "p.signif", method = "anova",
                       label.y = 1)
pdf('Fig3C',width = 6,height = 4)
p+theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()



############m6a cnv################
library("ComplexHeatmap")
library("circlize")

up=cbind()
down=cbind()
vd1.sbs=rbind()
cl.nm=c()
for(abbr in cancer.code){
  cnv=read.csv(paste0('cnv/',abbr,'_cnv_anno.txt'),sep = '\t',stringsAsFactors = F,header = F)
  muti=unique(unlist(strsplit(setdiff(names(table(cnv[,7])),m6a.genes),',')))
  cnv.all=read.csv(paste0('cnv/',abbr,'_CNV_T.txt'),sep = '\t',stringsAsFactors = F)
  m6a.exp.all=read.csv(paste0('exp/',abbr,'.tpm.txt'),sep = '\t',stringsAsFactors = F,row.names = 1,check.names = F)
  exp.samples=colnames(m6a.exp.all)
  print(c(abbr,length(unique(cnv[,4]))))
  cl.nm=c(cl.nm,abbr)
  cnal=length(unique(cnv.all[,1]))
  all.mt=rbind()
  
  dat.al=rbind()
  for(g in m6a.genes){
    if(sum(muti==g)){
      inds=grep(g,cnv[,7])
    }else{
      inds=which(cnv[,7]==g)
    }
    all.mt=rbind(all.mt,
                 c(sum(cnv[inds,6]>0.3)/cnal,
                   sum(cnv[inds,6]< -0.3)/cnal))
    amplification.samp=cnv[which(cnv[inds,6]>0),4]
    deletion.samp=cnv[which(cnv[inds,6]< 0),4]
    diploid.samp=setdiff(cnv.all[,1],c(amplification.samp,deletion.samp))
    
    amplification.samp=intersect(amplification.samp,exp.samples)
    deletion.samp=intersect(deletion.samp,exp.samples)
    diploid.samp=intersect(diploid.samp,exp.samples)
    
    ind1=match(deletion.samp,exp.samples)
    ind2=match(diploid.samp,exp.samples)
    ind3=match(amplification.samp,exp.samples)
    ind4=grep('11$',exp.samples)
    sp=c(rep('deletion',length(ind1)),rep('diploid',length(ind2))
         ,rep('amplification',length(ind3)),rep('Normal',length(ind4)))
    dat=data.frame(Expression=as.numeric(m6a.exp.all[which(row.names(m6a.exp.all)==g), c(ind1,ind2,ind3,ind4)])
                   ,category=sp,type=rep(g,length(sp)))
    dat.al=rbind(dat.al,dat)  
  }
  dat.al=cbind(dat.al,fa=rep(abbr,nrow(dat.al)))
  vd1.sbs=rbind(vd1.sbs,dat.al)
  up=cbind(up,all.mt[,1])
  down=cbind(down,all.mt[,2])
}

colnames(vd1.sbs)=c('Expression','category','type','fa')
vd1.sbs=as.data.frame(vd1.sbs)
vd1.sbs[,1]=log2(1+as.numeric(as.character(vd1.sbs[,1])))

library(ggpubr)
p <- ggboxplot(vd1.sbs[which(vd1.sbs[,2]!='diploid'),], x="type", y="Expression", color = "category", 
               palette = "jco",scales = "free_x"
               , facet.by = "fa",ncol=2, short.panel.labs = T)
p=p+stat_compare_means(aes(group=category), label = "p.signif", method = "anova",
                       label.y = 10)
pdf('FigS1.pdf',width = 10,height = 8)
p+theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

colnames(up)=cl.nm
colnames(down)=cl.nm
row.names(up)=m6a.genes
row.names(down)=m6a.genes
RightOrder <- rev(rownames(up))
UpColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFFFF","#00FF00"))
DnColor <- colorRamp2(breaks = c(0, 1), colors = c("#FFFFFF","#0000FF"))
row_an <-  HeatmapAnnotation(type = m61.genes.type, 
                             show_annotation_name = F, 
                             col = list(type = c("Writers" = "#5AC9FA", "Readers" = "#FAC67A", "Erasers" = "#51B743")), 
                             show_legend = T, 
                             annotation_legend_param = list(title = "m6A group",nrow = 1), 
                             which = "row" 
)

DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                 unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
    
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                 unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    if(up[i, j]>0.1){
      grid.text(label=round(up[i, j],2),x=(x + 0.5*width),
                y=(y+ 0.5*height),just = c('right','top'))
    }
    if(down[i, j]>0.1){
      
      grid.text(label=round(down[i, j],2),x=(x - 0.5*width),
                y=(y-0.5*height),just = c('left','bottom'))
    }
  }
}

p1 <- Heatmap(up, column_title = "Copy number variation across Tumor of digestive system"
              , rect_gp = gpar(type = "none")
              , show_heatmap_legend = F
              , cluster_rows = F
              , cluster_columns = T, 
              left_annotation = row_an, 
              cell_fun = DiagFunc(up = up, down = down) 
) 
col_fun = colorRamp2(c(-1, 0, 1), c("#0000FF", "#FFFFFF", "#00FF00")) 
lgd <- Legend(title = "CNV Frequency", 
              col_fun = col_fun,
              at = c(-1,-0.5,0,0.5,1), 
              labels = c("1","Loss","0","Gain","1"),  
              direction = "horizontal" 
)
pdf('Fig1C',width = 8,height = 6)
draw(p1, annotation_legend_list = lgd, annotation_legend_side = "bottom",heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()

############m6a mutation###################
library(maftools)
all.mt=cbind()
all.tmb=c()
all.samp.num=c()
for(abbr in cancer.code){
  data<-read.maf(maf=paste0('mutation/',abbr,'.maf'),isTCGA=T)
  gene2id=select(hgu133plus2.db,keys=data@gene.summary$Hugo_Symbol,keytype = 'SYMBOL',columns = c('ENTREZID'))
  all.samp.num=c(all.samp.num,as.numeric(data@summary[3,2]))
  mt.a=data@gene.summary$total[match(m6a.genes,gene2id$SYMBOL)]/as.numeric(data@summary[3,2])
  mt.a[is.na(mt.a)]=0
  all.mt=cbind(all.mt,mt.a)
  all.tmb=c(all.tmb,as.numeric(data@summary[14,2]))
}
row.names(all.mt)=m6a.genes
colnames(all.mt)=cancer.code
annotation_col = data.frame(
  MutationCount = all.tmb, 
  SampleNum = all.samp.num,
  MeanMutationCount = all.tmb*1.0/all.samp.num
)
row.names(annotation_col) = colnames(all.mt)
m61.genes.type=c(rep('Writers',7),rep('Readers',11),rep('Erasers',2))
annotation_row = data.frame(
  Regulaters = m61.genes.type
)
row.names(annotation_row) = row.names(all.mt)

pheatmap::pheatmap(all.mt,display_numbers = T,annotation_col = annotation_col,annotation_row = annotation_row,cluster_rows = F
                   ,filename = 'Fig1B.pdf',width = 8,height = 8)
################m6a methylation############################
m6a.cpg.anno=read.csv('methy/m6A.CpG.Anno.txt',sep = '\t',stringsAsFactors = F)
g.fcs=cbind()
g.ps=cbind()
for(code in cancer.code){
  methy=read.csv(paste0('methy/',code,'.450k.txt'),sep = '\t',row.names = 1,stringsAsFactors = F)
  n.inds=grep('11$',colnames(methy))
  t.inds=grep('01$',colnames(methy))
  print(c(code,length(t.inds),length(n.inds)))
  g.fc=c()
  g.p=c()
  for(g in unique(m6a.cpg.anno[,5])){
    g.mt=apply(methy[match(intersect(m6a.cpg.anno[which(m6a.cpg.anno[,5]==g),1],row.names(methy)),row.names(methy)),],2,function(x){
      x=x[!is.na(x)]
      return(mean(x))
    })
    g.p=c(g.p,t.test(g.mt[t.inds],g.mt[n.inds])$p.value)
    g.fc=c(g.fc,mean(g.mt[t.inds])/mean(g.mt[n.inds]))
  }
  g.fcs=cbind(g.fcs,g.fc)
  g.ps=cbind(g.ps,g.p)
}
colnames(g.ps)=cancer.code
colnames(g.fcs)=cancer.code
row.names(g.ps)=unique(m6a.cpg.anno[,5])
row.names(g.fcs)=unique(m6a.cpg.anno[,5])
up=-log10(g.ps)
down=log2(g.fcs)

up=up[match(m6a.genes[-8],row.names(g.ps)),]
down=down[match(m6a.genes[-8],row.names(g.ps)),]

UpColor <- colorRamp2(breaks = c(min(up), 5), colors = c("#FFFFFF","#6f9a8d"))
DnColor <- colorRamp2(breaks = c(-1,0, 1), colors = c("#0000FF",'#FFFFFF',"#FF0000"))

row_an <-  HeatmapAnnotation(type = m61.genes.type[-8], 
                             show_annotation_name = F, 
                             col = list(type = c("Writers" = "#5AC9FA", "Readers" = "#FAC67A", "Erasers" = "#51B743")),
                             show_legend = T,  
                             annotation_legend_param = list(title = "m6A group",nrow = 1), 
                             which = "row" 
)

DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                 unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
    
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                 unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    if(up[i, j]>=1.3){
      txt="****"
      if(up[i, j]>=1.3&up[i, j]<2){
        txt='*'
      }else if(up[i, j]>=2&up[i, j]<3){
        txt='**'
      }else if(up[i, j]>=3&up[i, j]<4){
        txt='***'
      }
      grid.text(label=txt,x=(x + 0.5*width),
                y=(y+ 0.5*height),just = c('right','top'))
    }
  }
}

p1 <- Heatmap(up, column_title = "Methylation across Tumor of digestive system"
              , rect_gp = gpar(type = "none")
              , show_heatmap_legend = F
              , cluster_rows = F
              , cluster_columns = T, 
              left_annotation = row_an, 
              cell_fun = DiagFunc(up = up, down = down) 
) 
col_fun = colorRamp2(c(-1, 0, 1), c("#0000FF", "#FFFFFF", "#FF0000")) 
lgd <- Legend(title = "log2(fold change)", 
              col_fun = col_fun, 
              at = c(-1,-0.5,0,0.5,1), 
              labels = c("-1","-0.5","0","0.5","1"),  
              direction = "horizontal" 
)
col_fun2 = colorRamp2(c(min(up), 5), c("#FFFFFF", "#6f9a8d")) 
lgd2 <- Legend(title = "-log10(p value)", 
               col_fun = col_fun2, 
               at = c(0,1,2,3,4,5), 
               labels = c('0',"1","2","3","4",">5"),  
               direction = "horizontal" 
)
pdf('Fig1D.pdf',width = 8,height = 6)
draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom",heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()
###############################Clinical##################################################
library("survival")
library("survminer") 
library("survivalROC")
coxFun=function(dat){
  colnames(dat)=c('time','status','gene')
  fmla <- as.formula("Surv(time, status) ~gene")
  cox <- coxph(fmla, data =dat)
  p <- summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
clinical.os=read.csv('PMC6066282-TCGA-CDR-clinical.txt',sep = '\t',stringsAsFactors = F)
head(clinical.os)

all.cancers.cox=rbind()
al.ps.Pfi=cbind()
al.hrs.Pfi=cbind()
al.ps.dfi=cbind()
al.hrs.dfi=cbind()
al.ps.dss=cbind()
al.hrs.dss=cbind()
al.ps.os=cbind()
al.hrs.os=cbind()

for(cd in cancer.code){
  tpm=read.csv(paste0('exp/',cd,'.tpm.txt')
               ,row.names = 1,stringsAsFactors = F,sep = '\t',check.names = F)
  tpm=tpm[,grep('-01$',colnames(tpm))]
  com.smp=intersect(colnames(tpm),paste0(clinical.os$bcr_patient_barcode,'-01'))
  tpm=tpm[,match(com.smp,colnames(tpm))]
  cl.os=clinical.os[match(com.smp,paste0(clinical.os$bcr_patient_barcode,'-01')),]
  cox.relt=t(apply(tpm, 1, function(x){
    return(c(coxFun(data.frame(cl.os$PFI.time,cl.os$PFI,as.numeric(x))),
             coxFun(data.frame(cl.os$DFI.time,cl.os$DFI,as.numeric(x))),
             coxFun(data.frame(cl.os$DSS.time,cl.os$DSS,as.numeric(x))),
             coxFun(data.frame(cl.os$OS.time,cl.os$OS,as.numeric(x)))
    ))
  }))
  colnames(cox.relt)=c(paste0('PFI_',c('p.value','HR','Low 95%CI','High 95%CI')),
                       paste0('DFI_',c('p.value','HR','Low 95%CI','High 95%CI')),
                       paste0('DSS_',c('p.value','HR','Low 95%CI','High 95%CI')),
                       paste0('OS_',c('p.value','HR','Low 95%CI','High 95%CI')))
  all.cancers.cox=rbind(all.cancers.cox,cbind(cox.relt,Cancer=rep(cd,20)))
  al.ps.Pfi=cbind(al.ps.Pfi,cox.relt[,1])
  al.hrs.Pfi=cbind(al.hrs.Pfi,cox.relt[,2])
  al.ps.dfi=cbind(al.ps.dfi,cox.relt[,5])
  al.hrs.dfi=cbind(al.hrs.dfi,cox.relt[,6])
  al.ps.dss=cbind(al.ps.dss,cox.relt[,9])
  al.hrs.dss=cbind(al.hrs.dss,cox.relt[,10])
  al.ps.os=cbind(al.ps.os,cox.relt[,13])
  al.hrs.os=cbind(al.hrs.os,cox.relt[,14])
  
}
colnames(al.ps.Pfi)=cancer.code
colnames(al.hrs.Pfi)=cancer.code
colnames(al.ps.dfi)=cancer.code
colnames(al.hrs.dfi)=cancer.code
colnames(al.ps.dss)=cancer.code
colnames(al.hrs.dss)=cancer.code
colnames(al.ps.os)=cancer.code
colnames(al.hrs.os)=cancer.code

write.table(cbind(Gene=row.names(all.cancers.cox),all.cancers.cox),file='all.m6a.clinical.cox.txt',sep = '\t',quote = F,row.names = F)

up_break=c(0, 5)
down_break=c(-0.4,0, 0.1)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("#0000FF",'#FFFFFF',"#FF0000")

up=-log10(al.ps.Pfi)
down=log2(al.hrs.Pfi)
title='PFI across Tumor of digestive system types'
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)

up=-log10(al.ps.dfi)
down=log2(al.hrs.dfi)
title='DFI across Tumor of digestive system types'
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)

up=-log10(al.ps.dss)
down=log2(al.hrs.dss)
title='DSS across Tumor of digestive system types'
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)

up=-log10(al.ps.os)
down=log2(al.hrs.os)
title='OS across Tumor of digestive system types'
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)

up=-log10(cbind(al.ps.Pfi,al.ps.dfi,al.ps.dss,al.ps.os))
down=log2(cbind(al.hrs.Pfi,al.hrs.dfi,al.hrs.dss,al.hrs.os))
title='Clinical across Tumor of digestive system types'
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)

plotMutiHeatmap=function(up,down,up_break,up_colors,down_break,down_colors,title){
  UpColor <- colorRamp2(breaks = up_break, colors = up_colors)
  DnColor <- colorRamp2(breaks = down_break, colors = down_colors)
  
  row_an <-  HeatmapAnnotation(type = m61.genes.type, 
                               show_annotation_name = F, 
                               col = list(type = c("Writers" = "#5AC9FA", "Readers" = "#FAC67A", "Erasers" = "#51B743")), 
                               show_legend = T,  
                               annotation_legend_param = list(title = "m6A group",nrow = 1), 
                               which = "row" 
  )
  col_an <-  HeatmapAnnotation(type = c(rep('PFI',7),rep('DFI',7),rep('DSS',7),rep('OS',7)), 
                               show_annotation_name = F, 
                               col = list(type = c("PFI" = "#ffe3ed", "DFI" = "#dcffcc", "DSS" = "#beebe9",'OS'="#9be3de")), 
                               show_legend = T,  
                               annotation_legend_param = list(title = "Clinical group",nrow = 1), 
                               which = "column" 
  )
  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                   unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                   unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
      if(up[i, j]>=1.3){
        txt="***"
        if(up[i, j]>=1.3&up[i, j]<2){
          txt='*'
        }else if(up[i, j]>=2&up[i, j]<3){
          txt='**'
        }else if(up[i, j]>=3&up[i, j]<4){
          txt='***'
        }
        grid.text(label=txt,x=(x + 0.5*width),
                  y=(y+ 0.5*height),just = c('right','top'))
      }
      if(down[i, j]>0.1){
      }
    }
  }
  
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = F
                , cluster_columns = F, 
                left_annotation = row_an, 
                top_annotation = col_an,
                cell_fun = DiagFunc(up = up, down = down) 
  ) 
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "log2(HR)", 
                col_fun = col_fun, 
                at = c(-0.4,-0.2,0,0.1), 
                labels = c("-0.4","-0.2","0","0.1"),  
                direction = "horizontal" 
  )
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "-log10(p value)", 
                 col_fun = col_fun2, 
                 at = c(0,1,2,3,4,5), 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "horizontal"
  )
  
  draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
       ,heatmap_legend_side = "bottom", merge_legend = TRUE)
}

#########
library('hgu133plus2.db')

library('GSVA')
library(GSEABase)
n.m6a.genes=m6a.genes
n.m6a.genes[4]='VIRMA'

n.m6a.genes=n.m6a.genes[which(apply(cbind(al.ps.Pfi,al.ps.dfi,al.ps.dss,al.ps.os),1,function(x){
  return(sum(as.numeric(x)<0.05))
})>3)]

gs=GeneSet(setName='m6a_score', setIdentifier="101",geneIds=n.m6a.genes,SymbolIdentifier()) 
gsc <- GeneSetCollection(list(gs))
fl <- tempfile()
toGmt(gsc, fl)
c2m6a=getGmt(fl)

all.cancer.m6a.score=rbind()
for(code in cancer.code){
  m.exp=read.csv(paste0('exp/',code,'.symbol.tpm.txt'),sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
  ssGSEA.m6a <- gsva(as.matrix(m.exp), c2m6a,method='ssgsea',min.sz=1, max.sz=500, verbose=TRUE)
  all.cancer.m6a.score=rbind(all.cancer.m6a.score,cbind(Sample=colnames(ssGSEA.m6a),m6a_score=t(ssGSEA.m6a),Cancer=rep(code,length(ncol(ssGSEA.m6a)))))
}
clinical.os=read.csv('PMC6066282-TCGA-CDR-clinical.txt',sep = '\t',stringsAsFactors = F)
com.smp=intersect(all.cancer.m6a.score[,1],paste0(clinical.os$bcr_patient_barcode,'-01'))
all.cancer.m6a.score=all.cancer.m6a.score[match(com.smp,all.cancer.m6a.score[,1]),]
clinical.os.m6a=clinical.os[match(com.smp,paste0(clinical.os$bcr_patient_barcode,'-01')),]
plotKMCox=function(dat){
  colnames(dat)=c('time','status','groups')
  sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  print((sdf))
  p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  colKm=rainbow(length(sf$strata))
  plot(sf, mark.time = TRUE,col=colKm,xlab=paste("Survival time in day","\np=",round(p,5)),ylab = "Survival probabilities",main="Method Kaplan Meier")
  legend('topright',paste0(gsub('groups=','',names(sf$strata)),'(N=',sdf$n,')'), col = colKm,
         lty = c(1,1, 1, 1),lwd=c(1,1,1,1),merge = TRUE,cex = 0.8)
  return(p)
}
par(mfrow=c(2,4))
all.m6a.score.cox=rbind()
for(code in cancer.code){
  inds=which(all.cancer.m6a.score[,3]==code)
  x=as.numeric(all.cancer.m6a.score[inds,2])
  all.m6a.score.cox=rbind(
    all.m6a.score.cox,
    c(coxFun(data.frame(clinical.os.m6a$PFI.time[inds],clinical.os.m6a$PFI[inds],as.numeric(x))),
      coxFun(data.frame(clinical.os.m6a$DFI.time[inds],clinical.os.m6a$DFI[inds],as.numeric(x))),
      coxFun(data.frame(clinical.os.m6a$DSS.time[inds],clinical.os.m6a$DSS[inds],as.numeric(x))),
      coxFun(data.frame(clinical.os.m6a$OS.time[inds],clinical.os.m6a$OS[inds],as.numeric(x)))
    ))
  gp=ifelse(x>median(x),'H','L')
  plotKMCox(data.frame(clinical.os.m6a$PFI.time[inds],clinical.os.m6a$PFI[inds],gp))
}
row.names(all.m6a.score.cox)=cancer.code
cox.tb=cbind(round(all.m6a.score.cox[,2],2),paste0(round(all.m6a.score.cox[,3],2),'-',round(all.m6a.score.cox[,4],2))
             ,round(all.m6a.score.cox[,1],3),round(all.m6a.score.cox[,6],2),paste0(round(all.m6a.score.cox[,7],2),'-',round(all.m6a.score.cox[,8],2))
             ,round(all.m6a.score.cox[,5],3),round(all.m6a.score.cox[,10],2),paste0(round(all.m6a.score.cox[,11],2),'-',round(all.m6a.score.cox[,12],2))
             ,round(all.m6a.score.cox[,9],3),round(all.m6a.score.cox[,14],2),paste0(round(all.m6a.score.cox[,15],2),'-',round(all.m6a.score.cox[,16],2))
             ,round(all.m6a.score.cox[,13],3))
cox.tb=rbind(rep(c('HR','95% CI','P'),4),cox.tb)
write.table(cbind(row.names(cox.tb),cox.tb)
            ,file = 'm65.gsva.score.cox.txt',sep = '\t',quote = F,row.names = F)

code='PAAD'
inds=which(all.cancer.m6a.score[,3]==code)
x=as.numeric(all.cancer.m6a.score[inds,2])
rocplot(clinical.os.m6a$PFI.time[inds],clinical.os.m6a$PFI[inds],x,c(12*30,12*3*30,12*5*30))

test.score=scale(x)[,1]

fit <- survivalROC(Stime =clinical.os.m6a$OS.time[inds], status = clinical.os.m6a$OS[inds]
                   , marker = test.score, predict.time = 365*2, method = "KM")
cutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]

dat=data.frame(clinical.os.m6a$OS.time[inds],clinical.os.m6a$OS[inds],ifelse(test.score>cutoff,'H','L'))
colnames(dat)=c('Days','vital_status','group')
fit_validation <- survfit( Surv(Days, vital_status) ~ group,data = dat )

#BiocManager::install('survminer')
library(survminer)
ggsurvplot(fit_validation, data = dat,
           conf.int = TRUE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("red","#2E9FDF"),
           legend = "bottom",
           legend.title = "RiskGroup",
           legend.labs = c("High","Low"))

rocplot(clinical.os.m6a$PFI.time[inds],clinical.os.m6a$PFI[inds],test.score,c(12*30,12*3*30,12*5*30))

#######
m.exp=read.csv(paste0('exp/','PAAD','.tpm.txt')
               ,sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
code='PAAD'
inds=which(all.cancer.m6a.score[,3]==code)

m.exp=m.exp[,match(all.cancer.m6a.score[inds,1],colnames(m.exp))]
dim(m.exp)
m6a.genes
n1.m6a.genes=n.m6a.genes
n1.m6a.genes[4]='KIAA1429'
lst.exp=m.exp[match(n1.m6a.genes,row.names(m.exp)),]
cbind()
least.cpgs=row.names(lst.exp)
cox.dna.tcga.exp=cbind(time=clinical.os.m6a$OS.time[inds],status=clinical.os.m6a$OS[inds],t(lst.exp))
colnames(cox.dna.tcga.exp)=c('time', 'status',least.cpgs)
fmla <- as.formula(paste0("Surv(time, status) ~",paste0(least.cpgs,collapse = '+')))
cox <- coxph(fmla, data = as.data.frame((cox.dna.tcga.exp)))
summary(cox)
scox=step(cox)
cox.score=t(lst.exp[match(names(scox$coefficients),row.names(lst.exp)),])%*%coefficients(scox)
rocplot(clinical.os.m6a$PFI.time[inds],clinical.os.m6a$PFI[inds],cox.score,c(12*30,12*3*30,12*5*30))

bk = unique(c(seq(-3,3, length=50)))
cutoff
annotation_col = data.frame(
  m6a_score = as.numeric(all.cancer.m6a.score[inds,2])
)
rownames(annotation_col) = all.cancer.m6a.score[inds,1]
pheatmap::pheatmap(log2(1+m.exp[,order(as.numeric(all.cancer.m6a.score[inds,2]))])
                   ,color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
                   ,annotation_col = annotation_col
                   ,cluster_cols = F,scale = 'row',breaks = bk,show_colnames = F)

all.exp=read.csv(paste0('exp/','PAAD','.symbol.tpm.txt')
                 ,sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
inds=which(all.cancer.m6a.score[,3]=='PAAD')

emt=c('TGFBR1','SNAI2','ITGA5','FZD1','EDNRA','TGFBR2','TWIST1','SNAI1','WNT4','ITGB1')
emt.exp=all.exp[match(emt,row.names(all.exp)),]
emt.exp=emt.exp[,match(all.cancer.m6a.score[inds,1],colnames(emt.exp))]

pheatmap::pheatmap(log2(1+emt.exp[,order(as.numeric(all.cancer.m6a.score[inds,2]))]),cluster_cols = F
                   ,scale = 'row',breaks = bk,show_colnames = F,annotation_col = annotation_col
                   ,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pheatmap::pheatmap(log2(1+nos.target.exp[,order(as.numeric(all.cancer.m6a.score[inds,2]))]),cluster_cols = F
                   ,scale = 'row',breaks = bk,show_colnames = F,annotation_col = annotation_col
                   ,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

group=ifelse(test.score>cutoff,'H','L')

vd1.sbs=rbind()
for(i in 1:nrow(emt.exp)){
  vd1.sbs=rbind(vd1.sbs,cbind(Expression=as.numeric(emt.exp[i,]),category=group,type=rep(row.names(emt.exp)[i],length(group))
                              ,type=rep('EMT-related genes',length(group))))
}
colnames(vd1.sbs)=c('Expression','category','type','fa')
vd1.sbs=as.data.frame(vd1.sbs)
vd1.sbs[,1]=log2(1+as.numeric(as.character(vd1.sbs[,1])))
library(ggpubr)
p <- ggboxplot(vd1.sbs, x="type", y="Expression", color = "category", 
               palette = "jco",scales = "free_x"
               , facet.by = "fa",ncol=2, short.panel.labs = T)
p=p+stat_compare_means(aes(group=category), label = "p.signif", method = "anova",
                       label.y = 15)
p+theme(axis.text.x = element_text(angle = 30, hjust = 1))
####vd#####
vd.exp=read.csv('ICGC_data/ICGC_PAAD_exp.txt'
                ,sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
vd.cli=read.csv('ICGC_data/PACA-AU.tsv'
                ,sep = '\t',row.names = 1,stringsAsFactors = F,check.names = F)
vd.cli=vd.cli[match(colnames(vd.exp),row.names(vd.cli)),]
vd.m6a.genes=n.m6a.genes
vd.m6a.genes[4]='KIAA1429'
match(vd.m6a.genes,row.names(vd.exp))
gs=GeneSet(setName='m6a_score', setIdentifier="101",geneIds=vd.m6a.genes,SymbolIdentifier()) 
gsc <- GeneSetCollection(list(gs))
fl <- tempfile()
toGmt(gsc, fl)
c2m6a=getGmt(fl)
vd.ssGSEA.m6a <- gsva(as.matrix(vd.exp), c2m6a,method='ssgsea',min.sz=1, max.sz=500, verbose=TRUE)
annotation_col = data.frame(
  m6a_score = as.numeric(vd.ssGSEA.m6a[1,])
)
rownames(annotation_col) = colnames(vd.exp)
pheatmap::pheatmap(vd.exp[match(m6a.genes,row.names(vd.exp))
                          ,order(as.numeric(vd.ssGSEA.m6a[1,]))]
                   ,color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
                   ,annotation_col = annotation_col
                   ,cluster_cols = F,scale = 'row',breaks = bk,show_colnames = F)

plotKMCox(data.frame(vd.cli$donor_survival_time
                     ,ifelse(vd.cli$donor_vital_status=='deceased',1,0),cutree(clust$tree_col,2)))

vd.dat=data.frame(vd.cli$donor_survival_time
                  ,ifelse(vd.cli$donor_vital_status=='deceased',1,0)
                  ,as.numeric(vd.ssGSEA.m6a[1,]))
vd.dat=vd.dat[!is.na(vd.dat[,1]),]
vd.dat=vd.dat[as.numeric(vd.dat[,1])>30,]

vd.score=scale(vd.dat[,3])[,1]
coxFun(vd.dat)
plotKMCox(data.frame(vd.dat[,1],vd.dat[,2],ifelse(vd.score>-0.6099,'H','L')))

dat=data.frame(vd.dat[,1],vd.dat[,2],ifelse(vd.score>cutoff,'H','L'))
colnames(dat)=c('Days','vital_status','group')
fit_validation <- survfit( Surv(Days, vital_status) ~ group,data = dat )


#BiocManager::install('survminer')
library(survminer)
ggsurvplot(fit_validation, data = dat,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("red","#2E9FDF"),
           legend = "bottom",
           legend.title = "RiskGroup",
           legend.labs = c("High","Low"))

rocplot(vd.dat[,1],vd.dat[,2],vd.dat[,3],c(365,365*3,365*5))

########upset########
library(psych)
all.m6a_score.cor=rbind()
for(code in cancer.code){
  inds=which(all.cancer.m6a.score[,3]==code)
  m6a_score=all.cancer.m6a.score[inds,2]
  load(paste0('exp/',code,'.ssGSEA.kegg.RData'))
  ssGSEA.kegg=ssGSEA.kegg[,match(names(m6a_score),colnames(ssGSEA.kegg))]
  ct=corr.test(as.numeric(m6a_score),t(ssGSEA.kegg))
  m6a.cor=cbind(KEGG=row.names(ssGSEA.kegg),R=ct$r[1,],P=ct$p[1,],FDR=p.adjust(ct$p[1,]),rep(code,nrow(ssGSEA.kegg)))
  all.m6a_score.cor=rbind(all.m6a_score.cor,m6a.cor)
}
sig.cor=all.m6a_score.cor[abs(as.numeric(all.m6a_score.cor[,2]))>0.4&as.numeric(all.m6a_score.cor[,4])<0.01,]
write.table(sig.cor,file = 'all.m6a_score.cor.txt',sep = '\t',quote = F,row.names = F)
library('UpSetR')
library(UpSetR)
up_set_list=list()
for(code in cancer.code){
  up_set_list=append(up_set_list,list(as.character(sig.cor[sig.cor[,5]==code,1]))) 
}
colnames(up_set_list)
up_input=fromList(up_set_list)
colnames(up_input)=cancer.code
upset(up_input, order.by = "freq",nsets=8)

last.kegg=names(which(table(sig.cor[,1])>5))

sig.cor[sig.cor[,1]%in%names(which(table(sig.cor[,1])>5)),]

par(mfrow=c(4,7))
for(code in cancer.code){
  inds=which(all.cancer.m6a.score[,3]==code)
  m6a_score=all.cancer.m6a.score[inds,2]
  load(paste0('exp/',code,'.ssGSEA.kegg.RData'))
  ssGSEA.kegg=ssGSEA.kegg[match(last.kegg,row.names(ssGSEA.kegg)),match(names(m6a_score),colnames(ssGSEA.kegg))]
  ct=corr.test(as.numeric(m6a_score),t(ssGSEA.kegg))
  for(i in 1:nrow(ssGSEA.kegg)){
    plot(as.numeric(ssGSEA.kegg[i,]),as.numeric(m6a_score)
         ,pch=20,col='blue',ylab='m6a score'
         ,xlab=gsub('KEGG_','',row.names(ssGSEA.kegg)[i])
    )
    legend("bottomright", c(paste0('R=',round(ct$r[i],2))
                            ,'FDR<0.001')
           ,cex = 0.9)
  }
}
par(mfrow=c(1,7))
for(code in cancer.code){
  inds=which(all.cancer.m6a.score[,3]==code)
  m6a_score=all.cancer.m6a.score[inds,2]
  plot(hist(as.numeric(m6a_score)))
}

#########################################immue#############################
library(estimate)

gs=GeneSet(setName='m6a_score', setIdentifier="101",geneIds=n.m6a.genes,SymbolIdentifier()) 
gsc <- GeneSetCollection(list(gs))
toGmt(gsc, fl)
c2m6a=getGmt(fl)

all.est.m6a.score=rbind()
for(code in cancer.code){
  f1 <- tempfile()
  f2 <- tempfile()
  
  outputGCT(paste0('exp/',code,'.symbol.tpm.txt')
            , f1)
  estimateScore(f1,f2,platform='illumina')
  est.score=t(read.csv(f2,sep = '\t',row.names = 1,check.names = F,skip = 2))
  est.score=est.score[-1,]
  est.score=cbind(est.score,Cancer=rep(code,nrow(est.score)))
  rnames=gsub('\\.','-',row.names(est.score))
  #est.score=apply(est.score, 2, as.numeric)
  row.names(est.score)=rnames
  all.est.m6a.score=rbind(all.est.m6a.score,est.score)
}
dim(est.score)
write.table(cbind(Sample=row.names(all.est.m6a.score),all.est.m6a.score)
            ,file='all.est.score.txt'
            ,sep = '\t',quote = F,row.names = F)
all.est.m6a.score=all.est.m6a.score[match(all.cancer.m6a.score[,1],row.names(all.est.m6a.score)),]

esm.all.rs=cbind()
esm.all.ps=cbind()

for(code in cancer.code){
  inds=which(all.cancer.m6a.score[,3]==code)
  m6a_score=all.cancer.m6a.score[inds,2]
  immu_score1=as.numeric(all.est.m6a.score[inds,1])
  immu_score2=as.numeric(all.est.m6a.score[inds,2])
  immu_score3=as.numeric(all.est.m6a.score[inds,3])
  ct1=corr.test(as.numeric(m6a_score),immu_score1)
  ct2=corr.test(as.numeric(m6a_score),immu_score2)
  ct3=corr.test(as.numeric(m6a_score),immu_score3)
  esm.all.rs=cbind(esm.all.rs,c(ct1$r,ct2$r,ct3$r))  
  esm.all.ps=cbind(esm.all.ps,c(ct1$p,ct2$p,ct3$p))  
}
colnames(esm.all.rs)=cancer.code
colnames(esm.all.ps)=cancer.code
row.names(esm.all.rs)=colnames(all.est.m6a.score)[1:3]
row.names(esm.all.ps)=colnames(all.est.m6a.score)[1:3]


corrplot(esm.all.rs
         #,p.mat = res1$p
         , insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white"
         #, type = "upper"
         ,col=col3(100),tl.cex=0.8)

up=-log10(esm.all.ps)
down=esm.all.rs

up_break=c(0,5)
down_break=c(-0.5,0.5)
up_colors=c('#fff7f3','#fcc5c0')
down_colors=c('#addd8e','#ffffff')

title='Immue-m6a_score across Tumor of digestive system types'
UpColor <- colorRamp2(breaks = up_break, colors = up_colors)
DnColor <- colorRamp2(breaks = down_break, colors = down_colors)

DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                 unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
    
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                 unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
    if(up[i, j]>=1.3){
      txt="***"
      if(up[i, j]>=1.3&up[i, j]<2){
        txt='*'
      }else if(up[i, j]>=2&up[i, j]<3){
        txt='**'
      }else if(up[i, j]>=3&up[i, j]<4){
        txt='***'
      }
      grid.text(label=txt,x=(x + 0.5*width),
                y=(y+ 0.5*height),just = c('right','top'))
    }
    if(down[i, j]>0.1){
    }
  }
}

p1 <- Heatmap(up, column_title = title
              , rect_gp = gpar(type = "none")
              , show_heatmap_legend = F
              , cluster_rows = F
              , cluster_columns = F, 
              cell_fun = DiagFunc(up = up, down = down) 
) 
col_fun = colorRamp2(down_break, down_colors) 
lgd <- Legend(title = "log2(HR)", 
              col_fun = col_fun, 
              at = c(-0.4,-0.2,0,0.1), 
              labels = c("-0.4","-0.2","0","0.1"),  
              direction = "horizontal" 
)
col_fun2 = colorRamp2(up_break, up_colors) 
lgd2 <- Legend(title = "-log10(p value)", 
               col_fun = col_fun2, 
               at = c(0,1,2,3,4,5), 
               labels = c('0',"1","2","3","4",">5"),  
               direction = "horizontal" 
)

draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
     ,heatmap_legend_side = "bottom", merge_legend = TRUE)

########The end#########################################



