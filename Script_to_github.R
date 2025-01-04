library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(ggstream)
library(ggridges)
library(gghalves)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(ggrastr)
library(monocle)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(stringr)
library(patchwork)
library(scales)
library(data.table)
library(reshape)
library(reshape2)

##01 cell identification ----------------------------------------------------
seurat_rds <- readRDS(paste0(ana_dir,'integrated_sample.rds'))
DefaultAssay(seurat_rds) <- "integrated"
seurat_rds <- RunPCA(seurat_rds, verbose = FALSE)

pc.list <- pc.num <- c(1:18)
n.neighbors=10
min.dist=0.3
spread=1.1
resolution=0.8
seurat_rds <- FindNeighbors(seurat_rds, reduction = "pca",dims = pc.list,  assay= "integrated")
seurat_rds <- FindClusters(object = seurat_rds, resolution = resolution)
seurat_rds <- RunUMAP(seurat_rds, reduction = "pca",dims = pc.list,  assay= "integrated",
                       n.neighbors= n.neighbors,
                       min.dist= min.dist,
                       spread= spread
)


p4 <- DimPlot(seurat_rds, reduction = "umap", label = TRUE,pt.size = 0.1,
              group.by = 'seurat_clusters',raster=F)+
  theme(plot.background =element_rect(fill = NA,colour = NA),
        #panel.border = element_rect(fill = NA,colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_rect(fill = NA,colour = NA),
        #,
        plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold',color = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.line = element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        #legend.position = "none",
        plot.margin=unit(c(0,0,0,0),'lines')#up right down left
  )
ggsave(paste0(ana_dir,"clt_umap",".png"), plot = p4, width = 10, height =8)
ggsave(paste0(ana_dir,"clt_umap",".pdf"), plot = p4, width = 10, height =8)

p4 <- DimPlot(seurat_rds, reduction = "umap", label = TRUE,pt.size = 0.1,
              group.by = 'cell_name',raster=F)+
  theme(plot.background =element_rect(fill = NA,colour = NA),
        #panel.border = element_rect(fill = NA,colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_rect(fill = NA,colour = NA),
        #,
        plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold',color = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.line = element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        #legend.position = "none",
        plot.margin=unit(c(0,0,0,0),'lines')#up right down left
  )
ggsave(paste0(ana_dir,"cell_umap",".pdf"), plot = p4, width = 10, height =8)

p4 <- DimPlot(seurat_rds, reduction = "umap", label = TRUE,pt.size = 1,raster=T,raster.dpi = c(300, 300))+
  theme(plot.background =element_rect(fill = NA,colour = NA),
        #panel.border = element_rect(fill = NA,colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_rect(fill = NA,colour = NA),
        #,
        plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold',color = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.line = element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        #legend.position = "none",
        plot.margin=unit(c(0,0,0,0),'lines')#up right down left
  )
ggsave(paste0(ana_dir,"umap_rst",".pdf"), plot = p4, width = 10, height =8)

seurat_rds@meta.data$disease <- factor(seurat_rds@meta.data$disease,levels = c('Healthy','Disease'),ordered = T)
p4 <- DimPlot(seurat_rds, reduction = "umap",group.by = 'disease', 
              split.by = 'sample',ncol = 3,cols =  c(Healthy='#FF6347',Disease='#48D1CC'),
              
              label = F,pt.size = 1,raster=F)+
  theme(plot.background =element_rect(fill = NA,colour = NA),
        #panel.border = element_rect(fill = NA,colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_rect(fill = NA,colour = NA),
        #,
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold',color = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.line = element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0),'lines')#up right down left
  )
ggsave(paste0(ana_dir,ts,"cl_umap_sp",".png"), plot = p4, width = 14, height =8)


##01.2 cellmark----------------------------------------------------
Idents(seurat_rds) <- seurat_rds@meta.data$cell_name
DefaultAssay(seurat_rds) <- "RNA"
df <- FindAllMarkers(seurat_rds)
df$gene <- rownames(df)
df <- df[df$avg_log2FC>0.5&df$p_val_adj<=0.05,]
df <- df%>%arrange(cluster,desc(avg_log2FC))
write.csv(df,paste0(ana_dir,'cellmrk.csv'))

df <- read.csv(paste0(ana_dir,'cellmrk.csv'))
df$cluster <- factor(df$cluster,levels = c("Epi","EC","Fib",
                                           "Mye","TC","BC"))
df2 <- df[df$avg_log2FC>0.5&df$p_val_adj<=0.05,]
write.csv(df2,paste0(ana_dir,'marker_gn_tt.csv' ),row.names = F)
cellmarkers <- df2 
cellmarkers$gene <- apply(cellmarkers,1,function(vec){strsplit(vec['gene'],'\\.')[[1]][1] })

mark <- cellmarkers[cellmarkers$p_val_adj<0.01&cellmarkers$avg_log2FC>0.5,]
mark <- mark[order(mark$avg_log2FC,decreasing = T),]
gene_ls <- c()
for (cell in levels(mark$cluster)) {
  tmp <- mark[mark$cluster==cell,]
  gene_tmp <- tmp$gene[1:50]
  gene_ls <- c(gene_ls,gene_tmp)
}
write.csv(gene_ls,paste0(ana_dir,"top_mkr.csv"),row.names = F)
DefaultAssay(seurat_rds) <- "RNA"
cluster.averages <- AverageExpression(seurat_rds, return.seurat = TRUE)
gene_ls <- gene_ls[gene_ls%in%rownames(cluster.averages@assays$RNA@features)]
mtx <- cluster.averages@assays$RNA@layers$scale.data[,]
rownames(mtx) <- rownames(cluster.averages@assays$RNA@features)
mtx <- mtx[gene_ls,]
hmcols <- c(colorRampPalette(c('#6495ED', 'white'))(10),
            colorRampPalette(c('white', '#EE3B3B'))(10)
)
pheatmap::pheatmap(mtx,
                   cluster_rows = F,
                   cluster_cols =F,
                   scale = "row",#"none", "row", "column"
                   color = hmcols,
                   show_rownames = F,
                   border_color = NA, #"grey60",
                   fontsize = 14,
                   fontsize_row = 10, #fontsize,
                   fontsize_col = 20, #fontsize,
                   fontsize_number = 0.8* fontsize,
                   kmeans_k = NA,
                   cutree_rows = NA,
                   treeheight_row=0,
                   cutree_cols = NA,
                   legend_breaks = NA,
                   legend_labels = NA,
                   annotation = NA,
                   annotation_col= NA,
                   annotation_legend = TRUE,
                   drop_levels = TRUE,
                   #gaps_row = seq(1,nrow(heatmap_up)),
                   gaps_col = c(1:9),
                   gaps_row = c(seq(0,450,50)),
                   labels_row = NULL,
                   labels_col = NULL,
                   cellwidth = 30,
                   cellheight = 2,
                   revC=F,
                   filename=paste0(ana_dir,'mrk_htmp.pdf')
)

DotPlot<-DotPlot(
  seurat_rds,
  cols = c("lightgrey", "#D52126"),
  features= unique(c('AOC1','APOC4','ALPI','PECAM1','KDR','CDH5',
                     'LUM','DCN','PDGFRA','CD68','CD14','LYZ',
                     'CD3D','CD3E','CD3G','CD19','CD79A','CD79B')),
  dot.scale=10,
  col.min = 0)+
  ggplot2::theme(
    
    panel.border = element_rect(fill = 'transparent',colour ='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major=element_blank(),
    panel.background = element_rect(fill = 'transparent',colour = NA),
    plot.background = element_rect(fill = 'transparent',colour = NA),
    axis.text.x = element_text(size=15,angle=90,vjust = 0.5,hjust=1),
    axis.text.y = element_text(size=20,angle=0,vjust = 0.5,hjust=0),
    axis.title.x = element_blank(),
    axis.title.y=element_blank()
  )
ggsave(paste0(tmp_dir,"_cell_mrk_dot.pdf"), plot = DotPlot, width = 20, height = 4)
ggsave(paste0(tmp_dir,"_cell_mrk_dot.png"), plot = DotPlot, width = 20 , height = 4)
#01.3 cellprop------------------------------------------
meta <- df
meta$ag_sp <- paste0(meta$spp,'_',meta$sample)
cell_number=table(meta$cell_new,meta$spp)
cell.prop <- data.frame(cell_number)
colnames(cell.prop) <- c('celltype', 'disease', 'number')

prop <- cell.prop%>%dplyr::group_by(disease)%>%dplyr::mutate(num_tt=sum(number))
prop$prop <- prop$number/prop$num_tt
prop$disease <- factor(prop$disease,levels = c('Normal','Tumor'))
prop$celltype<- factor(prop$celltype,levels = c('Epithelial cells',
                                                'Endothelial cells',
                                                'Fibroblast cells',
                                                'Myeloid cells',
                                                'T cells',
                                                'B cells'))
prop$sp_cl <- paste0(prop$sample,'_',prop$celltype)


ggplot(prop,aes(x=celltype,y=prop,fill=disease),color='white')+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_fill_manual(values = c('#af4212','#1e2973'),labels=c("healthy","Tumor"))+
  theme(panel.border = element_rect(fill = 'transparent',colour = 'black'),
        panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_blank(),
        plot.background =element_blank()#,
        #plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold'),
        # axis.ticks = element_blank(),
        # axis.text = element_blank(),
        # axis.title.x =element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        # legend.position = "none"
  )
ggsave(paste0(workwd,'cell_prop.pdf'),width = 12,height = 6)


prop2 <- cast(prop,celltype~disease,value = 'prop')
prop2$FC <- prop2$Tumor/prop2$Normal
prop2$logFC <- log2(prop2$FC)
prop2$drc <- ifelse(prop2$logFC>0,'U','D')

ggplot(prop2,aes(x=celltype,y=logFC,fill=drc),color='white')+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('#33b2e8','#da2d20'),labels=c("Increased","Decreased"))+
  theme(panel.border = element_rect(fill = 'transparent',colour = 'black'),
        panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_blank(),
        plot.background =element_blank()#,
        #plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold'),
        # axis.ticks = element_blank(),
        # axis.text = element_blank(),
        # axis.title.x =element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        # legend.position = "none"
  )
ggsave(paste0(workwd,'cell_prop_FC.pdf'),width = 12,height = 6)


#02 DEG analysis--------------------------
seurat_rds@meta.data$cl_dse <- paste0(seurat_rds@meta.data$cell_name,'_',seurat_rds@meta.data$disease)
Idents(seurat_rds) <- seurat_rds@meta.data$cl_dse

#02.1 DEG calculate--------------------------
S2_S1_tmp <- NULL
for (cell in unique(seurat_rds@meta.data$cell_name) ) {
  tmp <- FindMarkers(seurat_rds, ident.1=paste0(cell, '_Disease'), ident.2=paste0(cell, '_Healthy'),logfc.threshold = 0.2)
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$group <- paste0('Disease_Healthy')
  S2_S1_tmp <- rbind(S2_S1_tmp, tmp)
  print(paste0(cell, ' is finished'))
}
S2_S1_tmp <- S2_S1_tmp[S2_S1_tmp$p_val_adj<0.05&abs(S2_S1_tmp$avg_log2FC)>0.25,]
write.csv(S2_S1_tmp,paste0(DEG_dir,'Dse_Hty.csv'),row.names = F)

df <- read.csv(paste0(DEG_dir,'Dse_Hty.csv'))
f$direction <- ifelse(df$avg_log2FC>0,'U','D')


df$celltype <- factor(df$celltype,levels = c('Epithelial cells',
                                             'Endothelial cells',
                                             'Fibroblast cells',
                                             'Myeloid cells',
                                             'T cells',
                                             'B cells'))
df <- df%>%group_by(celltype)%>%arrange(celltype,avg_log2FC)%>%mutate(rank = row_number())


ggplot(df,aes(x=rank,y=avg_log2FC))+
  geom_point(aes(fill=direction),shape=21,color='transparent')+
  facet_wrap(vars(celltype) ,nrow = 1,scales='free')+
  scale_fill_manual(values = c(U='#da2d20',D='#33b2e8'))+
  theme(panel.border = element_rect(fill = 'transparent',colour = 'black'),
        panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        panel.background =element_blank(),
        plot.background =element_blank()#,
        #plot.title = element_text(hjust = 0.5,vjust = 1,size=15,face = 'bold'),
        # axis.ticks = element_blank(),
        # axis.text = element_blank(),
        # axis.title.x =element_blank(),
        # axis.title.y =element_blank()#,
        # legend.text=element_blank(),
        # legend.title=element_blank(),
        # legend.position = "none"
  )
ggsave(paste0(DEG_dir,'DEG_sc.pdf'),width = 15,height = 3)

#02.1 GO display--------------------------
for (drc in c('U','D')) {
  df_GO <- read.csv(paste0(DEG_dir,'GO_',drc,'.csv'))
  df_GO <- melt(df_GO,id.vars = 'Description',variable.name ='celltype',value.name='logP' )
  if (drc=='U') {
    cor_slt <- '#da2d20'
  }else {cor_slt <- '#33b2e8'}
  df_GO$logP <- -df_GO$logP
  dat2 <- df_GO
  #dat2$logP <- -dat2$logP
  dat2$celltype <- factor(dat2$celltype,levels = c("Epi","EC","Fib",
                                                   "Mye","TC","BC") )
  dat3 <- cast(dat2,Description~celltype,value = 'logP')
  dat3[is.na(dat3)] <- 0
  dat4 <- as.data.frame(dat3)
  rownames(dat4) <- dat4$Description
  dat4[dat4>0] <- 1
  dat4[dat4<0] <- -1
  dat4 <- dat4%>%dplyr::mutate(pat_num=rowSums(as.data.frame(lapply(dplyr::select(dat4,2:c(ncol(dat4)) ), as.numeric)))   )
  identical(dat3$Description,rownames(dat4))
  dat3$pat_num1 <- dat4$pat_num
  dat3 <- dat3%>%dplyr::mutate(pat_num=rowSums(dplyr::select(.,2:c(ncol(dat3)-1) )),
                               sort_order = purrr::pmap_dbl(dplyr::select(.,2:c(ncol(dat3)-1)   ), ~ { sum(2^(length(c(...)) - which(c(...) != 0))) })
  )%>%
    arrange( desc(pat_num1),-sort_order,desc(pat_num))
  dat3 <- as.data.frame(dat3)
  dat2$Description <- factor(dat2$Description,levels = rev(dat3$Description))
  ggplot(dat2,aes(x=celltype,y=Description,color=logP))+
    geom_point(size=5)+
    scale_color_gradient2(low ='white' ,
                          high = cor_slt)+
    theme(panel.border = element_rect(fill = 'transparent',colour ='black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major=element_blank(),
          panel.background = element_rect(fill = NA,colour = NA),
          # plot.background = element_rect(fill = 'transparent',colour = NA),
          # axis.line.x = element_line(size=0.2, colour = "black"),
          axis.ticks= element_line( colour = "black"),
          # axis.ticks.x = element_line(size=0.2, colour = "black"),
          axis.text.y = element_text(size=0.2),
          axis.text.x=element_text(),
          # axis.title.x = element_text(size=8,color="black",face="plain"),
          axis.title=element_blank(),
          # #legend.text=element_blank(),
          # legend.title=element_blank(),
          #legend.position = "none",
          # plot.title = element_text(hjust = 0.5,size=12),4
          
          plot.margin=unit(c(0,0,0,0),'lines'))
  
  ggsave(paste0(DEG_dir,'GO_',drc,'.pdf'),width = 2.8,height = 4.5)
  
}

#03 monocle analysis--------------------------------------------
mncl <- subset(seurat_rdscell_name=='Epi')
exp_mat <- as.matrix(GetAssayData(mncl, slot='data'))
exp_mat <- exp_mat[rowSums(exp_mat)!=0,]
pd <- data.frame(mncl[[]])
fd <- data.frame(gene_short_name=rownames(exp_mat))
rownames(fd) <- fd$gene
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)

mncl= FindVariableFeatures(mncl)
varGenes = VariableFeatures(mncl)

ordering_genes = head(varGenes, 1000)
cds_obj <- newCellDataSet(as(exp_mat, "sparseMatrix"),
                          phenoData = pd,
                          featureData = fd)
cds_obj <- estimateSizeFactors(cds_obj)
cds_obj <- estimateDispersions(cds_obj)
cds_obj <- setOrderingFilter(cds_obj, ordering_genes)
cds_obj <- reduceDimension(cds_obj, max_components = 2,
                           method = 'DDRTree')
cds_obj <- orderCells(cds_obj)

saveRDS(cds_obj,paste0(monocle_dir,"Epithelial_monocle.rds"))


#03.1 traj_display-----------------------------
p=plot_cell_trajectory(cds_obj,color_by="Pseudotime")+
  scale_color_gradientn(colours = scales::viridis_pal(option = "D")(20))
ggsave(paste0(monocle_dir,'psedutime_traj.png'),plot=p,width = 4,height = 3.5)
ggsave(paste0(monocle_dir,'psedutime_traj.pdf'),plot=p,width = 4,height = 3.5)

p=plot_cell_trajectory(cds_obj,color_by = "spp")+
  facet_wrap(spp~.,ncol=1)+
  scale_color_manual(values= c('#FF6347','#48D1CC') )
ggsave(paste0(monocle_dir,'psedutime_dese.png'),plot=p,width = 4,height = 6)
ggsave(paste0(monocle_dir,'psedutime_dese.pdf'),plot=p,width = 4,height = 6)


plotdf=data.frame(Pseudotime=pData(cds_obj)$Pseudotime,cluster=pData(cds_obj)$cell_new,state=pData(cds_obj)$State,age=pData(cds_obj)$spp)
plotdf=subset(plotdf,plotdf$Pseudotime != 'Inf')
write.csv(plotdf,paste0(monocle_dir,'pse_time_prop.csv'),row.names = F)

plotdf <- read.csv(paste0(monocle_dir,'pse_time_prop.csv'))
plotdf$age <- factor(plotdf$age,levels = c('Normal','Tumor'))
plotdf$state <- factor(plotdf$state,levels = c('1','2','3'))
ggplot(plotdf, aes(x=Pseudotime,y=state,fill=age))+
  geom_density_ridges(alpha=1,#transparent degree
                      scale=2,panel_scaling=F,
                      #color='black',
                      na.rm=TRUE,
                      quantile_lines=TRUE,
                      vline_linetype=2,
                      vline_color='black',
                      quantile_fun=function(x,...)mean(x)
                      
  ) +
  scale_y_discrete("")+
  scale_fill_manual(values=c('#FF6347','#48D1CC'))+
  theme_minimal()+ 
  geom_rug(aes(color=plotdf$age))+
  theme(panel.grid = element_blank())
ggsave(paste0(monocle_dir,'psedutime_prop.png'),plot=p,width = 4,height = 3)
ggsave(paste0(monocle_dir,'psedutime_prop.pdf'),plot=p,width = 4,height = 3)

#03.2 traj_gene pattern-----------------------------
BEAM_res=BEAM(cds_obj,branch_point = 1,cores = 20)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]

tmp1=plot_genes_branched_heatmap(cds_obj[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4, 
                                 cores = 20,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 branch_colors = c("#979797", '#FF6347','#48D1CC'), 
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
)

pdf(paste0(monocle_dir,"branched_heatmap.pdf"),width = 5,height = 6)
p <- tmp1$ph_res
print(p)
dev.off()































