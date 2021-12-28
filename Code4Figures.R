require(tidyverse)
require(clusterProfiler)
require(amap)
require(cluster)
require(Cairo)
require(pheatmap)
require(corrplot)
require(e1071)
require(RColorBrewer)
require(stringr)
library(survival)
library(survminer)
library(mixOmics)
library(ggridges)
library(DESeq2)
library(glmnet)
library(purrr)
library(caret)
library(spls)

#data preparation
dds_os2006=read_tsv("Norm_count_table.tsv") %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("GENE")
dds_os2006=dds_os2006[which(apply(dds_os2006,1,max)>100),]
dds_os2006=dds_os2006[-caret::nearZeroVar(t(dds_os2006),saveMetrics=F),]
DDS_libnorm_100r_nzv_log1p_center_scale=log1p(dds_os2006) %>% t %>% scale(scale=T) %>% t



#BIODICA loading and computing of the ICs with the DDS_libnorm_100r_nzv_log1p_center_scale matrix

# To convert in shell output from BIODICA to load in R:
# for i in `ls *.xls`; do n=`echo $i | sed -e "s/.xls//"`; ssconvert -O 'separator=; format=raw quote= '  $i $n.txt; sed "s/,/./g" $n.txt >t ; mv t $n.txt ; done

#list=c(RNA_OS2006=79,CGH_OS2006=73,Nanostring_trainingSet=70,Nanostring_PredSet=96,RNA_TARGET=82,`CGH_OS2006&RNA_OS2006`=73,`Nanostring_trainingSet&RNA_OS2006`=70)
#upset(fromExpression(list), order.by = "freq")



#setwd("Your folder")
fname="DDS_libnorm_100r_nzv_log1p_center_scale"
AnalFile="Figures"
decomp_jade=c()
decomp_jade$S=read_delim(paste0(fname,"_ica_S.txt"),delim=";")
decomp_jade$M=read_delim(paste0(fname,"_ica_A.txt"),delim=";")
decomp_jade$S %>% column_to_rownames("PROBE") %>% as.data.frame ->decomp_jade$S
decomp_jade$M %>% column_to_rownames("SAMPLE") %>% as.data.frame ->decomp_jade$M
mat=read_tsv(paste0("./",fname,".tsv"))
mat %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("SYMBOL") -> mat_df
Annot3=read_tsv("./Annot3_wo.tsv")
Annot3 %>% as.data.frame %>% remove_rownames %>% column_to_rownames("Full_name") -> Annot3
load("./msigdb.Rdata")
load("./Immnuno_sig")
load("pfs_nob_rna.Rdata")
load("surv_nob_rna.Rdata")
load("glm_final_lasso_mat_df2_035.RData")
dir.create(AnalFile)

#Figure 1A ICs correlation matrix and dendograms for Patient and ICs
CairoPDF(paste0(AnalFile,"/pheatmap",fname,"_correlation"),width=50,height=20)
decomp_jade$M %>% as.data.frame %>% cor(method="pearson")  %>% corrplot(method="color",order="hclust",hclust.method="ward.D")
plot(hclust(amap::Dist(decomp_jade$M,method="kendall"),method="ward.D"),main="Patient Tumours")
plot(hclust(amap::Dist(t(decomp_jade$M),method="kendall"),method="ward.D"),main="Gene Modules")
dev.off()
#Enrichment have been calculated with the gene in table S1 within the gprofilR and msigdb website

#contributive geneset to ICs
decomp_jade$M %>% as.data.frame %>% as.data.frame%>% t%>% amap::Dist(method="kendall")  %>%hclust(method="ward.D") %>% .$order ->Morder
fdr="3Mean"
decomp_jade_lfdr_logic =apply(decomp_jade$S,2,function(x){abs((x-mean(x))/sd(x))>3})
gene_list_decomp_jade <- map2(as.data.frame(decomp_jade_lfdr_logic),colnames(as.data.frame(decomp_jade_lfdr_logic)),function(x,y){data.frame("gene"=rownames(decomp_jade$S)[x],"Sig"=as.data.frame(decomp_jade$S)[x,y])})


#Figure 1B : Simple co-expression network inference from contributive gene to ICs. gene_ICA2.5_TO_CPS.tsv was loaded in cytoscape and cluster of gene was annotated with REAECTOMEviz plugin
purrr::map(gene_list_decomp_jade,~ .x$gene ) %>% unlist %>% unique -> list
mat %>% dplyr::filter(SYMBOL %in% list) %>% as.data.frame %>% remove_rownames %>% column_to_rownames("SYMBOL") %>% Dist(method="pearson") ->dist
dist %>% as.matrix %>% as.data.frame  %>% rownames_to_column("GENE1") %>%  gather(key=GENE2,value=Dist,-GENE1) ->dist_gath
dist_gath %>% dplyr::filter(abs(Dist)<0.25 & GENE1!=GENE2)  %>% write_tsv(path="network2.5_test_pearson_025.tsv")

#Correspondance between Gene and associated ICs for network coloring in Cytoscape
map2(gene_list_decomp_jade,names(gene_list_decomp_jade),~ cbind(as.vector(.x$gene),.y)%>% as.data.frame) %>% bind_rows -> gene_ICA2.5_TO_CPS
gene_ICA2.5_TO_CPS %>% as.data.frame %>% dplyr::select(GENE=1,ICS=2) %>% group_by(GENE) %>% mutate(ICS=paste(ICS,collapse=",")) %>% distinct -> gene_ICA2.5_TO_CPS
write_tsv(gene_ICA2.5_TO_CPS,path="gene_ICA2.5_TO_CPS.tsv")

#Figure 2A: heatmap of ICs and related clinical ANNOTATIONs
anno_colors <- list(sex= c(Masculin="blue",Feminin="red"),
Meta_bis=c(Definite_metastases="grey52",Localised_Osteosarcoma="light blue"),
Meta_long=c(Localised="light blue",Metastases="black"),
htumclass=c("."="white",`Equal or greater than 10 cm`="darkgoldenrod4",`Less than 10 cm`="darkgoldenrod1"),
age=brewer.pal(9,"Blues"),
etat2=c("0"="light grey","1"="firebrick"),
rep_histo=c("."="white",GR="light blue",PR="darkorchid2"),
pub=c("."="white",Prépubertaire="yellow",Intrapubertaire="darkolivegreen3",Postpubertaire="darkolivegreen"),
relapse=c(N="light blue",R="firebrick2"),
sstyp_ana=c(Autre="light grey",`De surface de haut grade`="darkorchid1",Télangiectasique="lightpink", `Forme commune chondroblastique`="gold1",`Forme commune fibroblastique`="dodgerblue1",`Forme commune sans précision`="darkslategray1",`Forme commune ostéoblastique`="firebrick1")
)
CairoPDF(paste0(AnalFile,"/Heatmap_pearson_ward_stable_cps_M.pdf"),width=20,height=14)
decomp_jade$M[,c(50,49,48,47,46,32,44,45,37,41,43,40,33,39,35,36,30,29,3,38,42,25,22,21,34,14,28,20,4,31,26,18,15,27)] %>% t  %>% pheatmap(annotation=Annot3,clustering_method="ward.D",clustering_distance_cols=amap::Dist(t(.),method="pearson"),clustering_distance_rows=amap::Dist(.,method="pearson"))
dev.off()

#Two groups of tumors
dendo2=  decomp_jade$M[,c(50,49,48,47,46,32,44,45,37,41,43,40,33,39,35,36,30,29,3,38,42,25,22,21,34,14,28,20,4,31,26,18,15,27)]  %>% Dist(method="pearson") %>% hclust(method = "ward.D") %>% cutree(k=2)

#Figure 2B and 2C: OS and PFS curves between G1 and G2
CairoPDF(paste0(AnalFile,"/Surv_group2.pdf"))
ggsurvplot(surv_fit(surv_nob_rna_NO_MOR_SO ~ group,data.frame(group=dendo2)),pval=T,conf.int=T)
ggsurvplot(surv_fit(pfs_nob_rna_NO_MOR_SO ~ group,data.frame(group=dendo2)),pval=T,conf.int=T)
dev.off()

#Figure 2E: Code to estimate contribution of our stratification, PCA for the stable ICs
PC=prcomp(decomp_jade$M[,c(50,49,48,47,46,32,44,45,37,41,43,40,33,39,35,36,30,29,3,38,42,25,22,21,34,14,28,20,4,31,26,18,15,27)])
PCL=as.data.frame(PC$rotation)
PC=as.data.frame(cbind(PC$x,data_frame(Annot3,Group=dendo2)))


#Figure 2E: PCA projection Of ICs with different Annotations
CairoPDF(paste0(AnalFile,"/PCA_dendo2_etat2.pdf"),width=8, height=5)
ggplot()+geom_segment(data=PCL,aes(x=0,y=0,xend=PC1*2,yend=PC2*2),arrow=arrow(length=unit(0.1,"inches")),alpha=0.5) + geom_point(data=PC,aes(PC1,PC2,color=as.factor(Group),alpha=as.factor(etat2),size=2)) +  scale_alpha_discrete(range = c(0.2, 1))  +theme_bw()
ggplot()+geom_segment(data=PCL,aes(x=0,y=0,xend=PC1*2,yend=PC2*2),arrow=arrow(length=unit(0.1,"inches")),alpha=0.5) + geom_point(data=PC,aes(PC1,PC2,color=as.factor(Group),alpha=as.factor(Meta_long),size=2)) +  scale_alpha_discrete(range = c(0.2, 1))  +theme_bw()
ggplot()+geom_segment(data=PCL,aes(x=0,y=0,xend=PC1*2,yend=PC2*2),arrow=arrow(length=unit(0.1,"inches")),alpha=0.5) + geom_point(data=PC,aes(PC1,PC2,color=as.factor(Group),alpha=as.factor(relapse),size=2)) +  scale_alpha_discrete(range = c(0.2, 1))  +theme_bw()
ggplot()+geom_segment(data=PCL,aes(x=0,y=0,xend=PC1*2,yend=PC2*2),arrow=arrow(length=unit(0.1,"inches")),alpha=0.5) + geom_point(data=PC,aes(PC1,PC2,color=as.factor(Group),alpha=as.factor(rep_histo),size=2)) +  scale_alpha_discrete(range = c(0.2, 1))  +theme_bw()
dev.off()


#PLS discriminant analysis of ICs
list.keepX <- c(1:10,  seq(20, 30, 10))
tune.splsda.srbct <- tune.splsda(decomp_jade$M[,c(50,49,48,47,46,32,44,45,37,41,43,40,33,39,35,36,30,29,3,38,42,25,22,21,34,14,28,20,4,31,26,18,15,27)],as.factor(dendo2), ncomp = 6, validation = 'Mfold', folds = 5, progressBar = TRUE, dist = 'max.dist', measure = "BER",test.keepX = list.keepX, nrepeat = 10, cpus = 4,near.zero.var=T,max.iter=1000)
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]
plot(tune.splsda.srbct, col = color.jet(6))
splsda.srbct <- mixOmics::splsda(decomp_jade$M[,c(50,49,48,47,46,32,44,45,37,41,43,40,33,39,35,36,30,29,3,38,42,25,22,21,34,14,28,20,4,31,26,18,15,27)],as.factor(dendo2), ncomp = ncomp, keepX = select.keepX)
perf.srbct <- perf(splsda.srbct, validation = "Mfold", folds = 5,
               dist = 'max.dist', nrepeat = 10,
               progressBar = FALSE)
CairoPDF(paste0(AnalFile,"/PLSDA_figures.pdf"))
plotLoadings(splsda.srbct, comp = 1, title = 'Loadings on comp 1', contrib = 'max', method = 'mean')
dev.off()


#FIgure 3A has been produced in Cytoscape by loading the differential expression of genes between G1 and G2 generated with DESeq2 as generated for figure 3C

#Reloading of Network after annotation in cytoscape and estimation of G1/G2 logfolchange of sub-Networks
Network_nodes=read_csv("./network2.5_test_abspearson_025_Annotateeed_by_Cytoscape.csv")

#Figure 3B: Ridges
CairoPDF(paste0(AnalFile,"/RIDGES_module_reactome_network.pdf_FULL"),width=8,height=8);
Network_nodes %>% dplyr::filter(module %in%  c(1,4,13,6,2,5,18,51) | ICS %in% c("IC26","IC1")) %>% mutate(module=ifelse(ICS=="IC26","IC26",module)) %>% mutate(module=ifelse(ICS=="IC1","IC1",module))%>% mutate(module=as.factor(module)) %>% mutate(module=fct_relevel(module,rev(c("IC26","1","4","13","6","2","5","18","IC1","51")))) %>% ggplot(aes(x=log2FoldChange,y=as.factor(module),fill=..x..)) +geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.)  + scale_y_discrete(expand = c(0.01, 0))  + scale_fill_gradientn(colors=brewer.pal(7,"RdBu")) +  xlim(-2.5,2.5) +theme_ridges()
dev.off()

#Figure 3C
#DEGs by DESeq2
DDS_os2006=read_tsv("Norm_count_table.tsv") %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("GENE")
DESeqDataSetFromMatrix(countData=round(DDS_os2006),colData=data.frame(Annot3,dendo2=as.factor(dendo2)),design=~dendo2) -> DDS_os2006
DDS_os2006<-DESeq(DDS_os2006); resuls=results(DDS_os2006,tidy=T);

#Volcanoplot
CairoPDF(paste0(AnalFile,"/Vulcano_G1_vs_G2_BW.pdf"),width=10,height=8);resuls %>% dplyr::filter(row!="RHOXF2B") %>% ggplot(aes(log2FoldChange,-log10(padj))) + geom_point(alpha=0.3,color=ifelse(resuls[resuls$row!="RHOXF2B",]$padj<0.05,ifelse(resuls[resuls$row!="RHOXF2B",]$log2FoldChange>0,"blue","red"),"black")) + geom_text(aes(log2FoldChange,-log10(padj),label=ifelse(resuls[resuls$row!="RHOXF2B",]$padj<0.001,row,ifelse(resuls[resuls$row!="RHOXF2B",]$padj<0.05 & grepl("HIS",row),row,""))),size=2,vjust=-1.5) +theme_bw()
dev.off()

#cbind(apply(log1p(dds_df),1,mean),apply(log1p(dds_df),1,sd)) -> Mean_and_SD_100_log1p
#write_tsv(rownames_to_column(Mean_and_SD_100_log1p,"SYMBOL"),path="Mean_and_SD_100_log1p")

#Figure 4A : file Group2_W.txt was generated as described in the method section
CGH=read_tsv("./CGH_RNAclustGroups_20190125/pcut_0.01/_Group2_W/Group2_W.txt")

 CairoPDF(paste0(AnalFile,"/CGH_group2.pdf"),width=30,height=2)
 CGH %>% mutate(Chrom=as.numeric(gsub("Chr","",Chrom))) %>% ggplot(aes(Start,Group2.1.MedianL2R-Group2.2.MedianL2R,color=(-log10(RawP))))+ geom_point(aes(size=as.numeric(ifelse(RawP<0.01,4,2)),shape=ifelse(RawP<0.01,18,16)),alpha=0.5) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x") +scale_color_gradient(low="dark blue",high="red") +theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) + scale_size_identity()+scale_shape_identity()
 dev.off()

 CairoPDF(paste0(AnalFile,"/CGH_group2.pdf"),width=30,height=5)
CGH %>% dplyr::select(Group2.1.MedianL2R,Group2.2.MedianL2R,Chr,Start,RawP) %>% gather(key=Group,value=MedianL2R,-Chr,-Start,-RawP) %>% ggplot(aes(x=Start))+geom_area(aes(y=MedianL2R,fill=Group),alpha=0.5,position = 'identity') +geom_line(aes(y=MedianL2R,linetype=Group),size=0.5) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x")+theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) + scale_size_identity()+scale_shape_identity()
CGH %>% ggplot(aes(x=Start))  + geom_point(aes(y=Group2.2.MedianL2R-Group2.1.MedianL2R,color=(-log10(RawP)),size=as.numeric(ifelse(RawP<0.01,4,2)),shape=ifelse(RawP<0.01,18,16)),alpha=0.7) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x") +scale_color_gradient(low="dark blue",high="red") +theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) + scale_size_identity()+scale_shape_identity()
CGH %>% dplyr::select(Group2.1.MedianL2R,Group2.2.MedianL2R,Chr,Start,RawP) %>% gather(key=Group,value=MedianL2R,-Chr,-Start,-RawP) %>% ggplot(aes(x=Start))+geom_area(aes(y=MedianL2R,fill=Group),alpha=0.5,position = 'identity') +geom_line(aes(y=MedianL2R,linetype=Group),size=0.5) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x")+theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") + scale_size_identity()+scale_shape_identity()
CGH %>% ggplot(aes(x=Start)) +geom_point(aes(y=Group2.2.MedianL2R-Group2.1.MedianL2R,color=(-log10(RawP)),size=as.numeric(ifelse(RawP<0.01,4,2)),shape=ifelse(RawP<0.01,18,16)),alpha=0.7) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x") +scale_color_gradient(low="dark blue",high="red") +theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") + scale_size_identity()+scale_shape_identity()
CGH %>% ggplot(aes(x=Start,y=(-log10(RawP))))  + geom_area() +geom_hline(yintercept=2,linetype="dashed",color="red",size=1)+ facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x") +scale_color_gradient(low="dark blue",high="red") +theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") + scale_size_identity()+scale_shape_identity()
CGH %>% dplyr::select(Group2.1.MedianL2R,Group2.2.MedianL2R,Chr,Start,RawP) %>% gather(key=Group,value=MedianL2R,-Chr,-Start,-RawP) %>% ggplot(aes(x=Start))+geom_area(aes(y=MedianL2R,fill=Group),alpha=0.5,position = 'identity') +geom_line(aes(y=MedianL2R,linetype=Group),size=0.5) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x")+theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) + scale_size_identity()+scale_shape_identity()
  CGH %>% mutate(Chrom=as.numeric(gsub("Chr","",Chrom))) %>% ggplot(aes(x=Start))+geom_area(aes(y=Group2.1.MedianL2R),fill="red",alpha=0.5) +geom_area(aes(y=Group2.2.MedianL2R),fill="blue",alpha=0.5)+geom_line(aes(y=Group2.1.MedianL2R),size=0.5)+geom_line(aes(y=Group2.2.MedianL2R),size=0.5,linetype="dashed") + geom_point(aes(y=Group2.2.MedianL2R-Group2.1.MedianL2R,color=(-log10(RawP)),size=as.numeric(ifelse(RawP<0.01,4,2)),shape=ifelse(RawP<0.01,18,16)),alpha=0.7) +facet_grid(~ Chr,scale="free_x",margins=F,switch="x",space="free_x") +scale_color_gradient(low="yellow",high="red") +theme_bw() +theme(panel.spacing = unit(0, "lines"),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) + scale_size_identity()+scale_shape_identity()
dev.off()
  
  
#Figure 4B
Annot=read_tsv("./TARGET_OS_Discovery_ClinicalData_20170525.txt") %>% separate(`TARGET USI`,into=c("Meta","Study","Name"),sep="-")

files=dir(pattern="quantification.txt")
purrr::map(files,function(x){read_tsv(x )%>% add_column(FileName=rep(x,45476))}) %>% purrr::reduce(rbind)  %>% separate(FileName,into=c("Meta","Study","Name","F","format"),sep="-") %>% dplyr::select(-format,-F,-Meta,-Study)-> expression
Annot %>% inner_join(expression,by="Name") ->All_data
All_data %>% dplyr::select(ensembl_gene_id,tpm,Name) %>% group_by(ensembl_gene_id) %>% mutate(M=mean(log1p(tpm)),S=sd(log1p(tpm))) %>% ungroup %>% mutate(SC=(log1p(tpm)-M)/S)%>% dplyr::select(ensembl_gene_id,SC,Name) %>% spread(key=Name,value=SC)   %>% dplyr::filter(!is.na(`0A4HLD`))  %>% write_tsv(path="./DDS_log1p_center_scale_by_row_EST.tsv")

All_data %>% dplyr::select(ensembl_gene_id,est_counts,Name) %>% spread(key=Name,value=est_counts) %>% as.data.frame %>% remove_rownames %>% column_to_rownames("ensembl_gene_id") -> expression_EST
dds=DESeqDataSetFromMatrix(countData=round(expression_EST),colData=data.frame(Full_name=colnames(expression_EST)),design=~1)
dds=estimateSizeFactors(dds)
dds_df=as.data.frame(counts(dds,normalized=T ))

bitr(rownames(dds_df),fromType="GENEID",toType="SYMBOL",OrgDb="EnsDb.Hsapiens.v75",drop=T) ->mat_to_SYMBOL
read_tsv("./Mean_and_SD_100_log1p") -> Mean_and_SD_100_log1p
dds_df %>% rownames_to_column("GENEID") %>% as.tibble%>% inner_join(mat_to_SYMBOL,by="GENEID") %>% inner_join(Mean_and_SD_100_log1p,by="SYMBOL")  %>% gather(key=Sample,value=EST_norm,-SYMBOL,-V1,-V2,-GENEID) %>% mutate(EST_norm_log_scaled=(log1p(EST_norm)-V1)/V2)  %>% dplyr::select(Sample,EST_norm_log_scaled,SYMBOL,GENEID) %>% mutate(REF=paste0(SYMBOL,"_",GENEID)) %>% dplyr::select(Sample,EST_norm_log_scaled,REF) %>% spread(key=Sample,value=EST_norm_log_scaled) -> DDS #write_tsv("DDS_log1p_20_nzv_center_scaled_byOS2006_by_row_EST.tsv")
Symbol=read.table("Symbol.txt")

cv.glmnet(t(mat_df[rownames(mat_df) %in% Symbol$SYMBOL,]), as.factor(dendo2), type.measure="mse", alpha=0.35,family="binomial") ->glm_lasso_mat_df2_035
glmnet(t(mat_df[rownames(mat_df) %in% Symbol$SYMBOL,]), as.factor(dendo2),lambda=0.5139,  alpha=0.35,family="binomial") ->glm_final_lasso_mat_df2_035

DDS %>% separate(REF,into=c("SYMBOL","GENEID"),sep="_")%>% dplyr::select(-GENEID) %>% group_by(SYMBOL) %>% summarise_all(funs(mean)) %>% dplyr::filter(SYMBOL %in% rownames(coef(glm_final_lasso_mat_df2_035)))    %>% as.data.frame() %>% remove_rownames %>% column_to_rownames("SYMBOL") %>% t %>% scale  %>%
  predict.glmnet(glm_final_lasso_mat_df2_035,s=0.5139,newx=.,type="response")-> glm_final_lasso_mat_df2_035_predict
glm_final_lasso_mat_df2_035_predict_disc=data.frame(dendo=as.factor(ifelse(scale(glm_final_lasso_mat_df2_035_predict[1:82])>0,1,2)))
sfit=survfit(surv[1:82] ~ dendo, data=glm_final_lasso_mat_df2_035_predict_disc)
CairoPDF(paste0(AnalFile,"/Suvr_lim_3000_glmnet_035_1SE.pdf"))
ggsurvplot(sfit,glm_final_lasso_mat_df2_035_predict_disc,pval=T,conf.int=T,risk.table=T,xlim=c(0,3000))
dev.off()

#Figure 4C
mean_sd<- Mean_and_SD_100_log1p
Annot=read_tsv("OS2006_clinical_Nanostring.tsv")
All_data=read_tsv("Data_nanostring.tsv")
RNA_TPM=read_tsv("./Norm_count_table.tsv")
RNA=colnames(read_tsv("./Norm_count_table.tsv"))

Annot[match(colnames(All_data_df),Annot$numenreg),]->Annot
All_data %>% as.data.frame %>% remove_rownames %>% column_to_rownames("GENE") -> All_data_df
Annot[Annot$numenreg %in% RNA,] -> Annot_RNA
RNA_TPM %>% gather(key=Patient,value=TPM,-GENE) %>% inner_join(Annot,by=c("Patient"="Full_name"))%>% dplyr::select(GENE,Patient=numenreg,TPM)%>% mutate(Patient=as.character(Patient))%>% distinct_all() %>% inner_join(gather(All_data,key=Patient,value=Count,-GENE),by=c("Patient","GENE") ) ->RNA_vs_NANO

RNA_vs_NANO %>% group_by(GENE) %>% mutate(cor=cor(TPM,Count,method="pearson"),pval=cor.test(TPM,Count,method="pearson")$p.value) %>% dplyr::select(GENE,cor,pval) %>% distinct %>% dplyr::filter(pval<0.00001 & cor>0) -> Cor_RNA_NANO_GENE

All_data_df[rownames(All_data_df) %in% Annot_RNA$numenreg,]
All_data %>% dplyr::select(any_of(c("GENE",Annot_RNA$numenreg))) %>% dplyr::filter(GENE %in% c(Cor_RNA_NANO_GENE$GENE,"CCDC34 iso 1")) %>%  dplyr::filter(!grepl("POS_",GENE ) & !grepl("NEG_",GENE )& !GENE %in% c("CNOT1","EIF4G2","SF1","SLC39A1","SURF4")) %>% as.data.frame %>% remove_rownames %>% column_to_rownames("GENE") -> train_data
Annot3 %>% rownames_to_column("Full_name") %>% mutate(dendo2=dendo2) ->Annot3_tb
Annot_RNA %>% inner_join(Annot3_tb,by=c("Full_name")) %>% dplyr::filter(numenreg %in% colnames(train_data)) %>% dplyr::select(dendo2)  ->train_annot

All_data %>% dplyr::select(-one_of(as.character(Annot_RNA$numenreg))) %>% dplyr::filter(GENE %in% c(Cor_RNA_NANO_GENE$GENE,"CCDC34 iso 1")) %>%  dplyr::filter(!grepl("POS_",GENE ) & !grepl("NEG_",GENE ) & !GENE %in% c("CNOT1","EIF4G2","SF1","SLC39A1","SURF4")) %>% as.data.frame %>% remove_rownames %>% column_to_rownames("GENE") -> test_data

Annot%>% dplyr::filter(numenreg %in% colnames(test_data)) %>% dplyr::select(etat2) %>% mutate(etat2=ifelse(etat2=="Vivant",0,1)) ->test_annot

All_data %>% dplyr::filter(GENE %in% c(Cor_RNA_NANO_GENE$GENE,"CCDC34 iso 1")) %>%  dplyr::filter(!grepl("POS_",GENE ) & !grepl("NEG_",GENE )& !GENE %in% c("CNOT1","EIF4G2","SF1","SLC39A1","SURF4")) %>% as.data.frame %>% remove_rownames %>% column_to_rownames("GENE") ->All_data_filter

cv.glmnet((t((train_data))), as.factor(train_annot$dendo2), type.measure="mse", intercept=T,alpha=0.8,family="binomial",nfolds=10,nlambda=10000)  ->glm_cv
glmnet((t((train_data))), as.factor(train_annot$dendo2),lambda=glm_cv$lambda.min,  alpha=0.8,family="binomial",intercept=T) ->glm_fin
coef((glm_fin)) %>%  as.matrix() %>% as.data.frame %>% rownames_to_column("GENE")  %>% arrange(s0)
cbind(predict(glm_fin,s=glm_fin$lambda,newx=(t((train_data))),type="class"),train_annot) %>% table
predict(glm_fin,s=glm_fin$lambda,newx=(t((test_data))),type="class") -> test_data_fit
test_data_fit_disc=data.frame(dendo=test_data_fit[,1])
test_data_fit_disc$dendo=ifelse(test_data_fit_disc$dendo==1,2,1)
Annot%>% dplyr::filter(numenreg %in% colnames(test_data))-> test_annot_all

OS_test=Surv(time=test_annot_all$del_dn_enregy *365.25 ,event=ifelse(test_annot_all$etat2=="Vivant",0,1))
PFS_test=Surv(time=test_annot_all$dnum_PFS*365.25,event=ifelse(test_annot_all$relapse==".",0,1))
pval2=surv_pvalue(surv_fit(OS_test ~ dendo,data.frame(group=test_data_fit_disc)))$pval
p1=ggsurvplot(surv_fit(OS_test ~ dendo,data.frame(group=test_data_fit_disc)),pval=T,conf.int=T,risk.table=T)
p2=ggsurvplot(surv_fit(PFS_test ~ dendo,data.frame(group=test_data_fit_disc)),pval=T,conf.int=T,risk.table=T)
CairoPDF(paste0(AnalFile,"/OS_PFS_pearson_P10_5_A0.8_IT.pdf"))
p1
p2
dev.off()

#Figure 4D
mat=read_tsv("TPM_gene_10_log1p_nzv.tsv")
read_tsv("./Mean_and_SD_100_log1p") -> mean_sd
mat %>% dplyr::filter(GENE %in% mean_sd$SYMBOL) ->mat_filt
mat_filt[,!colnames(mat_filt) %in% c("X205","X256","X303")] ->mat_filt
read_tsv("./Mappy_MoscatoPed_SampleAnnotations_2018_11_02.txt") -> Annot
load("./glm_final_lasso_mat_df2_035.RData")
mat_filt %>% inner_join(mean_sd,by=c("GENE"="SYMBOL")) %>% gather(key=Sample,value=TPM,-GENE,-V1,-V2) %>% mutate(Z=(TPM-V1)/V2,Sample=(Sample)) %>% dplyr::select(GENE,Sample,Z) %>% spread(key=Sample,value=Z) %>% dplyr::filter(GENE  %in%  rownames(coef(glm_final_lasso_mat_df2_035)))%>% full_join(data.frame(GENE=rownames(coef(glm_final_lasso_mat_df2_035))),by="GENE",fill=0) %>% dplyr::filter(!grepl("Intercept",GENE))%>% replace(.,is.na(.),0)%>% as.data.frame %>% remove_rownames %>% column_to_rownames("GENE") -> mat_filt_Z_df
predict.glmnet(glm_final_lasso_mat_df2_035,s=glm_final_lasso_mat_df2_035$lambda,newx=t(scale(as.matrix(mat_filt_Z_df)[order(rownames(mat_filt_Z_df)),])),type="response") ->mat_filt_Z_df_predict
mat_filt_Z_df_predict_disc=data.frame(dendo=as.factor(ifelse(scale(mat_filt_Z_df_predict)>0,2,1)))
CairoPDF(paste0(AnalFile,"/Barplot_predict_model_glmnet_035alpha.pdf"));
cbind(mat_filt_Z_df_predict,mat_filt_Z_df_predict_disc) %>% rownames_to_column("Sample") %>%mutate(Sample=as.numeric(Sample)) %>% inner_join(Annot,by=c("Sample"="Patient_id")) %>% as.tibble %>% dplyr::filter(`Tumor type`=="Osteosarcoma") %>% dplyr::select(dendo) %>% ggplot(aes(dendo)) +geom_bar(stats=count) ;
dev.off()


# Figure S1: Code to modelize CNV with RNA glmnet
ADMIX_RNA_df=read.table("ADMIX_RNA_df.txt")
M_df=decomp_jade$M 
      
      IC2CGH_glmnet=function(x,y){
        names(x)<-rownames(M_df)
        IC=data.frame(Sample=names(x),x)
        cv.glmnet(t(ADMIX_RNA_df), x[names(x) %in% colnames(ADMIX_RNA_df)], type.measure="mse", alpha=0.8,family="gaussian") ->glm_lasso_mat_df_IC_ADMIX
        glmnet(t(ADMIX_RNA_df), x[names(x) %in% colnames(ADMIX_RNA_df)],lambda=glm_lasso_mat_df_IC_ADMIX$lambda.min,  alpha=0.8,family="gaussian") ->glm_final_lasso_mat_df_IC_ADMIX
        coef((glm_final_lasso_mat_df_IC_ADMIX)) %>%  as.matrix() %>% as.data.frame %>% rownames_to_column("CNV")%>% dplyr::filter(s0!=0)->  IC2CGH_0.8
        IC2CGH_0.8 %>% arrange(s0) %>% tail(1) ->bCNV
        ADMIX_RNA_df[rownames(ADMIX_RNA_df)==bCNV$CNV,] %>% rownames_to_column("CNVs") %>% gather(key=Sample,value=LogFoldChange,-CNVs)  %>% inner_join(IC ,by="Sample") ->data
        tryCatch( cor.test(data$x,data$LogFoldChange,method="spearman"), error=function(e){return(NULL)}) ->cordata
        data %>% ggplot(aes(x,LogFoldChange)) + geom_point() + geom_smooth() + ggtitle(paste0(bCNV$CNV,"_vs_",y)) + xlab(paste0("rho:",signif(as.numeric(cordata$estimate),digits=4)," pval:",signif(as.numeric(cordata$p.value),digits=4)))
      }
      CairoPDF(paste0(AnalFile,"/IC2CGH_0.8_rho_loess.pdf"))
      map2(M_df,colnames(M_df),function(x,y){IC2CGH_glmnet(x,y)})
      dev.off()

# Figure S2:
      
#Figure done in cytoscape with the logfoldchange estimated with DESSeq2 between G1 and G2 
    


