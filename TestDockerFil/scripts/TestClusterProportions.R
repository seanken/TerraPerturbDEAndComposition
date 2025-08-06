library(Seurat)
library(aod)
library(dplyr)
library(tidyr)
library(purrr)
library(qs)

TestWithBetaBin <- function(meta,batch="orig.ident",cond="Assign",ref="NT_1",CellType="seurat_clusters",form=cbind(Num, Tot-Num) ~ Assign+batch)
{
  meta=meta[,c(cond,CellType,batch)]
  colnames(meta)=c("Assign","CT","batch")
  tab1<-meta %>% group_by(CT,Assign,batch) %>% summarise(Num=length(CT)) %>% as.data.frame()
  tab2<-meta %>% group_by(Assign,batch) %>% summarise(Tot=length(CT)) %>% as.data.frame()
  tab2=do.call(rbind,map(unique(meta$CT),function(x){tmp=tab2;tmp["CT"]=x;return(tmp)}))
  tab=right_join(tab1,tab2)
  tab[is.na(tab$Num),"Num"]=0
  
  out=split(tab,tab$CT)
  modelFitFunction=function(form,data,link){betabin(form,random=~1,data=data)}
  res=map(out,function(x){
    x["Assign"]=factor(x$Assign,levels=unique(meta$Assign))
    x["Assign"]=relevel(x$Assign,ref=ref)
    fit=tryCatch({modelFitFunction(form,data=x,link=link)},error=function(cond){print(cond);return(NULL)})
    if(is.null(fit))
    {
      return(NULL)
    }
    coef=NULL
    
    coef=summary(fit)@Coef
    
    #coef=summaryAOD(fit)@Coef
    coef=data.frame(coef)
    #print(head(coef))
    colnames(coef)[4]="pval"
    coef["Test"]=rownames(coef)
    
    return(coef)
    
  })
  
  names(res)=names(out)
  res[sapply(res, is.null)]=NULL
  
  for(i in names(res)){
    res[[i]]["CT"]=i
  }
  
  ret=do.call(rbind,res)
  
  ret=data.frame(ret)
  ret=ret[order(ret$pval),]
  ret=ret[grep("Assign",ret$Test),]
  ret["Test"]=sub("Assign","",ret$Test)
  ret["padj"]=p.adjust(ret[,"pval"],"fdr")
  rownames(ret)=NULL
  return(ret)
}


GetPerc<-function(seur,cond="Assign",CellType="CT",verbose=T)
{
  genes=rownames(seur@assays$RNA@counts)
  seur@meta.data["CT"]=seur@meta.data[,CellType]
  seur@meta.data["Assign"]=seur@meta.data[,cond]
  
  out=map(1:length(genes),function(i){
    x=genes[i]
    if(i %% 100 == 0 & verbose) {
      print(paste("Processing gene", i, "of", length(genes)))
    }
    meta=seur@meta.data
    meta["Gene"]=seur@assays$RNA@counts[x,]
    tab<-meta %>% group_by(CT,Assign) %>% summarise(Perc=100*mean(Gene>0)) %>% as.data.frame()
    tab["Gene"]=x
  })
  
  tab=do.call(rbind,out)
  tab=tab %>% pivot_wider(names_from=Assign,values_from=Perc,values_fill = 0)
  out=split(tab,tab$CT)
  return(out)
  
}




if(!interactive())
{
  args <- commandArgs(trailingOnly=T)
  #args <- commandArgs()
  #mat <- open_matrix_dir(dir =args[1])
  meta=qread(args[1])
  #mat=mat[,meta$Num==1]
  meta=meta[meta$Num==1,]
  mrk=TestWithBetaBin(meta,batch="orig.ident",cond="Assign",ref="NT_1",CellType="seurat_clusters",form=cbind(Num, Tot-Num) ~ Assign+batch)

  qsave(mrk,"results.cluster.comp.qs")
}


