library(glmGamPoi)
library(Seurat)
library(BPCells)
library(purrr)
library(tidyr)
library(dplyr)
library(qs)

TestPert<-function(seur=NULL,dat=NULL,meta=NULL,ctrl="NT_1",minCells=100,form=~Assign+batch,filter=.01)
{
  print("Prepare data")
  if(!is.null(seur)){
    meta=seur@meta.data
    dat=seur[["RNA"]]$counts
  }
  
  dat=dat[,meta$Num==1]
  meta=meta[meta$Num==1,]
  
  print("Get perts to test")
  vals=table(meta[,"Assign"])
  perts=names(vals)[vals>minCells]
  perts=perts[perts!=ctrl]
  perts=sample(perts)
  out=map(perts,function(x){
    print(x)
    gene=strsplit(x,"_")[[1]]
    gene=paste(gene[1:length(gene)-1],collapse="_")
    print(gene)
    meta2=meta[meta$Assign %in% c(ctrl,x),]
    dat2=as.matrix(dat[,rownames(meta2)])
    mn=rowMeans(dat2>0)
    dat2=dat2[mn>filter,]
    genes=rownames(dat2)
    meta2["Assign"]=factor(meta2[,"Assign"],levels=c(ctrl,x))
    mod_mat=model.matrix(form,meta2)
    coef=grep(x,colnames(mod_mat),value=T)[1]
    print(coef)
    print(dim(dat2))
    tic()
    fit=glmGamPoi::glm_gp(dat2, design = mod_mat, verbose=TRUE,subsample=T)
    toc()
    
    mrk=test_de(fit, coef, pval_adjust_method="BH");
    mrk["Gene"]=genes;mrk=mrk[order(mrk$pval),]
    print(head(mrk))
    print(sum(mrk$adj_pval<.05))
    print(mean(mrk$pval))
    print(quantile(mrk$pval))
    print(mrk[mrk$Gene==gene,])
    print("")
    return(mrk)
  })
  names(out)=perts
  return(out)
  
}


if(!interactive())
{
  args <- commandArgs(trailingOnly=T)
  mat <- open_matrix_dir(dir =args[1])
  meta=qread(args[2])
  mat=mat[,meta$Num==1]
  meta=meta[meta$Num==1,]
  #ctrl=args[2]
  ctrl="NT_1"
  mrk=TestPert(dat=mat,meta=meta,ctrl=ctrl,minCells=100,form=~Assign+orig.ident,filter=.01)
  qsave(mrk,"results.DE.qs")
}


