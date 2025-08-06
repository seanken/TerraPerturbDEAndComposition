library(nebula)
library(BPCells)
library(qs)

RunNebula<-function(dat,meta,form_fixed,sampleCol,cpc=.1,model="NBGMM",method="LN")
{
  print(dim(dat))
  print("Reorder")
  dat=dat[,order(meta[,sampleCol])]
  meta=meta[order(meta[,sampleCol]),]
  print("Run DE!")
  df = model.matrix(form_fixed, data=meta)
  print(head(df))
  print(head(meta[,sampleCol]))
  ##Need to make dat a matrix? or is BPCells ok?
  re = nebula(dat,meta[,sampleCol],pred=df,offset=meta$nCount_RNA,cpc=cpc,model=model,method=method)
  return(re)
}

##
##Gene_col is column with name of the gene targetted by the guide (or controls, etc)
##
RunNebula_pert<-function(seur=NULL,dat=NULL,meta=NULL,form_fixed=~Target+orig.ident,Gene_col="Target",sampleCol="Assign",cpc=.01,model="NBGMM",method="LN",ref="NT",minCells=100,toTest=c())
{
  if(!is.null(seur))
  {
    meta=seur@meta.data
    dat=seur[["RNA"]]$counts
  }
  meta["Target"]=map_chr(meta$Assign,function(x){s=strsplit(x,"_")[[1]];paste(s[1:(length(s)-1)],collapse="_")})
  targets=table(meta[,Gene_col])
  targets=names(targets[targets>minCells])
  
  if(!(ref %in% targets))
  {
    print("Not enough cells in reference group!")
    return(NULL)
  }
  if(length(toTest)>0)
  {
    targets=targets[targets %in% toTest]
  }
  targets=targets[targets!=ref]
  out=map(targets,function(x){
    print(paste("Running for",x))
    dat1=dat[,meta[,Gene_col] %in% c(x,ref)]
    meta1=meta[meta[,Gene_col] %in% c(x,ref),]
    meta1[Gene_col]=factor(meta1[,Gene_col],levels=c(ref,x))
    mrk=tryCatch({
      RunNebula(dat1,meta1,form_fixed,sampleCol,cpc=cpc,model=model,method=method)
    }, error=function(cond){
      print("Error!")
      print(cond)
      return(NULL)
    })
    return(mrk)
  })
  names(out)=targets
  out[sapply(out, is.null)]=NULL
  return(out)
  
}

if(!interactive())
{
  args <- commandArgs(trailingOnly=T)
  mat <- open_matrix_dir(dir =args[1])
  meta=qread(args[2])
  mat=mat[,meta$Num==1]
  meta=meta[meta$Num==1,]
  seur@meta.data["Target"]=map_chr(seur@meta.data$Assign,function(x){
    s=strsplit(x,"_")[[1]]
    paste(s[1:length(s)-1],collapse="_")
  })
  
  mrk=RunNebula_pert(seur=seur,form_fixed=~Target+orig.ident,Gene_col="Target",sampleCol="Assign",cpc=.01,model="NBGMM",method="LN",ref="NT",minCells=100,toTest=c())
  qsave(mrk,"results.DE.gene.qs")
}



