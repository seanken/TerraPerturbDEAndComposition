version 1.0

workflow GetDEandComposition {
  input {
    String inputDir="gs://fc-secure-b42fb9b0-04ed-4260-9c28-aa1274233114/Neur/toShare_BP" ##qs file with Seurat object and 
    File metaQS="gs://fc-secure-b42fb9b0-04ed-4260-9c28-aa1274233114/Neur/meta.human.qs"
    String mode="DE"
  }

  if(mode=="DE") {
    call GetDEByGuide {
      input:
        input_Dir= inputDir,
        metaQS=metaQS,
    }
  }

  if(mode=="DEGene") {
    call GetDEByGene {
      input:
        input_Dir= inputDir,
        metaQS=metaQS,
    }
  }

  if(mode=="Composition") {
    call GetComposition {
      input:
          metaQS=metaQS,
    }
  }

  output{
    File? de_by_guide = GetDEByGuide.de_by_guide
    File? de_by_gene = GetDEByGene.de_by_gene
    File? composition = GetComposition.composition
  }

}


task GetDEByGuide {
  input {
    String input_Dir
    File metaQS
  }

#gsutil cp gs://fc-secure-b42fb9b0-04ed-4260-9c28-aa1274233114/RScripts/GetDE.guide.R /app/GetDE.guide.R
  command <<<
    gsutil cp ~{input_Dir} dirBP
    Rscript /app/GetDE.guide.R dirBP ~{metaQS}
  >>>

  output {
    File de_by_guide = "results.DE.qs"
  }

  runtime{
        docker: "seanken/perturbr:latest"
        zones: "us-central1-b"
        memory: "80G"
        disks: "local-disk 100 HDD"
  }
}

task GetDEByGene {
  input {
    String input_Dir
    File metaQS
  }
    #gsutil cp gs://fc-secure-b42fb9b0-04ed-4260-9c28-aa1274233114/RScripts/GetDE.Gene.R /app/GetDE.Gene.R
  command <<<
    gsutil cp ~{input_Dir} dirBP
    Rscript /app/GetDE.Gene.R dirBP ~{metaQS}
  >>>

  output {
    File de_by_gene = "results.DE.gene.qs"
  }

  runtime{
        docker: "seanken/perturbr:latest"
        zones: "us-central1-b"
        memory: "80G"
        disks: "local-disk 100 HDD"
  }
}

task GetComposition {
  input {
    File metaQS
  }

    #gsutil cp gs://fc-secure-b42fb9b0-04ed-4260-9c28-aa1274233114/RScripts/TestClusterProportions.R /app/TestClusterProportions.R
  command <<<
    Rscript /app/TestClusterProportions.R ~{metaQS}
  >>>

  output {
    File composition = "results.cluster.comp.qs"
  }

  runtime{
        docker: "seanken/perturbr:latest"
        zones: "us-central1-b"
        memory: "80G"
        disks: "local-disk 100 HDD"
  }

}