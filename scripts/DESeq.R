require(data.table)
require(tximport)
require(rhdf5)
require(DESeq2)
require(parallel)
require(optparse)
require(tidyverse)

deseq_run <- function(basedir,transcriptome){

tx_im <- fread(transcriptome,
               col.names = c("filepath",
                             "replicate",
                             "condition"),
               header = FALSE) %>%
  mutate(
    condition = as.factor(condition),
    condition = fct_relevel(condition, "Control"))

tx_obj <- tximport(files=tx_im$filepath,type="kallisto",txOut=TRUE) # We use txOut = true because our transcripts are genes. This isn't true for eukaryotes
rownames(tx_obj$counts) <- as.character(rownames(tx_obj$counts)) # Idk if you need to do this but you can

coldata <- data.frame(
  condition=tx_im$condition,
  replicate=tx_im$replicate)

rownames(coldata) = paste0(tx_im$condition,tx_im$rep)
colnames(tx_obj$counts) = paste0(tx_im$condition,tx_im$rep)
tx_obj_orig_counts <- as.data.frame(tx_obj$counts)
colnames(tx_obj_orig_counts) <- paste0(colnames(tx_obj_orig_counts), "_RAW")
tx_obj$counts <- tx_obj$counts + 1


dds <- DESeqDataSetFromTximport(tx_obj,coldata, ~condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=FALSE) 
print(resultsNames(dds))
counter = 1
for(rep in unique(coldata$condition)){
  if(rep != "Control"){
    print(rep)
    condz <- paste0("condition_",rep,"_vs_Control")
    res_temp <- results(dds, name = condz)  %>%
      as.data.frame() %>%
      mutate(condition = rep) %>%
      cbind(counts(dds,normalized = TRUE)%>%
              as.data.frame(), tx_obj_orig_counts) %>%
      rownames_to_column(var = "sequence")

    if(counter == 1){
      res <- res_temp
      counter = 2
    } else {
      res <- rbind(res, res_temp)
    }
    
  }
}

return(res)

}


option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name", metavar="character"),
  make_option(c("-c", "--contrast"), type="character", default=NULL, 
              help="contrast file name", metavar="character")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Need input directory", call.=FALSE)
} else if (is.null(opt$contrast)){
  print_help(opt_parser)
  stop("Need input directory", call.=FALSE)
}



