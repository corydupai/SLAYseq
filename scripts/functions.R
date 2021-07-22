library(seqinr)
library(tidyverse)

degenerate_trasncriptome <- function(deg_sequence,
                                     experiment_dir){
  
  deg_sequence <- toupper(deg_sequence)
  dna_dict <- c(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "U" = "U",
    "W" = "AT",
    "S" = "CG",
    "M" = "AC",
    "K" = "GT",
    "R" = "AG",
    "Y" = "CT",
    "B" = "CGT",
    "D" = "AGT",
    "H" = "ACT",
    "V" = "ACG",
    "N" = "ACGT"
  )
  
  deg_sequence <- str_split(deg_sequence,
                            pattern = "") %>%
    unlist()
  
  deg_sequence <- lapply(deg_sequence,
                        function(x,dna_dict){
                          dna_dict[x == names(dna_dict)] %>%
                            str_split("") %>%
                            unlist
                        }, dna_dict = dna_dict )%>%
    expand.grid() %>%
    unite("Pasted",
          sep = "") %>%
    unique() %>%
    unlist() %>%
    unname()
  return(deg_sequence)

  
}


make_adapaters_fa <- function(trim_sequence,
                                     experiment_dir){
  
  rem_sequence <- str_split(trim_sequence,
                            pattern = ",") %>%
    unlist() %>%
    str_remove_all("\\s")%>%
    toupper()
  
  write.fasta(sequences = as.list(rem_sequence),
              names = rem_sequence,
              file.out = paste0(experiment_dir,"/setup_files/remove_sequences.fasta"))
  return(paste0(experiment_dir,"/setup_files/remove_sequences.fasta"))
  
}