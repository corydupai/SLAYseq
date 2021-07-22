library(shiny)
library(shinyFiles)
library(tidyverse)
library(fs)
library(rhandsontable)
library(shinyjs)

source("scripts/functions.R")

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

# Define UI for application that draws a histogram
ui <- 
  # dashboardPage(
  
  # IN SHINY
  # 1) Pick a folder
  # 2) Select fastqs? (Make row)
  # 3) Contrast/.csv file (Row with 3?)
  # 4) Index setup? (Make Row)
  #   A) Fake transcriptome via degenerate codons 
  #   B) Upload transcriptome fasta
  #   C) Transcriptome from fastqs
  # RUN AND WAIT
# 5) Fastqc
# 6) Run Kallisto
# 7) Run DESeq
#
# Further adds
# 1) Just DESeq


navbarPage(
  "Navbar!",
  tabPanel(
    "Complete analysis",
    useShinyjs(),
    mainPanel(
      id = "full_pan",
      fluidRow(
        column(6,
               h4("Choose directory"),
               shinyDirButton("maindir", "Filepath", "Pick project folder")),
        
        hidden(
          div(id = "PH",
              column(6,
                     h4("Choose Fastq files"),
                     shinyFilesButton("fastqs_exist", "Fastq select", "Please select a file", 
                                      multiple = TRUE, viewtype = "detail"))))),
        h2(""),
      rHandsontableOutput("contrasts"),
      hidden(actionButton("act_con", "Create contrast file")),
      hidden(
        div(id = "options_sect",
            
            h2("---------------------------"),
            fluidRow(
              column(4, radioButtons("tr_select",
                                     "Transcriptome options:",
                                     choices = 
                                       c("Use existing" = "Exists",
                                         "Upload my own" = "Upload",
                                         "Degenerate codon sequence (≤ 10M permutations)" = "Deg" ,
                                         "Generate from Fastq files (> 10M permutations)" = "Gen"),
                                     selected = character(0))
                     ,
                     hidden(div(id = "Use Existing", 
                                h2(""),
                                shinyFilesButton("tr_exists", "Transcriptome select", "Please select a file", 
                                                 multiple = FALSE, viewtype = "detail"))),
                     
                     hidden(div(id = "Upload Transcriptome", 
                                h2(""),
                                fileInput("savetr", "Save transcriptome", "Save file as...",
                                          multiple = FALSE))),
                     hidden(div(id = "Degenerate Codons",
                                h2("Degenerate Codons"),
                                textInput("IUPAC",
                                          label = "Add degenerate codon sequence",
                                          value = "NNB", width = NULL,
                                          placeholder = NULL),
                                actionButton("gen_tr", "Add")))),
              column(6, 
                     radioButtons("trim_seqs",
                                  "Trim sequences 3' or 5' sequences from reads?",
                                  choices = 
                                    c("Trim" = "Trim",
                                      "Don't Trim" = "NoTrim"),
                                  selected = character(0)),
                     hidden(textInput("trimmer",
                                      "Sequence(s) to trim (separate with commas)",
                                      placeholder="AAA,GGG")),
                     hidden(actionButton("trimb", "Submit sequence(s)")))))),
      hidden(div(id = "AB",
                 h2("---------------------------"),
                 column(4,offset = 4,
                        actionButton("analyze", "Start analysis")))
      )
    )
  ),
  tabPanel("Differential expression only",
           shinyFilesButton("DE", "File select", "Please select a file to download", multiple = TRUE, viewtype = "detail"),
           tags$div('Make this some text')
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  volumes <- c(Home = "/stor/work/Wilke/cdupai")
  shinyFileChoose(input, "DE", roots = volumes, session = session)
  shinyFileChoose(input, "fastqs_exist", roots = volumes, session = session)
  shinyFileChoose(input, "tr_exists", roots = volumes, session = session)
  shinyDirChoose(input, "maindir", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  # trim_me <- NULL
  trim_me <- eventReactive(c(
    input$trim_seqs,
    input$trimb),{
      
      if(isTruthy(input$trimb)){
        make_adapaters_fa(input$trimmer, path_tr())
        
      } else if(!isTruthy(input$trim_seqs)){
        return(NULL)
      } else if(input$trim_seqs == "Trim"){
        show("trimmer")
        return(NULL)
      } else {
        "Continue"
      }
      
      
    })
  
  observeEvent(c(
    input$tr_select,
    input$trim_seqs,
    input$trimb), {
      req(input$tr_select)
      if(is.null(trim_me())){
        return()
      }
      
      hide("Upload Transcriptome")
      hide("Degenerate Codons")
      hide("Use Existing")
      
      if(input$tr_select =="Exists"){
        show("Use Existing")
      } else if(input$tr_select == "Upload"){
        show("Upload Transcriptome")
      } else if(input$tr_select == "Deg"){
        show("Degenerate Codons")
      } else if(isTruthy(input$trimb)){
        show("AB")
      } else if(input$tr_select == "Gen"){
        show("AB")
      }
    })
  
  DF <- NULL
  DF <- eventReactive( c(
    input$fastqs_exist,
    input$savefastqs
  ), 
  {
    if(isTruthy(input$fastqs_exist)){
      fastq_df <- parseFilePaths(volumes, input$fastqs_exist)

    } else if(isTruthy(input$savefastqs)) {
      fastq_df <- input$savefastqs
      dir.create(paste0(path_tr(),"/fastqs/"),
                 recursive = TRUE)
      file.copy(input$savefastqs$datapath,
                paste0(path_tr(),"/fastqs/",fastq_df$name),
                overwrite = T)
      
    } else {
      return(NULL)
    }
    
    DF <- fastq_df %>%
      mutate(Condition = "Control",
             Read_pair = "1",
             Replicate = "",
             Subfile = "1",
             name = str_replace_all(name,
                                    "\\..*gz",
                                    ""))%>%
             # datapath = paste0(path_tr(),"/fastqs/",name)) %>%
      select(Name = name,
             Condition,
             Replicate,
             Read_pair,
             Subfile,
             Filepath = datapath)
    show("act_con")
    rhandsontable(DF, useTypes = FALSE, stretchH = "all")
  })
  
  output$contrasts <- renderRHandsontable(DF())
  
  out_cons <- eventReactive(input$act_con,{
    dir.create(paste0(path_tr(),"/setup_files/"),
               recursive = TRUE)
    ht <- hot_to_r(input$contrasts) %>%
      # mutate(Filepath = paste0(path_tr(),"/fastqs/",Name)) %>%
      select(Name,
             Filepath,
             Condition,
             Replicate,
             Read_pair,
             Subfile)
    
    write_tsv(ht,
              paste0(path_tr(),"/setup_files/contrasts.tsv"),
              col_names = FALSE)
    paste0(path_tr(),"/setup_files/contrasts.tsv")
    
  })
  observeEvent(out_cons(),{
    show("options_sect")
  })
  
  # observeEvent(input$combine_select,{
  #   show("tr_select")
  #   show("trim_seqs")
  # })
  
  observeEvent(input$trimmer,{
    if(input$trimmer!=""){
      show("trimb") 
    }
  })
  
  transcriptome_file <-
    eventReactive(
      c(input$tr_select,
        input$gen_tr,
        input$savetr,
        input$tr_exists),
      {
        req(input$maindir)
        
        if(input$tr_select =="Deg"){
          req(input$gen_tr)
          seqz <- str_split(input$IUPAC,",") %>%
            unlist() %>%
            as.list()
          seqz_out <- 
            lapply(seqz,
                 degenerate_trasncriptome,
                 path_tr()) %>%
            unlist() %>%
            as.list()
          write.fasta(sequences = seqz_out,
                      names = seqz_out,
                      file.out = paste0(path_tr(),"/setup_files/transcriptome.fasta"))
          # degenerate_trasncriptome(input$IUPAC,
          #                          path_tr())
          f_out <- paste0(path_tr(),"/setup_files/transcriptome.fasta")
        } else if(input$tr_select =="Upload") {
          req(input$savetr)
          file.copy(input$savetr$datapath,
                    paste0(path_tr(),"/setup_files/transcriptome.fasta"),
                    overwrite = T)
          f_out <- paste0(path_tr(),"/setup_files/transcriptome.fasta")
          
        } else if(input$tr_select =="Exists") {
          req(input$tr_exists)
          
          f_out <- parseFilePaths(volumes, input$tr_exists)$datapath %>%
            unlist()
        } else if(input$tr_select == "Gen"){
          f_out <- "GENERATE"
        }
        
        return(f_out)
      })
  
  path_tr <- eventReactive(input$maindir, {
    if(is.numeric(input$maindir)) return(NULL)
    val_out <- paste(c(volumes,
                       input$maindir$path[input$maindir$path != ""] %>% unlist()), 
                     collapse = "/")
    val_out
  })
  
  observeEvent(path_tr(),{
    show("PH")
    showTab(inputId = "fastq_panel", 
            target = "f_select")
  })

  observeEvent(transcriptome_file(),{
    if(length(transcriptome_file()) > 0){
      if(input$tr_select != "Gen"){
        show("AB")
      }
    }
  })
  
  
  observeEvent(input$analyze,{
    req(transcriptome_file)
    # print(transcriptome_file())
    bash_commands <- paste0("bash /stor/home/cdubois/R_analyses/RNAseq_pipeline/RNAseq/scripts/master.sh -c ",
                            out_cons(),
                            " -d ",
                            path_tr(),
                            " -x ",
                            trim_me(),
                            " -t ",
                            transcriptome_file(),
                            " &"
                            )
    system(bash_commands,
           wait = FALSE)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

