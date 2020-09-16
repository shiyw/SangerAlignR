library(shiny)
library(knitr)
library(Biostrings)
library(msa)
library(msaR)
library(sangerseqR)
library(DECIPHER)
source('polypeakfunctions.R')
## 安装msa，msaR及tinytex
# install.packages("msaR")
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("msa")
# install.packages('tinytex')
# tinytex::install_tinytex()
# # to uninstall TinyTeX, run tinytex::uninstall_tinytex() 

examplefile <- "example.ab1"
exampleref <- "TGGCAAGAGAGCGACAGTCAGTCGGACTTACGAGTTGTTTTTACAGGCGCAATTCTTTTTTTAGAATATTATACATTCATCTGGCTTTTTGGATGCACCGATGAGAGATCCAGTTTTCACAGCGAACGCTATGGCTTATCACCCTTTTCACGCGCACAGGCCGGCCGACTTTCCCATGTCAGCTTTCCTTGCGGCGGCTCAACCTTCGTTCTTTCCAGCGCTCACTTTACCACCGGGTCTCAGTAAACCGCTGGCGGATCATGCGCTCTCCGGTGCGGCTGAAGCTGGTTTACACGCGGCGCTTGGACATCACCACCAGGCGGCTCATCTGCGCTCTTTCAAGGGTCTCGAGCCAGAGGAGGATGTTGAGGACGATCCTAAAGTTACATTAGAAGCTAAGGAGCTTTGGGATCAATTCCACAAAATTGGAACAGAAATGGTCATCACTAAATCAGGAAGGTAAGGTCTTTACATTATTTAACCTATTGAATGCTGCATAGGGTGATGTTATTATATTACTCTGCGAAGAGTTGGGTCTATTTTATCGTAAAATATACTTTACATTATAAAATATTGCTCGGTTAAAATTCAGATGTACTGGATGCTGACATAGCATCGAAGCCTCT"

log <- "access.log"

shinyServer(function(input, output, session) {
  inputdata <- reactive({ 
    if(input$example) {
      return(makeBaseCalls(readsangerseq(examplefile), input$ratio))
    } else if(!is.null(input$seq)) {
      return(makeBaseCalls(readsangerseq(input$seq$datapath), input$ratio))
    } else return(NULL)
  })
  refseq <- reactive(
    if(input$example) cleanstring(exampleref)
    else cleanstring(input$ref)
  )
  h <- reactive({
    if (!is.null(input$seq) | input$example) {
      figheight(inputdata(), input$trim5, input$trim3, width=input$x, showtrim=input$showtrim)
    }
  }) 
 
  output$fileUploaded <- reactive({return(!is.null(input$seq) | input$example)})
  output$chromatogram <- renderPlot(
    chromatogram(inputdata(), showcalls="both", trim3=input$trim3, 
                 trim5=input$trim5, width=input$x, showtrim=input$showtrim, showhets=TRUE, cex.base=2), height=reactive(h())
  )
  output$h <- reactive(h())  
  outputdata <- reactive(alignchromatogram(inputdata(), 
                                           trim=input$trimref, 
                                           refseq=refseq(), 
                                           trim5=input$trim5, 
                                           trim3=input$trim3, 
                                           block.width=80
                                           )
                         )
  output$refseq <- renderText(
    if (nchar(input$ref) > 0 | input$example) {outputdata()$refseq})
  output$altseq <- renderText(
    if (nchar(input$ref) > 0 | input$example) {outputdata()$altseq})
  output$alignment <- renderText(
    if (nchar(input$ref) > 0 | input$example) {
      gsub(pattern="^.*#={39}(.+?)#-{39}.*$",
           replacement="\\1",
           x=outputdata()$alignment)
      
  })
  
  output$header <- renderText(
    if (nchar(input$ref) > 0 | input$example) {
      gsub(pattern="(^.+)#\\n#\\n#={39}.+$",
           replacement="\\1",
           x=outputdata()$alignment)
    }
  )
  observe({
    c(input$ratio, input$trim5, input$trim3)
    if(!is.null(input$seq) | input$example == TRUE)
      updateTabsetPanel(session, "maintabset", selected = "Chromatogram")
  })
  
  observe({
    inputref <- input$ref
    inputseq <- input$seq
    if(inputref != "" & !is.null(inputseq)) {
      updateTabsetPanel(session, "maintabset", selected = "Results")
    }
  })
  
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  output$downloadData <- downloadHandler(
    filename = function() paste(input$seq$name, "_report.pdf", sep=""),
    content = function(con) {
       pdfname = knit2pdf(input="peakparser_report.Rnw")
      file.copy(pdfname, con)
    }
  )
  output$msa <- renderMsaR({
    align_output = alignchromatogram(inputdata(), 
                                           trim=input$trimref, 
                                           refseq=refseq(), 
                                           trim5=30, 
                                           trim3=30, 
                                           block.width=80
                                           )
    seq1 <- align_output$altseq
    seq2 <- align_output$refseq
    seq3 <- refseq()
    seqs <- c(seq1, seq2, seq3)
    names(seqs) <- c("alle1", "alle2", "ref")
    seq_alignment <- msa(DNAStringSet(seqs))
    library(stringr)
    left_cut_n <- median(str_length(str_extract(seq_alignment@unmasked, "^-*")))
    right_cut_n <- median(str_length(str_extract(seq_alignment@unmasked, "-*$")))
    seq_alignment@unmasked[[1]] <- str_replace(seq_alignment@unmasked[[1]], paste0("^.{",left_cut_n,"}"), "")
    seq_alignment@unmasked[[2]] <- str_replace(seq_alignment@unmasked[[2]], paste0("^.{",left_cut_n,"}"), "")
    seq_alignment@unmasked[[3]] <- str_replace(seq_alignment@unmasked[[3]], paste0("^.{",left_cut_n,"}"), "")
    seq_alignment@unmasked[[1]] <- str_replace(seq_alignment@unmasked[[1]], paste0(".{",right_cut_n,"}$"), "")
    seq_alignment@unmasked[[2]] <- str_replace(seq_alignment@unmasked[[2]], paste0(".{",right_cut_n,"}$"), "")
    seq_alignment@unmasked[[3]] <- str_replace(seq_alignment@unmasked[[3]], paste0(".{",right_cut_n,"}$"), "")
    # msaPrettyPrint(seq_alignment, output="pdf", showNames="none", showLogo="top",
    #                logoColors="rasmol", shadingMode="similar",
    #                showLegend=FALSE, askForOverwrite=FALSE)
    writeXStringSet(as(unmasked(seq_alignment), "XStringSet"), file="aln.fasta")
    BrowseSeqs(seq_alignment@unmasked,
               htmlFile = paste("www/alignment.html", sep = ""),
               openURL = FALSE,
               colorPatterns = TRUE,
               highlight = NA,
               colWidth = Inf)
    msaR("aln.fasta", seqlogo=TRUE,conservation=TRUE)
  }
)
  output$downloadData2 <- downloadHandler(
    filename = function() paste(input$seq$name, "_aln.fasta", sep=""),
    content = function(con) {
      alignment_file_name = "aln.fasta"
      file.copy(alignment_file_name, con)
    }
  )
  #logging
  observe({
    if(!is.null(inputdata())) {
      isolate({
        alog <- file(log, "a")
        cat(format(Sys.time(), "%Y-%b-%d_%H-%M-%S"), file=alog)
        cat("\t", file=alog)
        if(input$example == TRUE) cat("Opened Example", file=alog)
        else cat(input$seq$name, file=alog)
        cat("\n", file=alog)
        close(alog)
      })
    }
  })
  output$ssRversion <- renderText(as.character(packageVersion("sangerseqR")))
})