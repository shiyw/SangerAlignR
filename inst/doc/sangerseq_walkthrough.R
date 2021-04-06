## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
  BiocStyle::latex()

## ----knitr, eval=TRUE, include=FALSE---------------------------------------
  library(knitr, quietly=TRUE)
  library(sangerseqR, quietly=TRUE)
  library(Biostrings, quietly=TRUE)
  opts_chunk$set(tidy=TRUE)

#modified from: 
#https://github.com/yihui/knitr-examples/blob/master/077-wrap-output.md
hook_output = knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
        n <- 90
        x <- knitr:::split_lines(x)
        # any lines wider than n should be wrapped
        if (any(nchar(x) > n)) {
           x <- gsub(sprintf('(.{%d})', n), "\\1\n## ", x)
        }
        
    hook_output(x, options)
})

## --------------------------------------------------------------------------
hetab1 <- read.abif(system.file("extdata", 
                                "heterozygous.ab1", 
                                package="sangerseqR"))
str(hetab1, list.len=20)

## --------------------------------------------------------------------------
homoscf <- read.scf(system.file("extdata", 
                                "homozygous.scf", 
                                package="sangerseqR"))
str(homoscf)

## --------------------------------------------------------------------------
#from a sequence file object
homosangerseq <- sangerseq(homoscf)

#directly from the file
hetsangerseq <- readsangerseq(system.file("extdata", 
                                          "heterozygous.ab1", 
                                          package="sangerseqR"))
str(hetsangerseq)

## --------------------------------------------------------------------------
#default is to return a DNAString object
Seq1 <- primarySeq(homosangerseq)
reverseComplement(Seq1)

#can return as string
primarySeq(homosangerseq, string=TRUE)

## --------------------------------------------------------------------------
chromatogram(hetsangerseq, width=200, height=2, trim5=50, trim3=100, 
             showcalls='both', filename="chromatogram.pdf")

## --------------------------------------------------------------------------
hetcalls <- makeBaseCalls(hetsangerseq, ratio=0.33)
hetcalls

## --------------------------------------------------------------------------
chromatogram(hetcalls, width=100, height=2, trim5=50, trim3=100, 
             showcalls='both', filename="chromatogram2.pdf")

## --------------------------------------------------------------------------
ref <- subseq(primarySeq(homosangerseq, string=TRUE), start=30, width=500)
hetseqalleles <- setAllelePhase(hetcalls, ref, trim5=50, trim3=300)
hetseqalleles

## --------------------------------------------------------------------------
pa <- pairwiseAlignment(primarySeq(hetseqalleles)[1:400], 
                        secondarySeq(hetseqalleles)[1:400], 
                        type="global-local")
writePairwiseAlignments(pa)

