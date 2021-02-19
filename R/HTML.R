#' HTML_HandleEndOfLine add a end of a line
#'
#' @param EndOfTheLine At what position you are at the end of this line
#' @param LineWidth The length of the line
#'
#' @return charater chain with
#' @export
#'
#' @examples
#' HTML_HandleEndOfLine(EndOfTheLine=277,LineWidth=60)
HTML_HandleEndOfLine <- function(EndOfTheLine, LineWidth) {
  if ((EndOfTheLine %% LineWidth) == 0) {
    BeginingOfTheLine <- EndOfTheLine - LineWidth
  } else{
    ##Case where the modulo is lesser than 10
    # 966
    #
    if ((EndOfTheLine %% LineWidth) < 10) {
      return("")
    }
    BeginingOfTheLine <- EndOfTheLine - EndOfTheLine %% LineWidth
  }
  ThicksValues <- seq(BeginingOfTheLine + 10, EndOfTheLine, 10)
  NbBlocks <- length(ThicksValues)
  OneChunk <- paste0(c(rep("&nbsp;", 9), "|"), collapse = "")
  LineChunks <- paste0(rep(OneChunk, NbBlocks), collapse = "&nbsp;")
  LineValues <- formatC(x = ThicksValues, width = 10, flag = " ")
  LineValues <- gsub(" ", "&nbsp;", LineValues)
  LineValues <- paste(LineValues, collapse = "&nbsp;")
  ReturnLine <-
    paste(
      "</br><mark class='position'>",
      LineChunks,
      "</br>",
      LineValues,
      "</br></br></mark>",
      sep = ""
    )
  return(ReturnLine)
}


#' UpdateSequence Transform a AA sequence to a HTML ready to print string
#'
#' @param IR Integer Ranges with all necessay informations
#' @param PROTEIN_AA AA sequence
#'
#'
#' @export
#'
#' @return a string with span tag
#'
UpdateSequence<-function(IR,PROTEIN_AA){
  INFOS<-data.frame(IR)
  PROTEIN_AA[INFOS$start]<-paste(INFOS$Start.Prefix, PROTEIN_AA[INFOS$start], INFOS$Start.Suffix, sep="")
  PROTEIN_AA[INFOS$end]<-paste(INFOS$End.Prefix, PROTEIN_AA[INFOS$end], INFOS$End.Suffix, sep="")
  return(PROTEIN_AA)
}


#' HTML_Pretty_Sequence Prepare an HTML ready
#'
#' @param Sequence String of the protein AA Sequence
#' @param Peptides data.frame of the peptides attributed to this protein
#' @param Rescue Not Used for now
#' @param PerfectIsoformName Name of the Isoform as observed in the
#' @param LineWidth Number of AA before a new line (will not include the space inside)
#'
#' @return HTML string of the sequence embedded with features included as span
#' @import IRanges
#' @import GenomicRanges
#' @export
#'
#' @examples
#' library(Biostrings)
#' PB_FilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' ProteinBank<-readAAStringSet(
#'   PB_FilePath
#' )
#' Sequence <- as.vector(ProteinBank["Example"])
#' Peptides <- data.frame(
#'  Start = c(2, 2, 240),
#'  End = c(10, 11, 253),
#'  Pep = c(1.9437e-06, 0.00014334, 0.0001654),
#'  Proteins = c(
#'   "strngt.11487.1 Q8WUM4-3 strngt.11487.3 strngt.11487.2",
#'   "strngt.11487.1 Q8WUM4-3 strngt.11487.3 strngt.11487.2",
#'   "strngt.11487.2"),
#'  Sequence = c("ATFISVQLK", "ATFISVQLKK", "YFYFQEVFPVLAAK ")
#' )
#' Sequence <- HTML_Pretty_Sequence(Sequence,Peptides)
HTML_Pretty_Sequence <- function(Sequence,
                                 Peptides,
                                 Rescue = NULL,
                                 PerfectIsoformName = NULL,
                                 LineWidth = 60) {
  #Protein IntegerRange
  LocalProteinSequence <-
    unlist(strsplit(x = as.vector(Sequence), split = ""))
  LocalProteinLength <- length(LocalProteinSequence)
  Protein.IR <- IRanges::IRanges(start = 1, end = LocalProteinLength)
  GenomicRanges::mcols(Protein.IR) <- data.frame(
    Class = "Protein",
    Start.Prefix = "<div style=\"font-family:monospace;font-size:16px;height:100%;overflow-y:scroll;\">",
    Start.Suffix = "",
    End.Prefix = "",
    End.Suffix = paste(
      HTML_HandleEndOfLine(LocalProteinLength,
                           LineWidth),
      "</div>"
    )
  )

  if (!is.null(PerfectIsoformName)) {
    UNIPROTIsoformName <-
      gsub(".*\\|([^\\|]*)-([^\\|]+)\\|.*",
           "\\1-\\2",
           PerfectIsoformName)
    #SplicingVariations IntegerRange
    SplicingVariations.IR <- GetRangesFromUniprot(UNIPROTIsoformName)
    ##Case Where the original sequence ends after the isoform sequence
    # Q69JNO-3 VSP_058760
    if (sum(!SplicingVariations.IR %within% Protein.IR)) {
      end(SplicingVariations.IR)[end(SplicingVariations.IR) > LocalProteinLength] <-
        LocalProteinLength
    }
  } else{
    SplicingVariations.IR <- GetEmptySpanRanges()
  }

  #Newline IntegerRange
  if (LocalProteinLength > LineWidth) {
    NewLines <- seq(LineWidth, LocalProteinLength, LineWidth)
    NewLines <- NewLines[!NewLines %in% LocalProteinLength]
    NewLines.AfterChar <-
      unlist(lapply(NewLines, function(x)
        HTML_HandleEndOfLine(x, LineWidth)))
    NewLines.IR <- IRanges::IRanges(start = NewLines, end = NewLines)
    GenomicRanges::mcols(NewLines.IR) <- data.frame(
      Class = "Newline",
      Start.Prefix = "",
      Start.Suffix = "",
      End.Prefix = "",
      End.Suffix = NewLines.AfterChar
    )
  } else{
    NewLines.IR <- IRanges::IRanges()
    NewLines <- c()
  }

  #Block IntegerRange
  BlocksPosition <- seq(10, LocalProteinLength, 10)
  BlocksPosition <- BlocksPosition[!BlocksPosition %in% NewLines]
  BlocksPosition.IR <-
    IRanges::IRanges(start = BlocksPosition,
                     end = BlocksPosition)
  GenomicRanges::mcols(BlocksPosition.IR) <- data.frame(
    Class = "BlockEnd",
    Start.Prefix = "",
    Start.Suffix = "",
    End.Prefix = "",
    End.Suffix = "&nbsp;"
  )



  IDX_Shared <- grep(" ", Peptides$Proteins)
  IDX_Proteotypic <- grep(" ", Peptides$Proteins, invert = TRUE)

  #ProteotypicIntergerRange
  SharedPeptides.IR <- IRanges::IRanges(start = Peptides$Start[IDX_Shared],
                                        end = Peptides$End[IDX_Shared])
  GenomicRanges::mcols(SharedPeptides.IR) <- data.frame(
    Title = paste(
      "Peptide: ",
      Peptides$Sequence[IDX_Shared],
      "(",
      formatC(Peptides$Pep[IDX_Shared], digits = 2),
      ")\n",
      sep = ""
    ),
    Class = "Peptide"
  )

  #PeptidesIntergerRange
  Proteotypic.IR <- IRanges::IRanges(start = Peptides$Start[IDX_Proteotypic],
                                     end = Peptides$End[IDX_Proteotypic])
  GenomicRanges::mcols(Proteotypic.IR) <-
    data.frame(
      Title = paste(
        "Peptide: ",
        Peptides$Sequence[IDX_Proteotypic],
        "(",
        formatC(Peptides$Pep[IDX_Proteotypic], digits = 2),
        ")\n",
        sep = ""
      ),
      Class = "Peptide Proteotypic"
    )

  Peptides.IR <- c(SharedPeptides.IR, Proteotypic.IR)

  #Infos
  SpanRanges <- c(SplicingVariations.IR, Peptides.IR)
  SpanRanges <- CollapseSpanRanges(SpanRanges)
  DisjoinSpanRange <- SplitSpanRanges(SpanRanges, NewLines.IR)
  GenomicRanges::mcols(DisjoinSpanRange)$Start.Prefix <-
    paste(
      "<span class=\"",
      GenomicRanges::mcols(DisjoinSpanRange)$Class,
      "\" title=\"",
      GenomicRanges::mcols(DisjoinSpanRange)$Title,
      "\">",
      sep = ""
    )
  GenomicRanges::mcols(DisjoinSpanRange)$Start.Suffix <- ""
  GenomicRanges::mcols(DisjoinSpanRange)$End.Prefix <- ""
  GenomicRanges::mcols(DisjoinSpanRange)$End.Suffix <- "</span>"

  #Insert Newlines into Disjoin IntegerRange
  NewLines.IR.ToHandle <- NewLines.IR
  OverlapsWithNewline <-
    data.frame(findOverlaps(DisjoinSpanRange, NewLines.IR.ToHandle))
  while (dim(OverlapsWithNewline)[1] > 0) {
    ##Case where the a peptide end with a new line it created a span that did not close anywhere and in doing so make it all
    if (end(DisjoinSpanRange[OverlapsWithNewline$queryHits[1]]) == end(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]])) {

    } else{
      Before <- DisjoinSpanRange[OverlapsWithNewline$queryHits[1]]
      After <- Before
      end(Before) <-
        start(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]])
      start(After) <-
        end(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]]) + 1
      DisjoinSpanRange <-
        DisjoinSpanRange[-OverlapsWithNewline$queryHits[1]]
      DisjoinSpanRange <- c(DisjoinSpanRange, Before, After)
      DisjoinSpanRange <- sort(DisjoinSpanRange)
    }
    NewLines.IR.ToHandle <-
      NewLines.IR.ToHandle[-OverlapsWithNewline$subjectHits[1]]
    OverlapsWithNewline <-
      data.frame(findOverlaps(DisjoinSpanRange, NewLines.IR.ToHandle))
  }

  # mcols(DisjoinSpanRange)$Start.Prefix <- paste( "<span class=\"", mcols(DisjoinSpanRange)$Class,
  #                                                "\" title=\"", mcols(DisjoinSpanRange)$Title, "\">", sep = "")
  GenomicRanges::mcols(DisjoinSpanRange)$Start.Prefix <-
    paste(
      "<span class=\"",
      GenomicRanges::mcols(DisjoinSpanRange)$Class,
      "\" data-toggle=\"popover\" data-html=\"true\" data-trigger=\"hover click\" data-placement=\"bottom\" title=\"Features\" data-content=\"",
      GenomicRanges::mcols(DisjoinSpanRange)$Title,
      "\">",
      sep = ""
    )

  GenomicRanges::mcols(DisjoinSpanRange)$Start.Suffix <- ""
  GenomicRanges::mcols(DisjoinSpanRange)$End.Prefix <- ""
  GenomicRanges::mcols(DisjoinSpanRange)$End.Suffix <- "</span>"

  LocalProteinSequence <-
    UpdateSequence(BlocksPosition.IR, LocalProteinSequence)
  LocalProteinSequence <-
    UpdateSequence(DisjoinSpanRange, LocalProteinSequence)
  LocalProteinSequence <-
    UpdateSequence(NewLines.IR, LocalProteinSequence)
  LocalProteinSequence <-
    UpdateSequence(Protein.IR, LocalProteinSequence)
  LocalProteinSequence <- paste(LocalProteinSequence, collapse = "")
  return(LocalProteinSequence)
}


#' HTML_PrettyBlastSequence prepare and HTML ready string with some more informations
#'
#' @param Sequence String of the protein AA Sequence
#' @param Peptides data.frame of the peptides attributed to this protein
#' @param Rescue Not Used for now
#' @param LineWidth Number of AA before a new line (will not include the space inside)
#'
#' @return HTML string of the sequence embedded with features included as span
#' @import IRanges
#' @import GenomicRanges
#' @export
#'
#' @examples
#' library(Biostrings)
#' PB_FilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' ProteinBank<-readAAStringSet(
#'   PB_FilePath
#' )
#' Sequence <- as.vector(ProteinBank["Example"])
#' Peptides <- data.frame(
#'  Start = c(2, 2, 240),
#'  End = c(10, 11, 253),
#'  Pep = c(1.9437e-06, 0.00014334, 0.0001654),
#'  Proteins = c(
#'   "strngt.11487.1 Q8WUM4-3 strngt.11487.3 strngt.11487.2",
#'   "strngt.11487.1 Q8WUM4-3 strngt.11487.3 strngt.11487.2",
#'   "strngt.11487.2"),
#'  Sequence = c("ATFISVQLK", "ATFISVQLKK", "YFYFQEVFPVLAAK ")
#' )
#' Sequence <- HTML_Pretty_Sequence(Sequence,Peptides)
HTML_PrettyBlastSequence<-function(
  Sequence,
  Peptides,
  Rescue,
  LineWidth=60){

  #Protein IntegerRange
  LocalProteinSequence <- unlist(strsplit(x = as.vector(Sequence), split = ""))
  LocalProteinLength<-length(LocalProteinSequence)
  Protein.IR <- IRanges(start = 1, end = LocalProteinLength)
  mcols(Protein.IR) <- data.frame(
    Class = "Protein",
    Start.Prefix = "<div style=\"font-family:monospace;font-size:16px;height:100%;overflow-y:scroll;\">",
    Start.Suffix = "",
    End.Prefix = "",
    End.Suffix = paste(
      HTML_HandleEndOfLine(LocalProteinLength, LineWidth),
      "</div>"
    )
  )

  #Newline IntegerRange
  if(LocalProteinLength>LineWidth){
    NewLines<-seq(LineWidth,LocalProteinLength,LineWidth)
    NewLines<-NewLines[!NewLines %in% LocalProteinLength]
    NewLines.AfterChar<-unlist( lapply( NewLines, function(x) HTML_HandleEndOfLine(x,LineWidth)))
    NewLines.IR<-IRanges( start=NewLines, end=NewLines)
    mcols(NewLines.IR)<-data.frame(Class="Newline", Start.Prefix="",
                                   Start.Suffix="", End.Prefix="", End.Suffix=NewLines.AfterChar)
  }else{
    NewLines.IR<-IRanges()
    NewLines<-c()
  }

  #Block IntegerRange
  BlocksPosition<-seq(10, LocalProteinLength, 10)
  BlocksPosition<-BlocksPosition[!BlocksPosition %in% NewLines]
  BlocksPosition.IR <- IRanges(start = BlocksPosition, end = BlocksPosition)
  mcols(BlocksPosition.IR) <- data.frame(
    Class = "BlockEnd",
    Start.Prefix = "",
    Start.Suffix = "",
    End.Prefix = "",
    End.Suffix = "&nbsp;"
  )

  #Blast differences
  Blast.IR <- IRanges(start = Rescue$QueryStart,
                      end = Rescue$QueryEnd)
  mcols(Blast.IR) <- data.frame(
    Class = "Blast",
    Title = paste(
      Rescue$SubjectSeqId,
      " from:",
      Rescue$SubjectStart,
      " to:",
      Rescue$SubjectEnd,
      " of: ",
      Rescue$SubjectLength,
      sep = ""
    )
  )
  Diffs.IR<-CompareTwoSequences(as.vector(Rescue$QuerySequence),
                                as.vector(Rescue$SubjectSequence),
                                Rescue$QueryStart)
  Blast.IR<-c(Blast.IR, Diffs.IR)


  IDX_Shared <- grep(" ",Peptides$Proteins)
  IDX_Proteotypic <- grep(" ",Peptides$Proteins,invert = TRUE)

  #ProteotypicIntergerRange
  SharedPeptides.IR <- IRanges(start = Peptides$Start[IDX_Shared],
                               end = Peptides$End[IDX_Shared])
  mcols(SharedPeptides.IR) <- data.frame(
    Title = paste(
      "Peptide: ",
      Peptides$Sequence[IDX_Shared],
      "(",
      formatC(Peptides$Pep[IDX_Shared], digits = 2),
      ")\n",
      sep = ""
    ),
    Class = "Peptide"
  )

  #PeptidesIntergerRange
  Proteotypic.IR <- IRanges(start = Peptides$Start[IDX_Proteotypic],
                            end = Peptides$End[IDX_Proteotypic])
  mcols(Proteotypic.IR) <- data.frame(
    Title = paste(
      "Peptide: ",
      Peptides$Sequence[IDX_Proteotypic],
      "(",
      formatC(Peptides$Pep[IDX_Proteotypic], digits = 2),
      ")\n",
      sep = ""
    ),
    Class = "Peptide Proteotypic"
  )

  Peptides.IR<-c(SharedPeptides.IR,Proteotypic.IR)

  #Infos
  SpanRanges<-c(Blast.IR,Peptides.IR)
  SpanRanges<-CollapseSpanRanges(SpanRanges)
  DisjoinSpanRange<-SplitSpanRanges(SpanRanges,NewLines.IR)
  mcols(DisjoinSpanRange)$Start.Prefix <- paste("<span class=\"", mcols(DisjoinSpanRange)$Class, "\" title=\"",
                                                mcols(DisjoinSpanRange)$Title, "\">", sep = "")
  mcols(DisjoinSpanRange)$Start.Suffix <- ""
  mcols(DisjoinSpanRange)$End.Prefix <- ""
  mcols(DisjoinSpanRange)$End.Suffix <- "</span>"

  #Insert Newlines into Disjoin IntegerRange
  NewLines.IR.ToHandle<-NewLines.IR
  OverlapsWithNewline<-data.frame(findOverlaps(DisjoinSpanRange,NewLines.IR.ToHandle))
  while(dim(OverlapsWithNewline)[1]>0){
    ##Case where the a peptide end with a new line it created a span that did not close anywhere and in doing so make it all
    if (end(DisjoinSpanRange[OverlapsWithNewline$queryHits[1]])==end(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]])) {
    } else{
      Before <- DisjoinSpanRange[OverlapsWithNewline$queryHits[1]]
      After <- Before
      end(Before) <- start(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]])
      start(After) <- end(NewLines.IR.ToHandle[OverlapsWithNewline$subjectHits[1]]) + 1
      DisjoinSpanRange <- DisjoinSpanRange[-OverlapsWithNewline$queryHits[1]]
      DisjoinSpanRange <- c(DisjoinSpanRange, Before, After)
      DisjoinSpanRange <- sort(DisjoinSpanRange)
    }
    NewLines.IR.ToHandle<-NewLines.IR.ToHandle[-OverlapsWithNewline$subjectHits[1]]
    OverlapsWithNewline<-data.frame(findOverlaps(DisjoinSpanRange,NewLines.IR.ToHandle))
  }

  # mcols(DisjoinSpanRange)$Start.Prefix <- paste( "<span class=\"", mcols(DisjoinSpanRange)$Class,
  #                                                "\" title=\"", mcols(DisjoinSpanRange)$Title, "\">", sep = "")
  mcols(DisjoinSpanRange)$Start.Prefix <- paste("<span class=\"", mcols(DisjoinSpanRange)$Class,
                                                "\" data-toggle=\"popover\" data-html=\"true\" data-trigger=\"hover click\" data-placement=\"bottom\" title=\"Features\" data-content=\"", mcols(DisjoinSpanRange)$Title, "\">", sep = "")

  mcols(DisjoinSpanRange)$Start.Suffix <- ""
  mcols(DisjoinSpanRange)$End.Prefix <- ""
  mcols(DisjoinSpanRange)$End.Suffix <- "</span>"

  LocalProteinSequence<-UpdateSequence(BlocksPosition.IR, LocalProteinSequence)
  LocalProteinSequence<-UpdateSequence(DisjoinSpanRange, LocalProteinSequence)
  LocalProteinSequence<-UpdateSequence(NewLines.IR, LocalProteinSequence)
  LocalProteinSequence<-UpdateSequence(Protein.IR, LocalProteinSequence)
  LocalProteinSequence<-paste(LocalProteinSequence, collapse="")

  #Status
  PeptideInSAAV <- FALSE
  PeptideInInsertion <- FALSE
  PeptideInDeletion <- FALSE
  PeptideOutOfBlastSAAV <- FALSE
  SAAV <- FALSE
  Insertion <- FALSE
  Deletion <- FALSE
  if(sum(mcols(Blast.IR)$Class=="SAAV")>0){SAAV <- TRUE}
  if(sum(mcols(Blast.IR)$Class=="Insertion")>0){Insertion <- TRUE}
  if(sum(mcols(Blast.IR)$Class=="Deletion")>0){Deletion <- TRUE}
  if(sum(Peptides.IR %over% Blast.IR[mcols(Blast.IR)$Class=="SAAV"])>0){ PeptideInSAAV <- TRUE}
  if(sum(Peptides.IR %over% Blast.IR[mcols(Blast.IR)$Class=="Insertion"])>0){ PeptideInInsertion <- TRUE}
  if(sum(Peptides.IR %over% Blast.IR[mcols(Blast.IR)$Class=="Deletion"])>0){ PeptideInDeletion <- TRUE}
  if(sum(!Peptides.IR %over% Blast.IR)>0){ PeptideOutOfBlastSAAV <- TRUE}
  return(list(SequenceInHTML=LocalProteinSequence,
              SAAVxPep=PeptideInSAAV,
              SAAV=SAAV,
              InsertionxPep=PeptideInInsertion,
              Insertion=Insertion,
              DeletionxPep=PeptideInDeletion,
              Deletion=Deletion,
              OutOfBlast=PeptideOutOfBlastSAAV))
}

#' Alter a peptide sequence to have an harmonious presentation
#'
#' @param Sequences Vector of peptides Sequence
#' @param MaxNChar Maximum width size of the peptide
#' @return string vector as a HTML span to have shorter form
#' @export
#'
#' @examples
#' library("DT")
#' library("data.table")
#' P_FilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' Peptides<-fread(
#'  P_FilePath,
#'  integer64="double")
#' Peptides<-Peptides[,c("Sequence",
#' "Start position","End position")]
#' datatable(Peptides)
#' Peptides$Sequence<-HTML_VECTOR_Peptide2Popover(Peptides$Sequence)
#' datatable(Peptides,escape=FALSE,
#' options=list(
#'  drawCallback=JS(
#'  paste("$('[data-toggle=\"popover\"]').popover({container: 'body'});")
#'  )))
HTML_VECTOR_Peptide2Popover<-function(Sequences,MaxNChar=15){
  Too_Long_Sequence_Lengths <- nchar(Sequences) > 15
  Sequences[Too_Long_Sequence_Lengths] <- paste(
    "<span  data-toggle=\"popover\"
  data-trigger=\"hover click\" data-placement=\"right\"
  data-html=\"true\" data-delay={show: 500, hide: 100}
  title=\"Proteins\" data-content=\"",
    Sequences[Too_Long_Sequence_Lengths],"\">",
    base::substr(Sequences[Too_Long_Sequence_Lengths],
                 1,MaxNChar-3),
    "...</span>",sep="")
  return(Sequences)
}
