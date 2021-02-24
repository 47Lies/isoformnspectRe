#' Handle fold change
#' Will substitute 15 to a more extreme value
#'
#' @param FoldChanges vector of fold changes (numeric vector)
#'
#' @return vector of numeric values with absolute values set to 15
#' @export
#'
#' @examples
#' FC<-c(+Inf,9,3,-20,-Inf)
#' HandleFoldChange(FC)
HandleFoldChange <- function(FoldChanges) {
  FoldChanges[FoldChanges < -15] <- -15
  FoldChanges[FoldChanges > 15] <- +15
  return(FoldChanges)
}

#' substitute the peptide by a html span tag that will substitute a too long peptide
#'
#' @param Peptide vector of peptides
#'
#' @return vector of peptides to be used in a DT::data.table will need escape=TRUE option
#' @export
#'
#' @examples
#' peptides<-c("FEZQFAZRFSQJDFCGFSDGHREZQGVSVFDSVGFDS",
#' "RETGERVGDF")
#' TrimSequenceOutput(peptides)
TrimSequenceOutput <- function(Peptide) {
  NeoPeptides <- Peptide
  TooLong <- nchar(Peptide) > 20
  NeoPeptides[TooLong] <- paste("<span title='",
                                Peptide[TooLong],
                                "'>",
                                substr(Peptide[TooLong], 1, 17),
                                "...</Span>"
                                ,
                                sep = "")
  return(NeoPeptides)
}


#' Title
#'
#' @param LIST_SampleNames list of sample names
#' @param LIST_SparlinesOptions options for sparkline
#'
#' @return List of options for data table including sparklines
#'
#' @examples
#' SamplesNames <- list(
#' group_un=c("echantillon un","echantillon deux"),
#' group_deux=c("echantillon trois","echantillon quatre"))
LocalOptions <- function(LIST_SampleNames,
                         LIST_SparlinesOptions = list(
                           type = "bar",
                           height = 40,
                           width = 60,
                           highlightColor = "red",
                           chartRangeMin = 0,
                           chartRangeMax = 1,
                           tooltipFormat = '{{offset:offset}} {{value.2}}'
                         )) {
  N_Groups = length(LIST_SampleNames)
  SCROLLY = "300px"
  PAGING = TRUE
  ORDER = list(N_Groups + 2, 'asc')
  PAGE_LENGTH = 5
  DOM = 'frtip'
  PROTEINS_RENDER = JS(
    "function(data, type, row, meta) {",
    "return '<span  data-toggle=\"popover\" data-trigger=\"hover click\" data-placement=\"right\" data-html=\"true\" data-delay={show: 500, hide: 100} title=\"Proteins\" data-content=\"' + data.replace(/;/g,\"<br>\") + '\">' + data.split(';').length + '</span>'",
    "}"
  )
  
  COLDEFS <- list()
  for (i in 1:N_Groups) {
    Targets = i + 1
    Width = '10%'
    Render = JS(
      paste(
        "function(data, type, full){return '<span class=",
        names(LIST_SampleNames)[i],
        "spark>' + data + '</span>'}",
        sep = ""
      )
    )
    COLDEFS[[i]] <-
      list(targets = Targets,
           width = Width,
           render = Render)
  }
  DRAW_CALLBACK_TEXT <- paste("function(){")
  for (i in 1:N_Groups) {
    DRAW_CALLBACK_TEXT <- paste(
      DRAW_CALLBACK_TEXT,
      "$('.",
      names(LIST_SampleNames)[i],
      "spark:not(:has(canvas))').sparkline('html', {
        type:'",
      LIST_SparlinesOptions$type,
      "',
        height:'",
      LIST_SparlinesOptions$height,
      "',
        width:'",
      LIST_SparlinesOptions$width,
      "',
        highlightColor:'",
      LIST_SparlinesOptions$highlightColor,
      "',
        chartRangeMin:",
      LIST_SparlinesOptions$chartRangeMin,
      ",
        chartRangeMax:",
      LIST_SparlinesOptions$chartRangeMax,
      ",
        tooltipFormat: '",
      LIST_SparlinesOptions$tooltipFormat,
      "',
        tooltipValueLookups: {
          'offset': {",
      paste(paste0(
        0:(length(LIST_SampleNames[[i]]) - 1), ": '", LIST_SampleNames[[i]], "'"
      ), collapse = ","),
      "}
        },
      });",
      sep = ""
    )
  }
  DRAW_CALLBACK_TEXT <- paste(DRAW_CALLBACK_TEXT, "}", sep = "\n")
  
  OPTIONS <- list(
    scrollY = SCROLLY,
    paging = PAGING,
    pageLength = PAGE_LENGTH,
    order = ORDER,
    dom = DOM,
    columnDefs = COLDEFS,
    drawCallback = JS(DRAW_CALLBACK_TEXT)
  )
  return(OPTIONS)
}


utils::globalVariables(
  c(
    "ProteinBankFastaFilePath",
    "PeptidePath",
    "SampleDescriptionPath",
    "IntensityName",
    "SampleColumnName",
    "SampleGroupColumnName"
  )
)

#' shiny server
#'
#' @param input shiny input
#' @param output shiny input
#' @param session shiny session
#' @importFrom utils read.table
#' @importFrom stats formula
#' @importFrom shinyjqui orderInput
#' @return shiny server
#' @export
#'
#' @examples
#' library("isoformnspectRe")
#' if(interactive()){
#' GlobalPath<-system.file("extdata",
#'  "Global.R",
#'   package = "isoformnspectRe")
#' source(GlobalPath,local=TRUE)
#' shiny::shinyApp(UI,server)
#' }
server <- function(input,
                   output,
                   session) {
  progress <- shiny::Progress$new()
  progress$set(message = "Read MaxQuant peptides file")
  MaxQuantPeptides <- utils::read.table(PeptidePath,
                                        header = TRUE,
                                        sep = "\t",
                                        quote = "\"")
  ColnamesToKeep <- c(
    "Sequence",
    "Proteins",
    "Leading.razor.protein",
    "Start.position",
    "End.position",
    "Unique..Proteins.",
    "PEP",
    colnames(MaxQuantPeptides)[grep(paste0("^", gsub(" ", ".", IntensityName)),
                                    colnames(MaxQuantPeptides))],
    "UpdateProteins",
    "NProteins"
  )
  progress$set(message = "Selecting columns of interest")
  MaxQuantPeptides <-
    MaxQuantPeptides[, colnames(MaxQuantPeptides) %in% ColnamesToKeep]
  progress$set(message = "Determining proteotypic peptides")
  ProteotypicProteins <-
    as.vector(MaxQuantPeptides[MaxQuantPeptides[, "Unique..Proteins."] == "yes", "Leading.razor.protein"])
  progress$set(message = "Reading sample description file")
  SampleDescription <- read.table(
    SampleDescriptionPath,
    header = TRUE,
    sep = "\t",
    quote = "\""
  )
  progress$set(message = "Handling intensities")
  MaxQuantPeptides <-
    MaxQuantPeptides[rowSums(MaxQuantPeptides[, colnames(MaxQuantPeptides)[grep(gsub(" ", ".", IntensityName), colnames(MaxQuantPeptides))]]) > 0,]
  Intensities <-
    paste(gsub(" ", ".", IntensityName),
          as.vector(SampleDescription[, SampleColumnName]),
          sep = "")
  #cat(Intensities)
  MaxIntensities <- apply(MaxQuantPeptides[, Intensities],
                          1, max)
  MaxQuantPeptides[, Intensities] <-
    round(MaxQuantPeptides[, Intensities] / MaxIntensities, digits = 2)
  IntensitiesPrefix <- gsub(" ", ".", IntensityName)
  progress$set(message = "Extracting sample group")
  SamplesByGroup <-
    lapply(unique(SampleDescription[, SampleGroupColumnName]),
           function(x, SampleDescription) {
             as.vector(SampleDescription[SampleDescription[, SampleGroupColumnName] == x, SampleColumnName])
           }, SampleDescription = SampleDescription)
  names(SamplesByGroup) <-
    unique(SampleDescription[, SampleGroupColumnName])
  NonSamplesByGroup <-
    lapply(unique(SampleDescription[, SampleGroupColumnName]),
           function(x, SampleDescription) {
             as.vector(SampleDescription[SampleDescription[, SampleGroupColumnName] != x, SampleColumnName])
           }, SampleDescription = SampleDescription)
  names(NonSamplesByGroup) <-
    unique(SampleDescription[, SampleGroupColumnName])
  SampleList <- c(SamplesByGroup,
                  NonSamplesByGroup)
  progress$set(message = "Creating informations columns")
  NSampleGroups <-
    length(unique(SampleDescription[, SampleGroupColumnName]))
  Progression <- 1
  withProgress(message = "Creating informations columns",
               min = 0,
               max = 1,
               {
                 lapply(unique(SampleDescription[, SampleGroupColumnName]), function(x, SampleList) {
                   incProgress(1 / NSampleGroups,
                               detail = paste("Handling group", unique(SampleDescription[, SampleGroupColumnName])[[Progression]]))
                   InterestIntensities <-
                     paste(IntensitiesPrefix,
                           SamplesByGroup[[x]],
                           sep = "")
                   UninterestIntensities <-
                     paste(IntensitiesPrefix,
                           NonSamplesByGroup[[x]],
                           sep = "")
                   #cat(UninterestIntensities)
                   if (length(InterestIntensities) > 1) {
                     MaxQuantPeptides[, paste(x, "Infos", sep = "_")] <<-
                       apply(MaxQuantPeptides[, InterestIntensities], 1
                             ,      function(y) {
                               paste(y, sep = ",", collapse = ",")
                             })
                     MaxQuantPeptides[, paste(x, "Mean", sep = "_")] <<-
                       round(apply(MaxQuantPeptides[, InterestIntensities],
                                   1, mean), digits = 2)
                   } else{
                     MaxQuantPeptides[, paste(x, "Infos", sep = "_")] <<-
                       MaxQuantPeptides[, InterestIntensities]
                     MaxQuantPeptides[, paste(x, "Mean", sep = "_")] <<-
                       round(MaxQuantPeptides[, InterestIntensities],
                             digits = 2)
                   }
                   MaxQuantPeptides[, paste("Non", x, "Mean", sep = "_")] <<-
                     round(apply(MaxQuantPeptides[, UninterestIntensities],
                                 1, mean), digits = 2)
                   Progression <<- Progression + 1
                 }, SampleList = SampleList)
               })
  progress$set(message = "Handling fold change")
  lapply(unique(SampleDescription[, SampleGroupColumnName]), function(x) {
    MaxQuantPeptides[, paste("Log 2 ", x, "/Non ", x, sep = "")] <<-
      HandleFoldChange(round(log(
        MaxQuantPeptides[, paste(x, "Mean", sep = "_")] / MaxQuantPeptides[, paste("Non", x, "Mean", sep = "_")], 2
      )))
  })
  progress$set(message = "Filtering intensities")
  
  MaxQuantPeptides <- MaxQuantPeptides[MaxIntensities > 0, ]
  
  MaxQuantPeptides[, "ProteinToPrint"] <-
    as.vector(MaxQuantPeptides[, "Leading.razor.protein"])
  Blast_Idx <- grep("Blast=", MaxQuantPeptides[, "ProteinToPrint"])
  Isoforms_Idx <-
    intersect(grep("-[0-9]*\\|", MaxQuantPeptides[, "ProteinToPrint"]),
              grep("UNIPROT=", MaxQuantPeptides[, "ProteinToPrint"]))
  Regular_Bool <-
    !1:length(MaxQuantPeptides[, "ProteinToPrint"]) %in% c(Blast_Idx, Isoforms_Idx)
  
  progress$set(message = "Creating HTML Link")
  
  #https://stackoverflow.com/questions/39039424/how-to-link-a-local-file-to-an-html-query-in-the-shiny-ui-r
  #points to a file in a www repertory located in the app.R file
  MaxQuantPeptides[Regular_Bool, "ProteinToPrint"] <- paste(
    "<a href=\"RegularSkeleton/",
    gsub("[[:punct:]]",
         "_", MaxQuantPeptides[Regular_Bool, "Leading.razor.protein"]),
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[[:punct:]]",
      replacement = " ",
      MaxQuantPeptides[Regular_Bool, "ProteinToPrint"]
    ),
    "</a>",
    sep = ""
  )
  
  MaxQuantPeptides[Blast_Idx, "ProteinToPrint"] <- paste(
    "<a href=\"BlastProtein/",
    gsub("[.,=\\|:]",
         "_", MaxQuantPeptides[Blast_Idx, "Leading.razor.protein"])
    ,
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[.,=\\|:]",
      replacement = " ",
      MaxQuantPeptides[Blast_Idx, "ProteinToPrint"]
    ),
    "</a>",
    sep = ""
  )
  
  MaxQuantPeptides[Isoforms_Idx, "ProteinToPrint"] <- paste(
    "<a href=\"IsoformProtein/",
    gsub("[.,=\\|:]",
         "_", MaxQuantPeptides[Isoforms_Idx, "Leading.razor.protein"])
    ,
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[.,=\\|:]",
      replacement = " ",
      MaxQuantPeptides[Isoforms_Idx, "ProteinToPrint"]
    ),
    "</a>",
    sep = ""
  )
  progress$set(message = "Reformating sequence")
  MaxQuantPeptides[, "Sequence"] <-
    as.vector(MaxQuantPeptides[, "Sequence"])
  MaxQuantPeptides[, "Sequence"] <-
    TrimSequenceOutput(MaxQuantPeptides[, "Sequence"])
  progress$close()
  SelectedProteins <- shiny::reactive({
    mRNA_Regexp <- "^str"
    UNIPROT_Regexp <- "\\|[A-Z0-9-]+\\|"
    Canonical_Regexp <- "sp\\|[A-Z0-9]+\\|"
    Isoform_Regexp <- "sp\\|[A-Z0-9]+-[0-9]+\\|"
    TrEMBL_Regexp <- "tr\\|[A-Z0-9]+\\|"
    PerfectMatch_Regexp <- "UNIPROT="
    PartialMatch_Regexp <- "Blast="
    Match_Regexp <- "="
    
    PresentProteins <-
      unique(as.vector(MaxQuantPeptides[, "Leading.razor.protein"]))
    LRPs <-
      unique(as.vector(MaxQuantPeptides[, "Leading.razor.protein"]))
    NoBlast <- LRPs[grep("Blast=", LRPs, invert = TRUE)]
    Shorts.NoBlast <-
      gsub("\\|.*$", "", gsub(pattern = "^[^\\|]*\\|", "", NoBlast))
    Blasts <- LRPs[grep("Blast=", LRPs)]
    Shorts.Blasts <-
      gsub("\\|.*$", "", gsub(pattern = "^[^\\|]*\\|", "", Blasts))
    BlastsWithCounterparts <-
      Blasts[Shorts.Blasts %in% Shorts.NoBlast]
    
    Isoforms <- LRPs[intersect(grep(PerfectMatch_Regexp, LRPs),
                               grep(Isoform_Regexp, LRPs))]
    Canonical <- LRPs[intersect(grep(PerfectMatch_Regexp, LRPs),
                                grep(Canonical_Regexp, LRPs))]
    
    Shorts.CanonicalFromIsoforms <-
      gsub("-.*$", "", gsub(pattern = "^[^\\|]*\\|", "", Isoforms))
    Shorts.Canonical <-
      gsub("\\|.*$", "", gsub(pattern = "^[^\\|]*\\|", "", Canonical))
    IsoformsWithCounterpart <-
      Isoforms[Shorts.CanonicalFromIsoforms %in% Shorts.Canonical]
    
    mRNA_Prot <-
      PresentProteins[grep(mRNA_Regexp, PresentProteins, invert = !input$mRNA)]
    
    UNIPROT_Prot <-
      PresentProteins[grep(UNIPROT_Regexp, PresentProteins, invert = !input$UNIPROT)]
    
    SelectedProt <- intersect(mRNA_Prot,
                              UNIPROT_Prot)
    #cat("Selected Prot", length(SelectedProt), "\n")
    if (input$UNIPROT) {
      if (input$Bank == "All") {
        Bank_Prot <- UNIPROT_Prot
      } else if (input$Bank == "Canonical") {
        Bank_Prot <-
          PresentProteins[grep(Canonical_Regexp, PresentProteins)]
      } else if (input$Bank == "Isoform") {
        Bank_Prot <-
          PresentProteins[grep(Isoform_Regexp, PresentProteins)]
        if (input$Counterpart) {
          Bank_Prot <-
            Bank_Prot[Bank_Prot %in% IsoformsWithCounterpart]
        }
      } else if (input$Bank == "TrEMBL") {
        Bank_Prot <-
          PresentProteins[grep(TrEMBL_Regexp, PresentProteins)]
      }
      SelectedProt <- intersect(SelectedProt,
                                Bank_Prot)
    }
    if (input$mRNA) {
      if (input$Match == "Both") {
        Matching_Prot <-
          PresentProteins[grep(Match_Regexp, PresentProteins)]
      } else if (input$Match == "Perfect match") {
        Matching_Prot <-
          PresentProteins[grep(PerfectMatch_Regexp, PresentProteins)]
      } else if (input$Match == "Blast") {
        Matching_Prot <-
          PresentProteins[grep(PartialMatch_Regexp, PresentProteins)]
        if (input$Counterpart) {
          Matching_Prot <-
            Matching_Prot[Matching_Prot %in% BlastsWithCounterparts]
        }
      }
      SelectedProt <- intersect(SelectedProt,
                                Matching_Prot)
    }
    GrpMean <- paste(input$Group, "Mean", sep = "_")
    NonGrpMean <- paste("Non", input$Group, "Mean", sep = "_")
    GrpLFC <-
      paste("Log 2 ", input$Group, "/Non ", input$Group, sep = "")
    InfosNames <-
      paste(input$SampleGroupsColumns_order, "Infos", sep = "_")
    KOL <- c("Sequence",
             "ProteinToPrint",
             InfosNames,
             GrpMean,
             NonGrpMean,
             GrpLFC)
    #cat(KOL)
    #cat(colnames(MaxQuantPeptides))
    #cat(KOL[!KOL %in% colnames(MaxQuantPeptides)])
    if (input$Proteotypic) {
      MaxQuantPeptides[MaxQuantPeptides[, "Leading.razor.protein"] %in% SelectedProt &
                         MaxQuantPeptides[, "Unique..Proteins."] == "yes",
                       KOL]
    } else{
      MaxQuantPeptides[MaxQuantPeptides[, "Leading.razor.protein"] %in% SelectedProt,
                       KOL]
    }
  })
  output$Button12 <- shiny::renderUI({
    ##req(input$pep)
    actionButton("page_12", "Proceed to sample description")
  })
  output$SampleGroupOrder <- shiny::renderUI({
    ColumnNames <-
      unique(as.vector(SampleDescription[, SampleGroupColumnName]))
    orderInput(inputId = "SampleGroupsColumns",
               label = "Reorder sample group order",
               items = ColumnNames)
  })
  
  
  output$Peptides <- DT::renderDT(
    SelectedProteins(),
    rownames = FALSE,
    escape = FALSE,
    selection = list(
      mode = "multiple",
      selected = event_data("plotly_selected", priority = "event")$pointNumber +
        1,
      target = 'row'
    ),
    options = LocalOptions(LIST_SampleNames = SamplesByGroup[input$SampleGroupsColumns_order])
  )
  
  
  proxy = DT::dataTableProxy('Peptides')
  observeEvent(input$clear, {
    proxy %>% selectRows(NULL)
  })
  output$Group <- shiny::renderUI({
    req(input$SampleGroupsColumns_order)
    selectInput(
      "Group",
      label = "Group of interest:",
      choices = input$SampleGroupsColumns_order,
      selected = input$SampleGroupsColumns_order[1]
    )
  })
  output$Description <- shiny::renderUI({
    laius <- ""
    if (input$Proteotypic) {
      laius <-
        paste(
          laius,
          "The displayed peptides are specific of the protein. The peptides are called proteotypic",
          sep = "<br>"
        )
    }
    if (input$mRNA) {
      laius <-
        paste(laius,
              "We found a transcript for the present proteins.",
              sep = "<br>")
    }
    if (input$UNIPROT) {
      laius <-
        paste(laius,
              "We found an UNIPROT entry for the present proteins.",
              sep = "<br>")
    }
    if (input$UNIPROT) {
      if (input$Bank == "All") {
        laius <-
          paste(laius,
                "The found proteins came from any part of UNIPROT.",
                sep = "<br>")
      } else if (input$Bank == "Canonical") {
        laius <-
          paste(laius,
                "The found proteins came from SwissProt canonical.",
                sep = "<br>")
      } else if (input$Bank == "Isoform") {
        laius <-
          paste(laius,
                "The found proteins came from SwissProt isoform.",
                sep = "<br>")
      } else if (input$Bank == "TrEMBL") {
        laius <-
          paste(laius, "The found proteins came from TrEMBL.", sep = "\n")
      }
      if (input$mRNA) {
        if (input$Match == "Both") {
          laius <-
            paste(
              laius,
              "The translated protein from the mRNA is identical or share a similarity with the one in UNIPROT.",
              sep = "<br>"
            )
        } else if (input$Match == "Perfect match") {
          laius <-
            paste(
              laius,
              "The translated protein from the mRNA is identical with the one in UNIPROT.",
              sep = "<br>"
            )
        } else if (input$Match == "Blast") {
          laius <-
            paste(
              laius,
              "The translated protein from the mRNA share a similarity with the one in UNIPROT.",
              sep = "<br>"
            )
        }
      }
      
    }
    if (input$Counterpart) {
      if (input$Match == "Blast") {
        laius <-
          paste(laius,
                "The protein matched by the blast is also present in the study.",
                sep = "<br>")
      }
      if (input$Bank == "Isoform") {
        laius <-
          paste(laius,
                "The canonical form of the isoform is also present in the study.",
                sep = "<br>")
      }
    }
    HTML(laius)
  })
  output$scatterPlotly <- plotly::renderPlotly({
    SP <- SelectedProteins()
    p <- plot_ly(
      data = SP,
      x = formula(paste("~", input$Group, "_Mean", sep = "")),
      y = formula(paste("~Non_", input$Group, "_Mean", sep = "")),
      text = ~ Sequence,
      mode = "markers",
      type = "scatter",
      marker = list(opacity = 0.2, color = "black")
    )
    p <- layout(p, showlegend = FALSE)
    p <- toWebGL(p)
    s <- input$Peptides_rows_selected
    if (!length(s))
      return(p)
    add_trace(
      p,
      data = SP[s, , drop = FALSE],
      x = formula(paste("~", input$Group, "_Mean", sep = "")),
      y = formula(paste("~Non_", input$Group, "_Mean", sep = "")),
      type = "scatter",
      mode = "markers",
      marker = list(
        opacity = 1,
        color = "red",
        size = 10
      )
    )
  })
  session$onSessionEnded(function() {
    stopApp()
  })
}
