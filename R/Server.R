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
    "SampleGroupColumnName",
    ":=",
    "!!",
    "Leading razor protein",
    "Unique (Proteins)"
    
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
#' @rawNamespace import(data.table, except = c(shift))
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
  MaxIntensity <- NNonNullIntensity <- Protein <- Sequence <- NULL
  progress <- shiny::Progress$new()
  progress$set(message = "Read MaxQuant peptides file")
  # MaxQuantPeptides <- utils::read.table(PeptidePath,
  #                                       header = TRUE,
  #                                       sep = "\t",
  #                                       quote = "\"")
  options(datatable.integer64 = "numeric")
  MaxQuantPeptides <- data.table::fread(PeptidePath,
                                        header = TRUE,
                                        sep = "\t",
                                        quote = "\"")
  if(length(grep(" $",IntensityName))==0){
    IntensityName<-paste0(IntensityName," ")
  }
  Intensities <- colnames(MaxQuantPeptides)[grep(IntensityName,
                                                 colnames(MaxQuantPeptides))]
  ColnamesToKeep <- c(
    "Sequence",
    "Proteins",
    "Leading razor protein",
    "Start position",
    "End position",
    "Unique (Proteins)",
    "PEP",
    Intensities,
    "UpdateProteins",
    "NProteins"
  )
  ColnamesToKeep <-
    intersect(ColnamesToKeep, colnames(MaxQuantPeptides))
  progress$set(message = "Selecting columns of interest")
  MaxQuantPeptides <-
    MaxQuantPeptides[, ColnamesToKeep, with = FALSE]
  progress$set(message = "Determining proteotypic peptides")
  ProteotypicProteins <-
    unlist(MaxQuantPeptides[`Unique (Proteins)` == "yes", "Leading razor protein"], use.names = FALSE)
  
  progress$set(message = "Reading sample description file")
  SampleDescription <- data.table::fread(SampleDescriptionPath,
                                         sep = "\t",
                                         quote = "\"")
  progress$set(message = "Handling intensities")
  #MaxQuantPeptides <-
  #  MaxQuantPeptides[rowSums(MaxQuantPeptides[, colnames(MaxQuantPeptides)[grep(gsub(" ", ".", IntensityName), colnames(MaxQuantPeptides))]]) > 0,]
  MaxQuantPeptides[, NNonNullIntensity := sum(as.double(.SD) > 0), by =  seq_len(nrow(MaxQuantPeptides)), .SDcols =
                     Intensities]
  MaxQuantPeptides <- MaxQuantPeptides[NNonNullIntensity > 0]
  #MaxIntensities <- apply(MaxQuantPeptides[, Intensities],
  #                        1, max)
  MaxQuantPeptides[, MaxIntensity := max(as.double(.SD)), by =  seq_len(nrow(MaxQuantPeptides)), .SDcols =
                     Intensities]
  #MaxQuantPeptides[, Intensities] <-
  #  round(MaxQuantPeptides[, Intensities] / MaxIntensities, digits = 2)
  #apply(Intensities,function(x){
  #      MaxQuantPeptides[,(x):=round(.SD/MaxIntensity,2),.SDcols=x]})
  MaxQuantPeptides[, (Intensities) := lapply(.SD, function(x)
    round(x / MaxIntensity, 2)), .SDcols = Intensities]
  
  IntensitiesPrefix <- IntensityName
  progress$set(message = "Extracting sample group")
  
  GRP <-
    unlist(unique(SampleDescription[, SampleGroupColumnName, with = FALSE]), use.names =
             FALSE)
  
  SamplesByGroup <-
    lapply(GRP,
           function(x, SampleDescription) {
             unlist(SampleDescription[get(SampleGroupColumnName) == x,][, .SD, .SDcols = SampleColumnName], use.names = FALSE)
           }, SampleDescription = SampleDescription)
  names(SamplesByGroup) <- GRP
  NonSamplesByGroup <-
    lapply(GRP,
           function(x, SampleDescription) {
             unlist(SampleDescription[get(SampleGroupColumnName) != x,][, .SD, .SDcols = SampleColumnName], use.names = FALSE)
           }, SampleDescription = SampleDescription)
  names(NonSamplesByGroup) <- GRP
  progress$set(message = "Creating informations columns")
  NSampleGroups <- length(GRP)
  Progression <- 1
  withProgress(message = "Creating informations columns",
               min = 0,
               max = 1,
               {
                 lapply(GRP, function(x) {
                   incProgress(1 / NSampleGroups,
                               detail = paste("Handling group", GRP[Progression]))
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
                     MaxQuantPeptides[, paste0(x, "_Infos") := paste(.SD, collapse = ","), by =
                                        seq_len(nrow(MaxQuantPeptides)), .SDcols = InterestIntensities]
                     MaxQuantPeptides[, paste0(x, "_Mean") := round(mean(as.numeric(.SD)),2), by =
                                        seq_len(nrow(MaxQuantPeptides)), .SDcols = InterestIntensities]
                     if(length(UninterestIntensities) > 1){
                       MaxQuantPeptides[, paste0("Non_",x, "_Mean") := round(mean(as.numeric(.SD)),2), by =
                                        seq_len(nrow(MaxQuantPeptides)), .SDcols = UninterestIntensities]
                     }else{
                       MaxQuantPeptides[, paste0("Non_",x, "_Mean") := .SD, by =
                                          seq_len(nrow(MaxQuantPeptides)), .SDcols = UninterestIntensities]
                     }
                   } else{
                     MaxQuantPeptides[, paste0(x, "_Infos") := .SD, by =
                                        seq_len(nrow(MaxQuantPeptides)), .SDcols = InterestIntensities]
                     MaxQuantPeptides[, paste0(x, "_Mean") := .SD, by =
                                        seq_len(nrow(MaxQuantPeptides)), .SDcols = InterestIntensities]
                     if(length(UninterestIntensities) > 1){
                       MaxQuantPeptides[, paste0("Non_",x, "_Mean") := round(mean(as.numeric(.SD)),2), by =
                                          seq_len(nrow(MaxQuantPeptides)), .SDcols = UninterestIntensities]
                     }else{
                       MaxQuantPeptides[, paste0("Non_",x, "_Mean") := .SD, by =
                                          seq_len(nrow(MaxQuantPeptides)), .SDcols = UninterestIntensities]
                     }
                   }
                   MaxQuantPeptides[,paste("Log 2 ", x, "/Non ", x):=HandleFoldChange(
                     round(
                       log(get(paste0(x, "_Mean"))/get(paste0("Non_",x, "_Mean")),2)
                       ,2)
                   )]
                   Progression <<- Progression + 1
                 })
               })
  #progress$set(message = "Handling fold change")
  # lapply(unique(SampleDescription[, SampleGroupColumnName]), function(x) {
  #   MaxQuantPeptides[, paste("Log 2 ", x, "/Non ", x, sep = "")] <<-
  #     HandleFoldChange(round(log(
  #       MaxQuantPeptides[, paste(x, "Mean", sep = "_")] / MaxQuantPeptides[, paste("Non", x, "Mean", sep = "_")], 2
  #     )))
  # })

  progress$set(message = "Filtering intensities")
  
  MaxQuantPeptides <- MaxQuantPeptides[MaxIntensity > 0,]
  
  #MaxQuantPeptides[, "Protein"] <-
  #  as.vector(MaxQuantPeptides[, "Leading.razor.protein"])
  MaxQuantPeptides[,"Protein":=`Leading razor protein`]
  BoolBlast <- grepl("Blast=", MaxQuantPeptides$Protein)
  BoolIsoforms <-
    grepl("-[0-9]*\\|", MaxQuantPeptides$Protein) &
              grepl("UNIPROT=", MaxQuantPeptides$Protein)
  BoolRegular <- !BoolIsoforms & !BoolBlast

  progress$set(message = "Creating HTML Link")
  
  #https://stackoverflow.com/questions/39039424/how-to-link-a-local-file-to-an-html-query-in-the-shiny-ui-r
  #points to a file in a www repertory located in the app.R file
  MaxQuantPeptides[BoolRegular, Protein:=paste(
    "<a href=\"RegularSkeleton/",
    gsub("[[:punct:]]",
         "_",`Leading razor protein`),
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[[:punct:]]",
      replacement = " ",
      `Leading razor protein`
    ),
    "</a>",
    sep = ""
  )]
  MaxQuantPeptides[BoolBlast, Protein:=paste(
    "<a href=\"BlastProtein/",
    gsub("[[:punct:]]",
         "_",`Leading razor protein`),
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[[:punct:]]",
      replacement = " ",
      `Leading razor protein`
    ),
    "</a>",
    sep = ""
  )]
  MaxQuantPeptides[BoolIsoforms, Protein:=paste(
    "<a href=\"IsoformProtein/",
    gsub("[[:punct:]]",
         "_",`Leading razor protein`),
    ".html\" target=\"_blank\">",
    gsub(
      pattern = "[[:punct:]]",
      replacement = " ",
      `Leading razor protein`
    ),
    "</a>",
    sep = ""
  )]
  progress$set(message = "Reformating sequence")
  MaxQuantPeptides[, Sequence:=TrimSequenceOutput(Sequence)]
  progress$close()
  SelectedProteins <- shiny::reactive({
    req(input$Group)
    mRNA_Regexp <- "^str"
    UNIPROT_Regexp <- "\\|[A-Z0-9-]+\\|"
    Canonical_Regexp <- "sp\\|[A-Z0-9]+\\|"
    Isoform_Regexp <- "sp\\|[A-Z0-9]+-[0-9]+\\|"
    TrEMBL_Regexp <- "tr\\|[A-Z0-9]+\\|"
    PerfectMatch_Regexp <- "UNIPROT="
    PartialMatch_Regexp <- "Blast="
    Match_Regexp <- "="
    
    PresentProteins <-
      unique(unlist(MaxQuantPeptides[, "Leading razor protein"]))
    LRPs <-
      unique(unlist(MaxQuantPeptides[, "Leading razor protein"]))
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
    GrpMean <- base::paste(as.vector(input$Group), "Mean", sep = "_")
    NonGrpMean <- base::paste("Non", as.vector(input$Group), "Mean", sep = "_")
    GrpLFC <-
      base::paste("Log 2 ", as.vector(input$Group), "/Non ", as.vector(input$Group), sep = " ")
    InfosNames <-
      base::paste(as.vector(input$SampleGroupsColumns_order[["text"]]), "Infos", sep = "_")
    KOL <- c("Sequence",
             "Protein",
             InfosNames,
             GrpMean,
             NonGrpMean,
             GrpLFC)

    #cat(KOL)
    #cat(colnames(MaxQuantPeptides))
    #cat(KOL[!KOL %in% colnames(MaxQuantPeptides)])
    if (input$Proteotypic) {
      return(MaxQuantPeptides[`Leading razor protein` %in% SelectedProt &
                         `Unique (Proteins)` == "yes", KOL, with = FALSE])
    } else{
      return(MaxQuantPeptides[`Leading razor protein` %in% SelectedProt,
                       KOL,with=FALSE])
    }
  })
  output$Button12 <- shiny::renderUI({
    ##req(input$pep)
    actionButton("page_12", "Proceed to sample description")
  })
  output$SampleGroupOrder <- shiny::renderUI({
    ColumnNames <-
      unique(unlist(SampleDescription[,.SD,.SDcols=SampleGroupColumnName],use.names=FALSE))
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
    options = LocalOptions(LIST_SampleNames = SamplesByGroup[input$SampleGroupsColumns_order$text])
  )
  
  
  proxy = DT::dataTableProxy('Peptides')
  observeEvent(input$clear, {
    proxy %>% selectRows(NULL)
  })
  output$Group <- shiny::renderUI({
    req(input$SampleGroupsColumns_order)
    #cat(names(input$SampleGroupsColumns_order))
    #cat(unlist(input$SampleGroupsColumns_order,names=FALSE))
    selectInput(
      inputId="Group",
      label = "Group of interest:",
      choices = unlist(input$SampleGroupsColumns_order$text,use.names = FALSE),
      multiple=FALSE
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
