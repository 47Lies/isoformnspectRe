#' Compute the elements that represent the blast configuration for a interactive plotly graph
#'
#' @param QuerySeqId Name of the query (submitted sequence)
#' @param QueryLength Length of the query (full length not the mapped part)
#' @param QueryStart Beggining postition of the matching part in the query
#' @param QueryEnd Ending position of the matching part in the query
#' @param SubjectSeqId Name of the subject (reference sequence)
#' @param SubjectStart Beggining postition of the matching part in the subject
#' @param SubjectEnd Ending position of the matching part in the subject
#' @param SubjectLength Length of the query (full length not the mapped part)
#' @param Length Length of the matching part between subject and query
#' @param Height height of the widget
#' @param Step fraction of the gap between the shared parts and the specific parts
#' @param LineDots inneer graduation set in AA
#' @rawNamespace import(plotly, except = c(slice,config))
#' @import viridisLite
#' @return List of Dataframes for a plotly graph
#' @export
#'
#' @examples
#' library(plotly)
#' library(viridisLite)
#' BlastPalette<-viridis(3)
#' R_FilePath <- system.file("extdata",
#'  "blast_rescue.example.txt",
#'   package = "isoformnspectRe")
#' Rescue<-read.table(
#'  R_FilePath,
#'  header=TRUE)
#'  Infos <- Rescue[1,]
#'  DFs <- BlastCoordPlotly(Infos$QuerySeqId, Infos$QueryLength, Infos$QueryStart,
#'  Infos$QueryEnd, Infos$SubjectSeqId, Infos$SubjectStart, Infos$SubjectEnd,
#'  Infos$SubjectLength, Infos$Length)
#'  names(BlastPalette) <- c(as.vector(Infos$QuerySeqId) ,
#'  paste( as.vector(Infos$QuerySeqId), "<br>", as.vector(Infos$SubjectSeqId), sep = "" ),
#'  as.vector(Infos$SubjectSeqId))
#'  P<-plotly::plot_ly()
#'  P <- plotly::add_trace(P,
#'  x = DFs$MatchBlock$x,
#'  y = DFs$MatchBlock$y,
#'  color = I(as.vector(BlastPalette[2])), type = "scatter", mode = "lines",
#'  fill = "toself", name = paste( as.vector(Infos$QuerySeqId), "<br>",
#'  as.vector(Infos$SubjectSeqId), sep = ""), hoverinfo="name",hoverlabel=list(namelength=-1))
#'  P
#'  P <- plotly::add_trace(P, x = DFs$SubjectBlock$x, y = DFs$SubjectBlock$y,
#'  color = I(as.vector(BlastPalette[3])), type = "scatter", mode = "lines", fill = "toself",
#'  name = Infos$SubjectSeqId, hoverinfo="name", hoverlabel=list(namelength=-1))
#'  P
#'  P <- plotly::add_trace(P, x = DFs$QueryBlock$x, y = DFs$QueryBlock$y,
#'  color = I(as.vector(BlastPalette[1])), type = "scatter", mode = "lines", fill = "toself",
#'   name = Infos$SubjectSeqId,
#'  hoverinfo="name", hoverlabel=list(namelength=-1))
#'  P
#'  P<-plotly::add_segments(P, x=~x, y=~y, xend=~xend, yend=~yend,
#'  hoverinfo="text", color=I(BlastPalette[3]), text=~text,
#'  data=DFs$SubjectConnectors)
#'  P
#'  P<-plotly::add_segments(P, x=~x, y=~y, xend=~xend, yend=~yend,
#'  hoverinfo="text", color=I(BlastPalette[1]), text=~text,
#'  data=DFs$QueryConnectors)
#'  P
#'  P<-plotly::hide_legend(P)
#'  P
BlastCoordPlotly<-function(QuerySeqId,
                           QueryLength,
                           QueryStart,
                           QueryEnd,
                           SubjectSeqId,
                           SubjectStart,
                           SubjectEnd,
                           SubjectLength,
                           Length,
                           Height=0.5,
                           Step=10,
                           LineDots=100){
  Gap<-round(Length/Step)
  #Alignement part shared between the query and the subject
  MatchBlock <- data.frame(
    x = c(1,Length,Length,1,1,NA),
    y = c(0,0,Height,Height,0,NA)
  )
  QueryBlock<-data.frame(
    x=numeric(),
    y=numeric()
  )
  SubjectBlock<-data.frame(
    x=numeric(),
    y=numeric()
  )
  #If there no specific part, should stay empty otherwise will contain
  # the connectors between every chunks
  QueryConnectors<-data.frame(
    x=numeric(),
    y=numeric(),
    xend=numeric(),
    yend=numeric(),
    text=character()
  )
  SubjectConnectors<-data.frame(
    x=numeric(),
    y=numeric(),
    xend=numeric(),
    yend=numeric(),
    text=character()
  )

  #Query's alignement begin after the begining of the query
  # i.e. some part of the query is specific before the alignement
  # i.e. N terminal part of the tested protein is not retrievev in alignement
  if(QueryStart!=1) {
    Xs<-c(-QueryStart - Gap + 1,
          -Gap,
          -Gap,
          -QueryStart - Gap + 1,
          -QueryStart - Gap + 1,
          NA)
    Ys <- c(1,
            1,
            1 + Height,
            1 + Height,
            1,
            NA)
    QueryBlock <- rbind(
      QueryBlock,
      data.frame(
        x = Xs,
        y = Ys
      )
    )
    QueryConnectors<-rbind(
      QueryConnectors,
      data.frame(
        x = -Gap,
        y = 1+Height/2,
        xend = 1,
        yend = Height/2,
        text = QuerySeqId
      )
    )
  }
  #Query's alignement ends before the end of the query
  # i.e. some part of the query is specific after the alignement
  # i.e. C Terminal part of the tested protein not retrieved in alignement
  if(QueryEnd != QueryLength) {
    Xs<-c(Length + Gap,
          QueryLength - QueryEnd + Length + Gap,
          QueryLength - QueryEnd + Length + Gap,
          Length + Gap,
          Length + Gap,
          NA)
    Ys <- c(1,
            1,
            1 + Height,
            1 + Height,
            1,
            NA)
    QueryBlock <- rbind(
      QueryBlock,
      data.frame(
        x = Xs,
        y = Ys
      )
    )
    QueryConnectors <- rbind(
      QueryConnectors,
      data.frame(
        x = Length + Gap,
        y = 1+Height/2,
        xend = Length,
        yend = Height/2,
        text = QuerySeqId
      )
    )
  }
  #Subject's alignement ends before the end of the Subject
  # i.e. some part of the subject is specific after the alignement
  # i.e. C terminal part of reference protein not retrieved in alignement
  if (SubjectEnd  != SubjectLength) {
    Xs<-c(Length + Gap,
          SubjectLength - SubjectEnd + Length + Gap,
          SubjectLength - SubjectEnd + Length + Gap,
          Length + Gap,
          Length + Gap,
          NA)
    Ys <- c(-1,
            -1,
            -1 + Height,
            -1 + Height,
            -1,
            NA)
    SubjectBlock <- rbind(
      SubjectBlock,
      data.frame(
        x = Xs,
        y = Ys
      )
    )
    SubjectConnectors <- rbind(
      SubjectConnectors,
      data.frame(
        x = Length + Gap,
        y = -1 + Height/2,
        xend = Length,
        yend = Height/2,
        text = SubjectSeqId
      ))
  }
  #Subject's alignement after the begining of the Subject
  # i.e. some part of the subject is specific before the alignement
  # i.e. N terminal part of reference protein not retrieved in alignement
  if (SubjectStart  != 1) {
    Xs<-c(-SubjectStart-Gap+1,
          -Gap,
          -Gap,
          -SubjectStart-Gap+1,
          -SubjectStart-Gap+1,
          NA)
    Ys <- c(-1,
            -1,
            -1 + Height,
            -1 + Height,
            -1,
            NA)
    SubjectBlock <- rbind(
      SubjectBlock,
      data.frame(
        x = Xs,
        y = Ys
      )
    )
    SubjectConnectors <- rbind(
      SubjectConnectors,
      data.frame(
        x = -Gap,
        y = -1 + Height/2,
        xend = 1,
        yend = Height/2,
        text = SubjectSeqId
      )
    )
  }
  return(
    list(
      MatchBlock = MatchBlock,
      QueryBlock = QueryBlock,
      SubjectBlock = SubjectBlock,
      SubjectConnectors = SubjectConnectors,
      QueryConnectors = QueryConnectors
    )
  )
}
