#' LIST_ComputeDatatableOptions Insert data table options for
#' automatic inclusion of sparkline javascript code. Need some very secific DATA
#'
#'
#' @param LIST_SampleNames List of sample names that will be included for the sparkline
#' @param type if it has to be included in a skeleton part or if it has to be included in the app
#' @param LIST_SparlinesOptions options to be included in the sparkline
#'
#' @return List of options that will be included in a DT::data.table
#' @export
#' @importFrom htmlwidgets getDependency
#' @importFrom htmlwidgets JS
#'
#'
#' @examples
#' library(data.table)
#' library(DT)
#' library(sparklines)
#' library(htmlwidgets)
#' P_FilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' Peptides<-fread(
#'  P_FilePath,
#'  integer64="double")
#' SampleNames<-list("siNT"=c("siNT_48-1",
#' "siNT_48-2",
#' "siNT_48-3",
#' "siNT_72-1",
#' "siNT_72-2",
#' "siNT_72-3"),
#' "siSF3B1"=c("siSF3B1-1",
#' "siSF3B1-2",
#' "siSF3B1-3",
#' "siSRSF1-1",
#' "siSRSF1-2",
#' "siSRSF1-3"))
#' Samples <- unlist(SampleNames)
#' MaxIntensities <- apply(Peptides[,paste("Intensity",Samples),with=FALSE],1,max)
#' Peptides[,
#'  (paste("Intensity",Samples)) := lapply(.SD, function(x) x/MaxIntensities),
#'   .SDcols = paste("Intensity",Samples)]
#' Peptides[,"siNT Infos" := paste(.SD, collapse = ","), .SDcols =
#' paste("Intensity",SampleNames[["siNT"]]), by = Sequence]
#' Peptides[,"siSF3B1 Infos" := paste(.SD, collapse = ","), .SDcols =
#' paste("Intensity",SampleNames[["siSF3B1"]]), by = Sequence]
#' PeptideColumnOfInterest<-c("Proteins","Sequence",
#' "siNT Infos","siSF3B1 Infos",
#' "Start position","End position","PEP")
#' LocalDatatable <- DT::datatable(
#' Peptides[,..PeptideColumnOfInterest],
#' width = "100%",
#' options = LIST_ComputeDatatableOptions(LIST_SampleNames=SampleNames),
#' rownames = FALSE, escape = FALSE)
#' LocalDatatable$dependencies <- append(LocalDatatable$dependencies,
#'  htmlwidgets:::getDependency("sparkline"))
#' LocalDatatable
#'
LIST_ComputeDatatableOptions<-function(
  LIST_SampleNames,
  type="Skeleton",
  LIST_SparlinesOptions=list(
    type="bar",
    height=40,
    width=60,
    highlightColor="red",
    chartRangeMin=0,
    chartRangeMax=1,
    tooltipFormat='{{offset:offset}} {{value.2}}')){
  N_Groups=length(LIST_SampleNames)
  if(type=="Skeleton"){
    SCROLLY="750px"
    PAGING=FALSE
    ORDER=list(N_Groups+2,'asc')
    PAGE_LENGTH = 20
    DOM='frtip'
    PROTEINS_RENDER=htmlwidgets::JS(
      "function(data, type, row, meta) {",
      "return '<span  data-toggle=\"popover\" data-trigger=\"hover click\" data-placement=\"right\" data-html=\"true\" data-delay={show: 500, hide: 100} title=\"Proteins\" data-content=\"' + data.replace(/;/g,\"<br>\") + '\">' + data.split(';').length + '</span>'",
      "}")
  }else{
    SCROLLY=NULL
    PAGING=TRUE
    PAGE_LENGTH = 5
    DOM='frtip'
    ORDER=list(N_Groups=3,'desc')
    PROTEINS_RENDER=JS(
      "function(data, type, row, meta) {",
      "return '<span onclick=\"GetCopy(this)\" title=\"' + data + '\">' + data.substr(0, 10)+'...' + '</span>'",
      "}")

  }
  #LIST_SampleNames=SamplesByGroup
  COLDEFS<-list(
    list(targets=c(0),render=PROTEINS_RENDER)
  )
  for (i in 1:N_Groups) {
    Targets = i + 1
    Width='10%'
    Render = JS(
      paste(
        "function(data, type, full){return '<span class=",names(LIST_SampleNames)[i],"spark>' + data + '</span>'}",
        sep = ""
      )
    )
    COLDEFS[[i + 1]] <- list(targets = Targets, width=Width, render = Render)
  }
  DRAW_CALLBACK_TEXT<-paste("function(){")
  for(i in 1:N_Groups){
    DRAW_CALLBACK_TEXT<-paste(
      DRAW_CALLBACK_TEXT,
      "$('.",names(LIST_SampleNames)[i],"spark:not(:has(canvas))').sparkline('html', {
        type:'",LIST_SparlinesOptions$type,"',
        height:'",LIST_SparlinesOptions$height,"',
        width:'",LIST_SparlinesOptions$width,"',
        highlightColor:'",LIST_SparlinesOptions$highlightColor,"',
        chartRangeMin:",LIST_SparlinesOptions$chartRangeMin,",
        chartRangeMax:",LIST_SparlinesOptions$chartRangeMax,",
        tooltipFormat: '",LIST_SparlinesOptions$tooltipFormat,"',
        tooltipValueLookups: {
          'offset': {",
      paste(paste0(0:(length(LIST_SampleNames[[i]])-1), ": '", LIST_SampleNames[[i]], "'"), collapse = ","),
      "}
        },
      });",sep="")
  }
  DRAW_CALLBACK_TEXT<-paste(DRAW_CALLBACK_TEXT,"$('[data-toggle=\"popover\"]').popover({container: 'body'});}",sep="\n")

  OPTIONS<-list(
    scrollY = SCROLLY,
    paging = PAGING,
    pageLength = PAGE_LENGTH,
    order = ORDER,
    dom=DOM,
    columnDefs=COLDEFS,
    #columnDefs=list(list(targets=c(2,3,4,5),render=JS("function(data, type, full){return '<span class=spark>' + data + '</span>'}"))),
    # #drawCallback=JS("function(){
    #   $('.WNTspark:not(:has(canvas))').sparkline('html', {
    #     type:'bar',
    #     height:'40',
    #     width:'60',
    #     highlightColor:'red',
    #     chartRangeMin:0,
    #     chartRangeMax:1});
    #   $('[data-toggle=\"popover\"]').popover({container: 'body'});}")
    #columnDefs=COLDEFS,
    drawCallback=JS(DRAW_CALLBACK_TEXT)
  )
  return(OPTIONS)
}
