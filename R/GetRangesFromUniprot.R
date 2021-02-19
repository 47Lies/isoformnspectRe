#' GetRangesFromUniprot From a Uniprot isoform entry get the variations from the canonical sequence. The start and end coordonate are isoform sequence based.
#'
#' @param UNIPROT_ISOFORM_ENTRY Isoform entry name typically a capital letter
#'
#' @return IRanges object with each line representing the variant splicing of the Isoform. The class and the Title metadata columns are ready for a HTML span inclusion.
#' @import httr
#' @importFrom IRanges IRanges
#' @importFrom IRanges shift
#' @importFrom GenomicRanges mcols
#' @import xml2
#' @export
#' @examples
#' VariationRanges<-GetRangesFromUniprot("O94875-10")
#' VariationRanges
GetRangesFromUniprot <-
  function(UNIPROT_ISOFORM_ENTRY = "O94875-10") {
    #https://stackoverflow.com/questions/32019566/r-xml-parse-for-a-web-address
    URL <-
      paste("https://www.uniprot.org/uniprot/",
            UNIPROT_ISOFORM_ENTRY,
            ".xml",
            sep = "")
    #internet sometimes the URL won't be accessible
    #so retry a few times before crash
    #
    XMLFromUniprot <-
      httr::RETRY("GET", URL, times = 5, pause_min = 5)
    ## Apparently content from httr will return an xml2 handler, thank you Hadley
    UniprotXMLDoc <- httr::content(XMLFromUniprot)
    #https://rdrr.io/cran/xml2/man/xml_ns_strip.html
    xml2::xml_ns_strip(UniprotXMLDoc)
    SpliceVariantNodes <-
      xml2::xml_find_all(UniprotXMLDoc, "//feature[@type='splice variant']")
    SpliceVariantList <- xml2::as_list(SpliceVariantNodes)
    names(x = SpliceVariantList) <-
      xml2::xml_attr(SpliceVariantNodes, "id")
    IsoformInfoXPath <-
      paste("//isoform[id='",
            UNIPROT_ISOFORM_ENTRY,
            "']/sequence",
            sep = "")
    LocalIsoformNode <-
      xml2::xml_find_first(UniprotXMLDoc, IsoformInfoXPath)
    CanonicalSequence <-
      xml2::xml_text(xml2::xml_find_first(UniprotXMLDoc, "//sequence[@length]"))
    LocalIsoformVariations <-
      xml2::xml_attr(LocalIsoformNode, "ref")
    LocalIsoformVariations <-
      unlist(strsplit(split = " ", LocalIsoformVariations))
    #https://rdrr.io/cran/xml2/man/as_list.html
    #attributes, didn't know it even exists until I meet it
    #https://www.dummies.com/programming/r/how-to-play-with-attributes-in-r/
    #Inelegant but can't do better since GenomicRanges::IRanges do not support unlist for now
    #https://support.bioconductor.org/p/109137/

    #Q93008-1 case no variation since it's the canonical form but still considered as isoform
    if (sum(is.na(LocalIsoformVariations)) != 0) {
      return(GetEmptySpanRanges())
    }

    VariationsRanges <- IRanges::IRanges()
    TrashToAvoidPrint <- lapply(LocalIsoformVariations,
                                function(x) {
                                  #x<-"VSP_047405"
                                  B <-
                                    as.numeric(attr(SpliceVariantList[[x]]$location$begin, "position"))
                                  E <-
                                    as.numeric(attr(SpliceVariantList[[x]]$location$end, "position"))
                                  #Insertion Case where a single base is substitute by a whole segement, so no begin node no start node but a position in the xml
                                  # variant "VSP_035381" in Q99592-2 case
                                  if (length(B) == 0 &&
                                      length(E) == 0) {
                                    B <-
                                      as.numeric(attr(SpliceVariantList[[x]]$location$position, "position"))
                                    E <- B
                                  }
                                  Class <- "VariantSplicing"
                                  ShiftInducedByVariation <-
                                    E - B + 1
                                  Original <-
                                    SpliceVariantList[[x]]$original[[1]]
                                  if (is.null(Original)) {
                                    Original <- substr(x = CanonicalSequence,
                                                       start = B,
                                                       stop = E)
                                  }
                                  Variation <-
                                    SpliceVariantList[[x]]$variation[[1]]
                                  ShiftCompensatedByVariation <- 0
                                  if (is.null(Variation)) {
                                    Variation <- ""
                                    E = B
                                  }
                                  ShiftCompensatedByVariation <-
                                    nchar(Variation)
                                  E = B + ShiftCompensatedByVariation
                                  Range <- IRanges::IRanges(
                                    start = B,
                                    end = E,
                                    names = x,
                                    Class = Class,
                                    Original = Original,
                                    Variation = Variation,
                                    Loss = ShiftInducedByVariation,
                                    Compensation = ShiftCompensatedByVariation
                                  )
                                  VariationsRanges <<-
                                    c(VariationsRanges, Range)
                                  return(NULL)
                                })
    VariationShift <-
      cumsum(
        -GenomicRanges::mcols(VariationsRanges)$Loss + GenomicRanges::mcols(VariationsRanges)$Compensation
      )
    VariationShift <- c(0, VariationShift)
    VariationShift <- VariationShift[-length(VariationShift)]
    VariationsRanges <-
      IRanges::shift(VariationsRanges, VariationShift)
    GenomicRanges::mcols(VariationsRanges)$Title <-
      paste(
        "Variation:",
        GenomicRanges::mcols(VariationsRanges)$Original,
        "->",
        GenomicRanges::mcols(VariationsRanges)$Variation,
        "\n",
        sep = ""
      )
    GenomicRanges::mcols(VariationsRanges)[, c("Original", "Variation", "Loss", "Compensation")] <-
      NULL
    return(VariationsRanges)
  }


#' GetEmptySpanRanges Get an empty IRange
#'
#' @return span an empty Genomic Range with class and title attribute for HTML span inclusion
#' @import IRanges
#' @import GenomicRanges
#' @export
#' @examples
#' EmptyRange<-GetEmptySpanRanges()
#' EmptyRange
GetEmptySpanRanges <- function() {
  Span <- IRanges::IRanges()
  GenomicRanges::mcols(Span)$Class = character()
  GenomicRanges::mcols(Span)$Title = character()
  return(Span)
}


#' CollapseSpanRanges Handle a set of range to have only one range per position. This will allow us to use the
#'
#' @param SpanRanges IRanges that will be used for HMTL formatting
#'
#' @return IRanges with class and title information for HTML span tag
#' @import GenomicRanges
#' @importFrom IRanges disjoin
#' @export
#' @examples
#' library(GenomicRanges)
#'
#' VariationRanges<-GetRangesFromUniprot("O94875-10")
#' VariationRanges
#' DuplicatedVariationRanges<-VariationRanges
#' mcols(DuplicatedVariationRanges)$Title<-paste("Duplicated ",mcols(DuplicatedVariationRanges)$Title)
#' VariationRanges<-c(VariationRanges,DuplicatedVariationRanges)
#' Collapsed<-CollapseSpanRanges(VariationRanges)
#' VariationRanges
#' Collapsed
CollapseSpanRanges <- function(SpanRanges) {
  DisjoinSpanRange <- IRanges::disjoin(SpanRanges)
  GenomicRanges::mcols(DisjoinSpanRange)$Class <- ""
  GenomicRanges::mcols(DisjoinSpanRange)$Title <- ""
  OverlapWithVariation <-
    data.frame(GenomicRanges::findOverlaps(DisjoinSpanRange, SpanRanges))
  TrashToAvoidPrint <- apply(OverlapWithVariation, 1, function(x) {
    GenomicRanges::mcols(DisjoinSpanRange[x["queryHits"]])$Class <<-
      paste0(
        GenomicRanges::mcols(DisjoinSpanRange[x["queryHits"]])$Class,
        " ",
        GenomicRanges::mcols(SpanRanges[x["subjectHits"]])$Class
      )
    GenomicRanges::mcols(DisjoinSpanRange[x["queryHits"]])$Title <<-
      paste(
        GenomicRanges::mcols(DisjoinSpanRange[x["queryHits"]])$Title,
        GenomicRanges::mcols(SpanRanges[x["subjectHits"]])$Title,
        sep = "<br>"
      )
  })
  return(DisjoinSpanRange)
}


#' SplitSpanRanges Split span ranges by insert Ranges. The span ranges will be break in two if they overlap an insert IRanges. Typically for the endline.
#'
#' @param SpanRanges Information range to be split if need be
#' @param InsertRanges Range of the position to split. Shoudl be of width 1
#' @import GenomicRanges
#' @export
#' @return SpanRanges with splitted
#'
#' @examples
#' library(GenomicRanges)
#' VariationRanges<-GetRangesFromUniprot("O94875-10")
#' Insert<-IRanges(60,60)
#' SplitSpanRanges(VariationRanges,Insert)
SplitSpanRanges <- function(SpanRanges, InsertRanges) {
  Overlaps <-
    data.frame(GenomicRanges::findOverlaps(SpanRanges, InsertRanges))
  while (dim(Overlaps)[1] > 0) {
    ##Case where the a peptide end with a new line it created a span that did not close anywhere and in doing so make it all
    if (end(SpanRanges[Overlaps$queryHits[1]]) == end(InsertRanges[Overlaps$subjectHits[1]])) {

    } else{
      Before <- SpanRanges[Overlaps$queryHits[1]]
      After <- Before
      end(Before) <- start(InsertRanges[Overlaps$subjectHits[1]])
      start(After) <-
        end(InsertRanges[Overlaps$subjectHits[1]]) + 1
      SpanRanges <- SpanRanges[-Overlaps$queryHits[1]]
      SpanRanges <- c(SpanRanges, Before, After)
      SpanRanges <- sort(SpanRanges)
    }
    InsertRanges <- InsertRanges[-Overlaps$subjectHits[1]]
    Overlaps <- data.frame(GenomicRanges::findOverlaps(SpanRanges, InsertRanges))
  }
  return(SpanRanges)
}
