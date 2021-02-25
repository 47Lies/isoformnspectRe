#' CompareTwoSequences Compare two AA sequences in a blast results
#'
#' @param QuerySequence The query sequence in the blast result file, the query sequence correspond to the sequence of the protein bank aligned to the knowned protein bank (Uniprot for example).
#' @param SubjectSequence The subject sequence in the blast results file, the subject sequence correspond to the known protein (typically one of UNIPROT) on which the query have been aligned.
#' @param QueryStart In the case where the matching part of the query do not beggin with its N-Term extremity.
#'
#' @return IRanges with metadata corresponding to the differences between the two sequences, SAAV, Insertion, Deletion
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' R_FilePath <- system.file("extdata",
#'  "blast_rescue.example.txt",
#'   package = "isoformnspectRe")
#' Rescue<-read.table(
#'  R_FilePath,
#'  header=TRUE,
#'  stringsAsFactors = FALSE)
#' CompareTwoSequences(
#'  as.vector(Rescue$QuerySequence[1]),
#'  as.vector(Rescue$SubjectSequence[1]),
#'  as.vector(Rescue$QueryStart[1]))
CompareTwoSequences <- function(QuerySequence,
                                SubjectSequence,
                                QueryStart=12) {
  QuerySequenceVector <- unlist(
    strsplit(
      split  = character(0),
      x = QuerySequence))
  SubjectSequenceVector <- unlist(
    strsplit(
      split  = character(0),
      x = SubjectSequence))
  DeletionsInQuery <- which(QuerySequenceVector == "-")
  Differences <- GetEmptySpanRanges()
  if (length(DeletionsInQuery) > 0) {
    #dirty trick every consecuttive -  will be merged
    DeletionsInQuery <- IRanges::reduce(
      IRanges::IRanges(
        DeletionsInQuery,
        DeletionsInQuery))
    #get the deleted sequence from the subject
    Deletions.Sequences <- substring(
      SubjectSequence,
      first = IRanges::start(DeletionsInQuery),
      last = IRanges::end(DeletionsInQuery)
    )
    #Each deletion will create a shift in the coordinate that we have to compensate
    ShiftInducedByDeletions <- c(
      0,
      cumsum(
        IRanges::width( DeletionsInQuery)))
    ShiftInducedByDeletions <-
      ShiftInducedByDeletions[-length(ShiftInducedByDeletions)]
    DeletionsInQuery <-
      IRanges::shift(DeletionsInQuery, -ShiftInducedByDeletions)
    #Deletion IRanges in protein coordinate
    EdgesOfTheDeletionsInProtein <-
      IRanges::IRanges(
        start = IRanges::start(DeletionsInQuery) - 1,
        end = IRanges::start(DeletionsInQuery))
    GenomicRanges::mcols(EdgesOfTheDeletionsInProtein) <- data.frame(
      Class = "Deletion",
      Title = paste(
        "Deletion:",
        Deletions.Sequences,
        sep = ""))
    EdgesOfTheDeletionsInProtein <-
      IRanges::shift(
        EdgesOfTheDeletionsInProtein,
        QueryStart - 1)
    Differences <- c(
      Differences,
      EdgesOfTheDeletionsInProtein)
    SubjectSequenceVector <-
      SubjectSequenceVector[-which(QuerySequenceVector == "-")]
    QuerySequenceVector <-
      QuerySequenceVector[-which(QuerySequenceVector == "-")]
  }

  # Much more straightforward for insertion
  #
  Insertions <- which(SubjectSequenceVector == "-") + QueryStart - 1
  if (length(Insertions) > 0) {
    Insertions.IR <- IRanges::reduce(
      IRanges::IRanges(
        Insertions,
        Insertions))
    GenomicRanges::mcols(Insertions.IR)$Class <- "Insertion"
    GenomicRanges::mcols(Insertions.IR)$Title <- "Insertion"
    #add to our vector of ranges
    Differences <- c(
      Differences,
      Insertions.IR)
    #Since we don't want to mess up insertion with SNV, we are doing to replace the gap in the subject sequences
    SubjectSequenceVector[which(SubjectSequenceVector == "-")] <-
      QuerySequenceVector[which(SubjectSequenceVector == "-")]

  }
  #Since we discard the gap from
  SubjectSequenceWhithoutGap.Sequence <-
    paste(SubjectSequenceVector, collapse = "")
  Diffs <- which(SubjectSequenceVector != QuerySequenceVector)
  if (length(Diffs) > 0) {
    Diffs.IR <- IRanges::reduce(
      IRanges::IRanges(
        Diffs,
        Diffs))
    Mismatch.Sequences <- substring(
      SubjectSequenceWhithoutGap.Sequence,
      first = IRanges::start(Diffs.IR),
      last = IRanges::end(Diffs.IR)
    )
    GenomicRanges::mcols(Diffs.IR)$Class <- "SAAV"
    GenomicRanges::mcols(Diffs.IR)$Title <- paste(
      "SAAV:",
      Mismatch.Sequences,
      sep = "")
    Differences <- c(
      Differences,
      Diffs.IR)
  }
  return(Differences)
}
