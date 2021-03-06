% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetRangesFromUniprot.R
\name{GetRangesFromUniprot}
\alias{GetRangesFromUniprot}
\title{GetRangesFromUniprot From a Uniprot isoform entry get the variations from the canonical sequence. The start and end coordonate are isoform sequence based.}
\usage{
GetRangesFromUniprot(UNIPROT_ISOFORM_ENTRY = "O94875-10")
}
\arguments{
\item{UNIPROT_ISOFORM_ENTRY}{Isoform entry name typically a capital letter}
}
\value{
IRanges object with each line representing the variant splicing of the Isoform. The class and the Title metadata columns are ready for a HTML span inclusion.
}
\description{
GetRangesFromUniprot From a Uniprot isoform entry get the variations from the canonical sequence. The start and end coordonate are isoform sequence based.
}
\examples{
VariationRanges<-GetRangesFromUniprot("O94875-10")
VariationRanges
}
