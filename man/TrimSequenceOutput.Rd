% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Server.R
\name{TrimSequenceOutput}
\alias{TrimSequenceOutput}
\title{substitute the peptide by a html span tag that will substitute a too long peptide}
\usage{
TrimSequenceOutput(Peptide)
}
\arguments{
\item{Peptide}{vector of peptides}
}
\value{
vector of peptides to be used in a DT::data.table will need escape=TRUE option
}
\description{
substitute the peptide by a html span tag that will substitute a too long peptide
}
\examples{
peptides<-c("FEZQFAZRFSQJDFCGFSDGHREZQGVSVFDSVGFDS",
"RETGERVGDF")
TrimSequenceOutput(peptides)
}
