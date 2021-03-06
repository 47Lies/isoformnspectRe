% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Render.R
\name{RenderAnIsoform}
\alias{RenderAnIsoform}
\title{Title}
\usage{
RenderAnIsoform(
  IsoformProt,
  ProteinBankFastaFilePath,
  PepFilePath,
  SampleDescriptionFilePath,
  IntensityPrefix,
  SampleGroupName,
  SampleName,
  OutputDir = "./"
)
}
\arguments{
\item{IsoformProt}{Name of the protein as displayed in the peptide file, most likely the first word(no space) of the fasta protein bank name.}

\item{ProteinBankFastaFilePath}{Path to the fasta file that have been used for this analysis}

\item{PepFilePath}{Path to the peptide file}

\item{SampleDescriptionFilePath}{Path to the sample description file}

\item{IntensityPrefix}{Prefix of the intensity, keep the space in case of}

\item{SampleGroupName}{Name of the Sample group column in the Sample description file}

\item{SampleName}{Name of the Sample column in the Sample description file}

\item{OutputDir}{Name of the output directory, the isoform html file will be created in the isoform skelton subdirectory}
}
\description{
Title
}
\examples{
library(rmarkdown)
library("flexdashboard")
library("sessioninfo")
library("sparkline")
ProteinBankFastaFilePath <- system.file("extdata",
 "protein_bank.example.fasta",
  package = "isoformnspectRe")
PepFilePath <- system.file("extdata",
 "peptides.example.txt",
  package = "isoformnspectRe")
IntensityPrefix <- "Intensity "
Pep<-read.table(PepFilePath, header=TRUE,sep="\t",quote="\"")
ProteinOfInterest <- unique(Pep$Leading.razor.protein)
ProteinOfInterest <- ProteinOfInterest[grep("UNIPROT=",ProteinOfInterest)]
ProteinOfInterest <- ProteinOfInterest[grep("-",ProteinOfInterest)][1]
SampleDescriptionFilePath <- system.file("extdata",
 "SampleDescription.txt",
  package = "isoformnspectRe")
SampleGroupName <- "LooseSampleGroup"
SampleName <- "SampleName"
RenderAnIsoform(IsoformProt=ProteinOfInterest,
 ProteinBankFastaFilePath=ProteinBankFastaFilePath,
 PepFilePath=PepFilePath,
 SampleDescriptionFilePath=SampleDescriptionFilePath,
 IntensityPrefix=IntensityPrefix,
 SampleGroupName=SampleGroupName,
 SampleName=SampleName,
 OutputDir=tempdir())
}
