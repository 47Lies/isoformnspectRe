# isoformnspectRe
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/47Lies/isoformnspectRe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/47Lies/isoformnspectRe/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->


Companion package of the Iso-Spe pipeline

The package aims to handle the R dependancies of the Iso-Spe pipeline. The Nextflow philosophy is to not have any absolute path in its script and we needed to have some Rmarkdown somewhere, the R package seems to be the best solution.

# Installation

Install the developpement version from GitHub with:

```
devtools::install_github("47Lies/isoformnspectRe")
```

If you use renv package manager, be sure to have devtools installed before:

```
renv::install("47Lies/isoformnspectRe")
```

To update from renv, you may need to delete the previously installed version due to its presence in the cache and re install it:
```
renv::remove("isoformnspectRe")
renv::install("47Lies/isoformnspectRe")
```

# Usage

## Update peptides

When using MaxQuant, we observed that some peptides were wrongly attributed. For example, some peptides were solely attributed to the Blast protein (the one that was found from the mRNA and look like a protein from UNIPROT) form while it should have been also attributed to the Original(the protein from UNIPROT previously mentionned). We came with a solution to make a new search of the peptide in the protein bank.

The function have a multisession capacity
```
library("isoformnspectRe")
ProteinBankFastaFilePath <- system.file("extdata",
  "protein_bank.example.fasta",
  package = "isoformnspectRe")
PeptidePath <- system.file("extdata",
  "peptides.example.txt",
  package = "isoformnspectRe")
options(datatable.integer64 = "numeric")
Peptide <- data.table::fread(opt$Peptides, sep = "\t")
DATATABLE_AnnotateProteotypic(PeptideDataFrame, ProteinBank, Threads = 1)
```



## Render rmarkdown

Based on a files obtained through the Iso and Spe pipeline. There is 3 archetypes of markdowns if the mrna match an isoform, have a blast hit or is a perfect match.

```
library("isoformnspectRe")
ProteinBankFastaFilePath <- system.file("extdata", "protein_bank.example.fasta", package = "isoformnspectRe")
PepFilePath <- system.file("extdata", "peptides.example.txt", package = "isoformnspectRe")
IntensityPrefix <- "Intensity "
Pep<-read.table(PepFilePath, header=TRUE,sep="\t", quote="\"")
ProteinOfInterest <- unique(Pep$Leading.razor.protein)
ProteinOfInterest <- ProteinOfInterest[grep("UNIPROT=",ProteinOfInterest)]
ProteinOfInterest <- ProteinOfInterest[grep("-",ProteinOfInterest,invert=TRUE)][1]
SampleDescriptionFilePath <- system.file("extdata",  "SampleDescription.txt", package = "isoformnspectRe")
SampleGroupName <- "LooseSampleGroup"
SampleName <- "SampleName"
RenderARegular(RegularProt=ProteinOfInterest,
ProteinBankFastaFilePath=ProteinBankFastaFilePath,
PepFilePath=PepFilePath,
SampleDescriptionFilePath=SampleDescriptionFilePath,
IntensityPrefix=IntensityPrefix,
SampleGroupName=SampleGroupName,
SampleName=SampleName,
OutputDir=tempdir())
```

## ShinyApp

The included shiny app aims to provide a navigation tool across the mass spec MaxQuant peptides and have html link to each of the rmarkdown of the proteins.

You will need to provide some informations before being able to launch it:

**ProteinBankFastaFilePath** path to the protein bank fasta file produced by the Iso-Spe pipeline

**PeptidePath** path to the peptide file produced by the MaxQuant search engine, can be either the regular peptide file or the updated version. We recommand to use the updated version.

**IntensityName** the used intensity prefix in the maxquant peptide file

**SampleDescriptionPath** path to the sample description file. You might have some trouble with names using '-'.

**SampleColumnName** Column name for the sample in the sample description file

**SampleGroupColumnName** Column name for the sample group in the sample description file

```
library("isoformnspectRe")
ProteinBankFastaFilePath <- system.file("extdata",
  "protein_bank.example.fasta",
  package = "isoformnspectRe")
PeptidePath <- system.file("extdata",
  "peptides.example.txt",
  package = "isoformnspectRe")
SampleDescriptionPath <- system.file("extdata",
  "SampleDescription.txt",
  package = "isoformnspectRe")
IntensityName <- "Intensity "
SampleColumnName <- "SampleName"
SampleGroupColumnName <- "LooseSampleGroup"
shiny::shinyApp(UI,server)
```

