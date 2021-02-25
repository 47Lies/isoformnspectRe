# isoformnspectRe
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/47Lies/isoformnspectRe/workflows/R-CMD-check/badge.svg)](https://github.com/47Lies/isoformnspectRe/actions)
  [![R-CMD-check](https://github.com/47Lies/isoformnspectRe/actions/workflows/r.yml/badge.svg)](https://github.com/47Lies/isoformnspectRe/actions/workflows/r.yml)
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

## Render rmarkdown



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
ProteinBankFastaFilePath <- system.file("extdata", "protein_bank.example.fasta", package = "isoformnspectRe")
PeptidePath <- system.file("extdata", "peptides.example.txt", package = "isoformnspectRe")
SampleDescriptionPath <- system.file("extdata", "SampleDescription.txt", package = "isoformnspectRe")
IntensityName <- "Intensity "
SampleColumnName <- "SampleName"
SampleGroupColumnName <- "LooseSampleGroup"
shiny::shinyApp(UI,server)
```

