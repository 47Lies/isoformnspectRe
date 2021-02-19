#' Title
#'
#' @param RegularProt Name of the protein as displayed in the peptide file, most likely the first word(no space) of the fasta protein bank name.
#' @param ProteinBankFastaFilePath Path to the fasta file that have been used for this analysis
#' @param PepFilePath Path to the peptide file
#' @param SampleDescriptionFilePath Path to the sample description file
#' @param IntensityPrefix Prefix of the intensity, keep the space in case of
#' @param SampleGroupName Name of the Sample group column in the Sample description file
#' @param SampleName Name of the Sample column in the Sample description file
#' @param OutputDir Name of the output directory, the regular protein html file will be created in the isoform skelton subdirectory
#'
#' @export
#' @import rmarkdown
#' @import flexdashboard
#' @import GenomicFeatures
#' @import sessioninfo
#' @import sparklines
#' @import viridisLite
#'
#' @examples
#' library(rmarkdown)
#' library("flexdashboard")
#' library("sessioninfo")
#' ProteinBankFastaFilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' PepFilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' IntensityPrefix <- "Intensity "
#' Pep<-read.table(PepFilePath, header=TRUE,sep=",",quote="\"")
#' ProteinOfInterest <- unique(Pep$Leading.razor.protein)
#' ProteinOfInterest <- ProteinOfInterest[grep("UNIPROT=",ProteinOfInterest)]
#' ProteinOfInterest <- ProteinOfInterest[grep("-",ProteinOfInterest,invert=TRUE)][1]
#' SampleDescriptionFilePath <- system.file("extdata",
#'  "SampleDescription.txt",
#'   package = "isoformnspectRe")
#' SampleGroupName <- "LooseSampleGroup"
#' SampleName <- "SampleName"
#' RenderARegular(RegularProt=ProteinOfInterest,
#'  ProteinBankFastaFilePath=ProteinBankFastaFilePath,
#'  PepFilePath=PepFilePath,
#'  SampleDescriptionFilePath=SampleDescriptionFilePath,
#'  IntensityPrefix=IntensityPrefix,
#'  SampleGroupName=SampleGroupName,
#'  SampleName=SampleName,
#'  OutputDir=tempdir())
RenderARegular<-function(RegularProt,
                         ProteinBankFastaFilePath,
                         PepFilePath,
                         SampleDescriptionFilePath,
                         IntensityPrefix,
                         SampleGroupName,
                         SampleName,
                         OutputDir="./"){
  FileName <- gsub("[[:punct:]]",
                   "_",
                   RegularProt)
  Rmd <- system.file("rmd",
                     "RegularProteinSkeleton.V2.Rmd",
                     package = "isoformnspectRe")
  Css <- system.file("extdata",
                     "IsoAndSpe.css",
                     package = "isoformnspectRe")
  Logo <- system.file("extdata",
                     "IsoAndSpe.png",
                     package = "isoformnspectRe")
  Temp<-tempdir()
  if (!file.exists(paste0(OutputDir, "/RegularSkeleton/", FileName, ".html", sep = ""))) {
      rmarkdown::render(
      input = Rmd,
      output_file = paste0(FileName, ".html"),
      output_format = flex_dashboard(
        css = Css,
        logo = Logo,
        orientation = "columns"
      ),
      intermediates_dir = Temp,
      output_dir = paste0(OutputDir, "/RegularSkeleton"),
      params = list(
        ProteinBankFasta = ProteinBankFastaFilePath,
        PeptideResults = PepFilePath,
        ProteinOfInterest = RegularProt,
        SampleDescription =  SampleDescriptionFilePath,
        IntensityPrefix = IntensityPrefix,
        SampleGroupName = SampleGroupName,
        SampleName = SampleName),
      clean = TRUE,
      quiet = TRUE
    )
  }
  unlink(Temp)
  return(0)
}

#' Title
#'
#' @param IsoformProt Name of the protein as displayed in the peptide file, most likely the first word(no space) of the fasta protein bank name.
#' @param ProteinBankFastaFilePath Path to the fasta file that have been used for this analysis
#' @param PepFilePath Path to the peptide file
#' @param SampleDescriptionFilePath Path to the sample description file
#' @param IntensityPrefix Prefix of the intensity, keep the space in case of
#' @param SampleGroupName Name of the Sample group column in the Sample description file
#' @param SampleName Name of the Sample column in the Sample description file
#' @param OutputDir Name of the output directory, the isoform html file will be created in the isoform skelton subdirectory
#'
#' @export
#' @import rmarkdown
#' @import flexdashboard
#' @import GenomicFeatures
#' @import sessioninfo
#' @import sparklines
#' @import viridisLite
#'
#' @examples
#' library(rmarkdown)
#' library("flexdashboard")
#' library("sessioninfo")
#' ProteinBankFastaFilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' PepFilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' IntensityPrefix <- "Intensity "
#' Pep<-read.table(PepFilePath, header=TRUE,sep=",",quote="\"")
#' ProteinOfInterest <- unique(Pep$Leading.razor.protein)
#' ProteinOfInterest <- ProteinOfInterest[grep("UNIPROT=",ProteinOfInterest)]
#' ProteinOfInterest <- ProteinOfInterest[grep("-",ProteinOfInterest)][1]
#' SampleDescriptionFilePath <- system.file("extdata",
#'  "SampleDescription.txt",
#'   package = "isoformnspectRe")
#' SampleGroupName <- "LooseSampleGroup"
#' SampleName <- "SampleName"
#' RenderAnIsoform(IsoformProt=ProteinOfInterest,
#'  ProteinBankFastaFilePath=ProteinBankFastaFilePath,
#'  PepFilePath=PepFilePath,
#'  SampleDescriptionFilePath=SampleDescriptionFilePath,
#'  IntensityPrefix=IntensityPrefix,
#'  SampleGroupName=SampleGroupName,
#'  SampleName=SampleName,
#'  OutputDir=tempdir())
RenderAnIsoform<-function(IsoformProt,
                          ProteinBankFastaFilePath,
                          PepFilePath,
                          SampleDescriptionFilePath,
                          IntensityPrefix,
                          SampleGroupName,
                          SampleName,
                          OutputDir="./"){
  FileName <- gsub("[[:punct:]]",
                   "_",
                   IsoformProt)
  Rmd <- system.file("rmd",
                     "RegularProteinSkeleton.V2.Rmd",
                     package = "isoformnspectRe")
  Css <- system.file("extdata",
                     "IsoAndSpe.css",
                     package = "isoformnspectRe")
  Logo <- system.file("extdata",
                      "IsoAndSpe.png",
                      package = "isoformnspectRe")
  Temp<-tempdir()
  if (!file.exists(paste0(OutputDir, "/IsoformProtein/", FileName, ".html"))) {
    rmarkdown::render(
      input = Rmd,
      output_file = paste0(FileName, ".html"),
      output_format = flex_dashboard(
        css = Css,
        logo = Logo,
        orientation = "columns"
      ),
      intermediates_dir = Temp,
      output_dir = paste0(OutputDir, "/IsoformProtein"),
      params = list(
        ProteinBankFasta = ProteinBankFastaFilePath,
        PeptideResults = PepFilePath,
        ProteinOfInterest = IsoformProt,
        SampleDescription = SampleDescriptionFilePath,
        IntensityPrefix=IntensityPrefix,
        SampleGroupName=SampleGroupName,
        SampleName=SampleName
      ),
      clean = TRUE,
      quiet = TRUE
    )
  }
  unlink(Temp)
  return(0)
}


#' Title
#'
#' @param ProteinBankFastaFilePath Path to the fasta file that have been used for this analysis
#' @param PepFilePath Path to the peptide file
#' @param SampleDescriptionFilePath Path to the sample description file
#' @param IntensityPrefix Prefix of the intensity, keep the space in case of
#' @param SampleGroupName Name of the Sample group column in the Sample description file
#' @param SampleName Name of the Sample column in the Sample description file
#' @param OutputDir Name of the output directory, the isoform html file will be created in the isoform skelton subdirectory
#' @param BlastProt Name of the protein of interest
#' @param BlastFilePath Path to the Blast result file
#'
#' @export
#' @import rmarkdown
#' @import flexdashboard
#' @import GenomicFeatures
#' @import sessioninfo
#' @import sparklines
#' @import viridisLite
#'
#' @examples
#' library(rmarkdown)
#' library("flexdashboard")
#' library("sessioninfo")
#' ProteinBankFastaFilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' PepFilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' BlastFilePath <- system.file("extdata",
#'  "blast_rescue.example.txt",
#'   package = "isoformnspectRe")
#' IntensityPrefix <- "Intensity "
#' Pep<-read.table(PepFilePath, header=TRUE,sep=",",quote="\"")
#' ProteinOfInterest <- unique(Pep$Leading.razor.protein)
#' ProteinOfInterest <- ProteinOfInterest[grep("Blast=",ProteinOfInterest)][1]
#'
#' SampleDescriptionFilePath <- system.file("extdata",
#'  "SampleDescription.txt",
#'   package = "isoformnspectRe")
#' SampleGroupName <- "LooseSampleGroup"
#' SampleName <- "SampleName"
#' RenderABlast(BlastProt=ProteinOfInterest,
#'  ProteinBankFastaFilePath=ProteinBankFastaFilePath,
#'  PepFilePath=PepFilePath,
#'  SampleDescriptionFilePath=SampleDescriptionFilePath,
#'  IntensityPrefix=IntensityPrefix,
#'  SampleGroupName=SampleGroupName,
#'  BlastFilePath=BlastFilePath,
#'  SampleName=SampleName,
#'  OutputDir=tempdir())
RenderABlast<-function(BlastProt,
                       ProteinBankFastaFilePath,
                       PepFilePath,
                       SampleDescriptionFilePath,
                       IntensityPrefix,
                       SampleGroupName,
                       BlastFilePath,
                       SampleName,
                       OutputDir="./"){
  FileName <- gsub("[[:punct:]]",
                   "_",
                   BlastProt)
  Rmd <- system.file("rmd",
                     "BlastSkeleton.V2.Rmd",
                     package = "isoformnspectRe")
  Css <- system.file("extdata",
                     "IsoAndSpe.css",
                     package = "isoformnspectRe")
  Logo <- system.file("extdata",
                      "IsoAndSpe.png",
                      package = "isoformnspectRe")
  Temp<-tempdir()
  if(!file.exists(paste0(OutputDir, "/BlastProtein/", FileName, ".html"))) {
    rmarkdown::render(
      input = Rmd,
      output_file = paste0(FileName, ".html"),
      output_format = flex_dashboard(
        css = Css,
        logo = Logo,
        orientation = "columns"
      ),
      intermediates_dir = Temp,
      output_dir = paste0(OutputDir, "/BlastProtein"),
      params = list(
        ProteinBankFasta = ProteinBankFastaFilePath,
        PeptideResults = PepFilePath,
        ProteinOfInterest = BlastProt,
        SampleDescription =  SampleDescriptionFilePath,
        Blast = BlastFilePath,
        IntensityPrefix=IntensityPrefix,
        SampleGroupName=SampleGroupName,
        SampleName=SampleName
      ),
      clean = TRUE,
      quiet = TRUE
    )
  }
  return(0)

}
