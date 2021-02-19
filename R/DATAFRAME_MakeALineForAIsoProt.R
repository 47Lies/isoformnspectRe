#' Make a dataframe to be use in a future data table that will represent the protein, the peptides that match this protein, the other proteins matched with these peptides.
#' This function can handle the Canonical form and also the isoforms
#'
#'
#' @param LeadingRazorProtein Protein name as displayed in the MaxQuant seach file
#' @param DT_Peptides Peptide data table since we are using 'Leading razor protein' and not Leading.razor.protein
#' @param BS_ProteinBank AA StringSet protein bank.
#' @param ProteinWidth Length of the protein displayed width. expressed in AA
#'
#' @return Dataframe with every elements needed in the displayed data frame
#' Name name of the protein that will be displayed
#' N_Peptides Number of peptides that match this protein
#' PEPTIDES kable HTML ready to print table of the peptide that match this protein
#' N_Match_Proteins Number of differents proteins that match this set of peptides
#' Sequence HTML Ready to print sequence of the protein
#' @export
#' @rawNamespace import(data.table, except = shift)
#' @import knitr
#' @import kableExtra
#' @import DT
#'
#' @examples
#' library(data.table)
#' library(Biostrings)
#' library(DT)
#' PB_FilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' ProteinBank<-readAAStringSet(
#'   PB_FilePath
#' )
#' P_FilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' Peptides<-fread(
#'  P_FilePath,
#'  integer64="double")
#' names(ProteinBank) <- gsub(
#'  " .*$",
#'  "",
#'  names(ProteinBank))
#' LeadingRazorProtein <-
#' "strngt.5372.8,start=14,stop=904,frame=2,flag=1,UNIPROT=sp|Q86U42-2|PABP2_HUMAN"
#' DF<-DATAFRAME_MakeALineForAIsoProt(LeadingRazorProtein,Peptides,ProteinBank)
#' datatable(DF,escape=FALSE)
DATAFRAME_MakeALineForAIsoProt <- function(LeadingRazorProtein,
                                           DT_Peptides,
                                           BS_ProteinBank,
                                           ProteinWidth=60){
  PeptidesColsOfInterest <- NULL
  PeptidesColsOfInterest <-
    c("Start position", "End position", "PEP", "Proteins", "Sequence")
  #https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  PeptidesOfInterest <-
    DT_Peptides[DT_Peptides$`Leading razor protein` == LeadingRazorProtein,
                PeptidesColsOfInterest,
                with=FALSE]
  colnames(PeptidesOfInterest) <-
    c("Start", "End", "Pep", "Proteins", "Sequence")
  SequenceOfInterest <- as.vector(BS_ProteinBank[LeadingRazorProtein])
  names(SequenceOfInterest) <- NULL

  Name <- gsub(",", " ", LeadingRazorProtein)
  Name <- gsub("UNIPROT=", " ", Name)
  Name <- paste("<a href=\"https://www.uniprot.org/uniprot/",
                STRING_ProteinShortName(LeadingRazorProtein),
                "\" target=\"_blank\">",
                Name,
                "</a>",sep="")


  N_Peptides <- dim(PeptidesOfInterest)[1]

  List_Proteins <- strsplit(split = ";", x = PeptidesOfInterest$Proteins)
  List_Proteins <- lapply(List_Proteins, function(x) {
    unlist(lapply(x,
                  function(y) {
                    STRING_ProteinShortName(y)
                  }))
  })

  Total_MatchProteins <- unique(unlist(List_Proteins))

  Vector_Match <- unlist(lapply(List_Proteins,function(x){base::paste(x,collapse=" ")}))

  N_Match_Proteins <-length(Total_MatchProteins)

  VectorTotal_MatchProteins <- paste(Total_MatchProteins,collapse=" ")

  PeptidesOfInterest$Proteins <- Vector_Match
  PeptidesOfInterest$Pep <- base::formatC(PeptidesOfInterest$Pep, digits = 2)
  PeptidesOfInterest <- PeptidesOfInterest[order(PeptidesOfInterest$Start),]
  KABLE_PEPTIDES <-
    as.vector(
      kableExtra::scroll_box(
        kableExtra::kable_styling(
          knitr::kable(
            PeptidesOfInterest,
            row.names = FALSE)
          ),
        fixed_thead = T,
        height = "150px"
      )
    )
  Sequence <- HTML_Pretty_Sequence(
    PerfectIsoformName = LeadingRazorProtein,
    Sequence = SequenceOfInterest,
    Peptides = PeptidesOfInterest,
    LineWidth = ProteinWidth
  )
  return(data.frame(Name=Name,
                    N_Peptides=N_Peptides,
                    Match_Proteins=VectorTotal_MatchProteins,
                    N_Match_Proteins=N_Match_Proteins,
                    PEPTIDES=KABLE_PEPTIDES,
                    Sequence=Sequence))
}


#' Make a dataframe to be use in a future data table that will represent the protein, the peptides that match this protein, the other proteins matched with these peptides.
#' This function is dedicated to Blast forms, since it will compare two sequences and higlight differences.
#'
#' @param LeadingRazorProtein Protein name as displayed in the MaxQuant seach file
#' @param DT_Peptides Peptide data table since we are using 'Leading razor protein' and not Leading.razor.protein
#' @param BS_ProteinBank AA StringSet protein bank.
#' @param ProteinWidth Length of the protein displayed width. expressed in AA
#' @param Rescue Dataframe of a blast results
#' @return Dataframe with every elements needed in the displayed data frame
#' Name name of the protein that will be displayed
#' BlastHit name of the protein from the known protein bank that match the displayed protein
#' N_Peptides Number of peptides that match this protein
#' Match_Proteins Proteins match by this set of peptides
#' PEPTIDES kable HTML ready to print table of the peptide that match this protein
#' N_Match_Proteins Number of differents proteins that match this set of peptides
#' Sequence HTML Ready to print sequence of the protein
#' SAAV Boolean that indicate if this protein have an SAAV
#' SAAVxPep Boolean that indicate if this protein have peptide that display an SAAV
#' Insertion Boolean that indicate if this protein have an Insertion
#' InsertionxPep Boolean that indicate if this protein have a peptide that displayed a Insertion
#' Deletion Boolean that indicate if this protein have an deletion
#' DeletionxPep Boolean that indicate if this protein have a peptide that displayed a Deletion
#' OutOfBlast Boolean that indicate if this protein have a peptide that display a section out of the blast, i.e. a part that is specific
#' @export
#' @import knitr
#' @import kableExtra
#' @import DT
#'
#' @examples
#' library(data.table)
#' library(Biostrings)
#' library(DT)
#' PB_FilePath <- system.file("extdata",
#'  "protein_bank.example.fasta",
#'   package = "isoformnspectRe")
#' ProteinBank<-readAAStringSet(
#'   PB_FilePath
#' )
#' P_FilePath <- system.file("extdata",
#'  "peptides.example.txt",
#'   package = "isoformnspectRe")
#' Peptides<-fread(
#'  P_FilePath,
#'  integer64="double")
#' R_FilePath <- system.file("extdata",
#'  "blast_rescue.example.txt",
#'   package = "isoformnspectRe")
#' Rescue<-read.table(
#'  R_FilePath,
#'  header=TRUE)
#'
#' names(ProteinBank) <- gsub(
#'  " .*$",
#'  "",
#'  names(ProteinBank))
#' LeadingRazorProtein <-
#' "strngt.16858.3,start=259,stop=4446,frame=1,flag=LonguestORF,Origin=3',Blast=sp|Q15596|NCOA2_HUMAN"
#' DF <- DATAFRAME_MakeALineForABlastProt(
#' LeadingRazorProtein,
#' Peptides,
#' Rescue,
#' ProteinBank)
#' datatable(DF,escape=FALSE)
#' ToRepresent<-DF[,c(
#'  "Name",
#'  "BlastHit",
#'  "N_Peptides",
#'  "Match_Proteins",
#'  "N_Match_Proteins",
#'  "Sequence")]
#' datatable(ToRepresent,escape=FALSE)
#'
DATAFRAME_MakeALineForABlastProt<-function(
  LeadingRazorProtein,
  DT_Peptides,
  Rescue,
  BS_ProteinBank,
  ProteinWidth=60){
  SequenceOfInterest <- as.vector(BS_ProteinBank[LeadingRazorProtein])
  names(SequenceOfInterest) <- NULL
  QueryName<-gsub(",Blast.*$","",LeadingRazorProtein)
  MrnaName<-gsub(",.*$","",QueryName)
  LocalPep<-DT_Peptides[DT_Peptides$`Leading razor protein`==LeadingRazorProtein,]
  BlastHitName<-gsub("^.*Blast=","",LeadingRazorProtein)

  LocalRescue<-Rescue[Rescue$QuerySeqId==QueryName & Rescue$SubjectSeqId==BlastHitName,
                      c("QueryStart",
                        "QueryEnd",
                        "QuerySequence",
                        "SubjectSequence",
                        "SubjectSeqId",
                        "SubjectStart",
                        "SubjectEnd",
                        "SubjectLength")][1,]

  PeptidesColsOfInterest <- c("Start position", "End position", "PEP", "Proteins", "Sequence")
  PeptidesOfInterest <- DT_Peptides[DT_Peptides$`Leading razor protein` == LeadingRazorProtein,
                                    PeptidesColsOfInterest,
                                    with=FALSE]
  colnames(PeptidesOfInterest) <- c("Start", "End", "Pep", "Proteins", "Sequence")

  N_Peptides <- dim(PeptidesOfInterest)[1]

  List_Proteins <- strsplit(split = ";", x = PeptidesOfInterest$Proteins)
  List_Proteins <- lapply(List_Proteins, function(x) {
    unlist(lapply(x,
                  function(y) {
                    STRING_ProteinShortName(y)
                  }))
  })

  Total_MatchProteins <- unique(unlist(List_Proteins))

  Vector_Match <- unlist(
    lapply(
      List_Proteins,
      function(x){
        base::paste(x,collapse=" ")}))

  N_Match_Proteins <-length(Total_MatchProteins)

  VectorTotal_MatchProteins <- paste(Total_MatchProteins,collapse=" ")

  PeptidesOfInterest$Proteins <- Vector_Match
  PeptidesOfInterest$Pep <- formatC(PeptidesOfInterest$Pep, digits = 2)
  PeptidesOfInterest <- PeptidesOfInterest[order(PeptidesOfInterest$Start),]
  KABLE_PEPTIDES <-
    as.vector(
      kableExtra::scroll_box(
        kableExtra::kable_styling(
          knitr::kable(PeptidesOfInterest, row.names = FALSE)
          ),
        fixed_thead = T,
        height = "150px"))
  Sequence <- HTML_PrettyBlastSequence(
    Sequence = SequenceOfInterest,
    Peptides = PeptidesOfInterest,
    Rescue = LocalRescue,
    LineWidth = ProteinWidth
  )

  return(data.frame(Name=MrnaName,
                    BlastHit=BlastHitName,
                    N_Peptides=N_Peptides,
                    Match_Proteins=VectorTotal_MatchProteins,
                    N_Match_Proteins=N_Match_Proteins,
                    PEPTIDES=KABLE_PEPTIDES,
                    Sequence=Sequence$SequenceInHTML,
                    SAAV=Sequence$SAAV,
                    SAAVxPep=Sequence$SAAVxPep,
                    Insertion=Sequence$Insertion,
                    InsertionxPep=Sequence$InsertionxPep,
                    Deletion=Sequence$Deletion,
                    DeletionxPep=Sequence$DeletionxPep,
                    OutOfBlast=Sequence$OutOfBlast))
}
