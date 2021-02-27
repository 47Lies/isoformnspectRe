#' LIST_INDEX_VECTOR_GetMatchs
#' get the matchings index un the protein bank
#' @param Peptides vector of peptides in amino acid code (capital letters)
#' @param ProteinBank AAStringSet with their names
#' @param Threads Number of parallel sessions (numeric)
#'
#' @return LIST_INDEX_VECTOR_GetMatchs list of index vector each index of the list correspond to the index in the params Peptides vector. In each list element the vector correspond to the index of the protein that match this very particular pepetide.
#' @export
#' @import Biostrings
#' @import furrr
#' @importFrom future availableCores
#' @importFrom future plan
#' @importFrom future multisession
#' @importFrom future sequential
#' @import stringr
#' @examples
#' library(Biostrings)
#' Proteins<-c("MARFLTLCTWLLLLGPGLLATVRAECSQDCATCSYRLVRPADINFLACV",
#' "MECEGKLPSLKIWETCKELLQLSKPELPQDGTSTLRENSKPEESHLLAKRYGGF",
#' "MKRYGGFMKKMDELYPMEPEEEANGSEILAKRYGGF",
#' "MKKDAEEDDSLANSSDLLKELLETGDNRERSHHQDGSDNEEEVSKRYGGF",
#' "MRGLKRSPQLEDEAKELQKRYGGFMRRVGRPEWWMDYQKRYGGFLKRFAEALPSDEEGESYSKEVPEMEKRYGGFMRF")
#' names(Proteins) <-
#' c("pseudo prot A",
#'   "pseudo prot B",
#'   "pseudo prot B",
#'   "pseudo prot C",
#'   "pseudo prot D")
#' Proteins<-AAStringSet(Proteins)
#' Peptides<-c("RRRR","MK","PQDGTSTLRENS")
#' LIST_INDEX_VECTOR_GetMatchs(Peptides,Proteins,3)
LIST_INDEX_VECTOR_GetMatchs <-
  function(Peptides, ProteinBank, Threads = 1) {
    x <- as.vector(Peptides)
    #x<-c("rko","QTASLTAAYGQLSK","ELSIEMQDKDCQEASGHLESR","TSEFHLGLIEGPDKNK")
    x <- gsub("I|L", "(I|L)", x)
    x <- paste("(^|M|K|R)", x, "([^P]|$)", sep = "")
    wrkrs = Threads
    if (availableCores() > 1) {
      if (Threads > availableCores()) {
        wrkrs = availableCores()
      }
      future::plan(multisession, workers = wrkrs)
      options(future.globals.maxSize = 4 * 1024 * 1024 ** 2)
      MatchsIndex <- furrr::future_map(x,
                                       ~ stringr::str_which(string = as.vector(ProteinBank),
                                                            pattern = .x),
                                       .progress = TRUE)
      future::plan(future::sequential)
    } else{
      MatchsIndex <- lapply(x,
                            function(x) {
                              stringr::str_which(string = as.vector(ProteinBank),
                                                 pattern = x)
                            })
    }
    return(MatchsIndex)
  }


#' DATAFRAME_AnnotateProteotypic Since we observe that some proteins are not present in the Proteins field, we decieded to
#'
#' @param PeptideDataFrame MaxQuant peptide search as a data frame
#' @param ProteinBank AAStringSet of the protein fasta file that was used for the mass spectrometry analysis
#' @param Threads Number of parallel sessions (numeric)
#'
#' @return PeptideDataFrame Updated data frame with informations
#' @import tidyr
#' @rawNamespace import(data.table,except=c(shift))
#' @export
#' @examples 
#' library(isoformnspectRe)
#' ProteinBankFastaFilePath <- system.file("extdata",
#' "protein_bank.example.fasta",
#' package = "isoformnspectRe")
#' ProteinBank <- Biostrings::readAAStringSet(ProteinBankFastaFilePath)
#' PeptidePath <- system.file("extdata",
#'                            "peptides.example.txt",
#'                            package = "isoformnspectRe")
#' options(datatable.integer64 = "numeric")
#' Peptides <- data.table::fread(PeptidePath, sep = "\t")
#' DATATABLE_AnnotateProteotypic(Peptides, ProteinBank, Threads = 1)
DATATABLE_AnnotateProteotypic <-
  function(PeptideDataFrame, ProteinBank, Threads = 1) {
    Matchs <-
      LIST_INDEX_VECTOR_GetMatchs(PeptideDataFrame$Sequence, ProteinBank, Threads)
    PeptideDataFrame$UpdateProteins <- unlist(lapply(Matchs,
                                                     function(x) {
                                                       paste(names(ProteinBank)[x], collapse = ";")
                                                     }))
    PeptideDataFrame$NProteins <- unlist(lapply(Matchs, length))
    PeptideDataFrame$Old_Proteins <- PeptideDataFrame$Proteins
    PeptideDataFrame$Proteins <- PeptideDataFrame$UpdateProteins
    
    PeptideDataFrame$SplitedUpdateProteins <- PeptideDataFrame$UpdateProteins
    if(sum("data.table" %in% class(PeptideDataFrame))==0){
      PeptideDataFrame<-data.table::data.table(PeptideDataFrame)
    }
    #as many lines as proteins per peptides
    Peptide <-
      data.table::data.table(tidyr::unnest(
        PeptideDataFrame,
        SplitedUpdateProteins = strsplit(PeptideDataFrame$SplitedUpdateProteins, ";")
      ))
    #keeps only the peptides that match a leading razor proteins
    PeptideDataFrame <-
      PeptideDataFrame[PeptideDataFrame$SplitedUpdateProteins %in% PeptideDataFrame$`Leading razor protein`,]
    #Get localisation
    Positions <-
      matrix(stringr::str_locate(ProteinBank[PeptideDataFrame$SplitedUpdateProteins], gsub("I|L", "(I|L)", PeptideDataFrame$Sequence)), ncol =
               2)
    PeptideDataFrame[, c("Start position")] <- Positions[, 1]
    PeptideDataFrame[, c("End position")] <- Positions[, 2]
    PeptideDataFrame$"UpdatedSequence" <-
      substr(ProteinBank[PeptideDataFrame$SplitedUpdateProteins], PeptideDataFrame$"Start position", PeptideDataFrame$"End position")
    
    return(PeptideDataFrame)
  }
