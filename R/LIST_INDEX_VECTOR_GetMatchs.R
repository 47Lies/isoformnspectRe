#' LIST_INDEX_VECTOR_GetMatchs
#' get the matchings index un the protein bank
#' @param Peptides vector of peptides in amino acid code (capital letters)
#' @param ProteinBank AAStringSet with their names
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
#' LIST_INDEX_VECTOR_GetMatchs(Peptides,Proteins)
LIST_INDEX_VECTOR_GetMatchs<-function(Peptides,ProteinBank){
  x<-as.vector(Peptides)
  #x<-c("rko","QTASLTAAYGQLSK","ELSIEMQDKDCQEASGHLESR","TSEFHLGLIEGPDKNK")
  x<-gsub("I|L","(I|L)",x)
  x<-paste("(^|M|K|R)",x,"([^P]|$)",sep="")
  if(future::availableCores()>1){
    future::plan(future::multisession)
    options(future.globals.maxSize=4*1024*1024**2)
    MatchsIndex<-furrr::future_map(
      x,
      ~stringr::str_which(
        string=as.vector(ProteinBank),
        pattern=.x),
      .progress = TRUE)
    future::plan(future::sequential)
  }else{
    MatchsIndex<-lapply(x,
                        function(x){
                          stringr::str_which(
                            string=as.vector(ProteinBank),
                            pattern=x)
                        }
    )
  }
  return(MatchsIndex)
}


#' DATAFRAME_AnnotateProteotypic Since we observe that some proteins are not present in the Proteins field, we decieded to
#'
#' @param PeptideDataFrame MaxQuant peptide search as a data frame
#' @param ProteinBank AAStringSet of the protein fasta file that was used for the mass spectrometry analysis
#'
#' @return PeptideDataFrame Updated data frame with informations
#' @export
DATAFRAME_AnnotateProteotypic<-function(PeptideDataFrame,ProteinBank){
  Matchs<-LIST_INDEX_VECTOR_GetMatchs(PeptideDataFrame$Sequence,ProteinBank)
  PeptideDataFrame$UpdateProteins<-unlist(
    lapply(
      Matchs,
      function(x){
        paste(names(ProteinBank)[x],collapse=";")
      }
    )
  )
  PeptideDataFrame$NProteins<-unlist(lapply(Matchs,length))
  PeptideDataFrame$Old_Proteins<-PeptideDataFrame$Proteins
  PeptideDataFrame$Proteins<-PeptideDataFrame$UpdateProteins
  return(PeptideDataFrame)
}
