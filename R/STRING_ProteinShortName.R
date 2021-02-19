#' STRING_ProteinShortName Handle the name of a protein as presented in the Maxquant peptide file, to be presented in the kable
#'
#' @param ProteinFullName name as observed in the MaxQuant peptide file
#'
#' @return short name with the UNIPROT code if present. b is for blast equivalent.
#' @export
#'
#' @examples
#' LeadingRazorProtein <-
#' "strngt.5372.8,start=14,stop=904,frame=2,flag=1,UNIPROT=sp|Q86U42-2|PABP2_HUMAN"
#' STRING_ProteinShortName(LeadingRazorProtein)
STRING_ProteinShortName<-function(ProteinFullName){
  if(length(grep("UNIPROT=",ProteinFullName))==1){
    ShortName<-gsub("^.*UNIPROT=[a-z]+\\|([^\\|]*)\\|.*$","\\1",ProteinFullName)
  }else if(length(grep("Blast=",ProteinFullName))==1){
    ShortName<-gsub("^.*Blast=[a-z]+\\|([^\\|]*)\\|.*$","\\1",ProteinFullName)
    ShortName<-paste("b:",ShortName,sep="")
  }else if(length(grep("^str",ProteinFullName))==1){
    ShortName<-gsub("^([^,]),.*$","\\1",ProteinFullName)
  }else if(length(grep("^(sp|tr)",ProteinFullName))==1){
    ShortName<-gsub("^[^\\|]*\\|([^\\|]*)\\|.*$","\\1",ProteinFullName)
  }else{
    ShortName<-ProteinFullName
  }
  return(ShortName)
}
