#' @title blast double sequence alignment with Rideogram
#' @description Use the R language package Rideogram to display the results of blast double sequence alignment
#' @param df1 
#' @param df2 
#' @return residual matrices
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' re <- Two_Fold_Cross_Validation(a,b)
#' print(re)
#' }
#' @export
blast_RIdeogram <- function(df1,df2){
  df3<-df2[,c(3,7,8,9,10)]
  df3$fill<-ifelse(df3$V3>90,"0080cc",
                   ifelse(df3$V3<80,"0ab276","e64e60"))
  df3$Species_1<-1
  df3$Species_2<-1
  df4 <- df3 %>%
    select(Species_1,V7,V8,Species_2,V9,V10,fill)
  colnames(df4)<-colnames(synteny_dual_comparison)
  ideogram(karyotype =df1 ,
           synteny=df4,output = "blast_RIdeogram.svg")
  rsvg_pdf("blast_RIdeogram.svg","blast_RIdeogram.pdf")
}
NULL