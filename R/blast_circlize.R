
#' @title blast double sequence alignment results with circlize
#' @description Use the R language package circlize to visualize the blast double sequence alignment results
#' @param df 
#' @param df1 
#' @return a random sample of size
#' @export
blast_circlize <- function(df,df1){
  col<-RColorBrewer::brewer.pal(6,"Paired")
  circos.par("start.degree" = 130)
  circos.initialize(factors = df$chr,x=df$x)
  circos.trackPlotRegion(factors = df$chr,y=df$y,
                         panel.fun = function(x,y){
                           circos.axis()
                         },track.height = 0.1)
  highlight.sector(sector.index = "chloroplast",col=col[1])
  highlight.sector(sector.index = "mitochondrial",col=col[2])
  circos.text(x=70000,y=0.5,
              labels = "chloroplast",
              sector.index = "chloroplast")
  circos.text(x=220000,y=0.5,
              labels = "mitochondrial",
              sector.index = "mitochondrial",
              facing = "outside")
  col_fun = colorRamp2(c(70,90,100),
                       c("green", "yellow", "red"))
  for (i in 1:13){
    x<-sort(c(df1[i,8],df1[i,7]))
    y<-sort(c(df1[i,10],df1[i,9]))
    z<-df1[i,3]
    circos.link("chloroplast",x,"mitochondrial",y,
                col=add_transparency(col_fun(z)))
  }
  circos.clear()
  lgd_links = Legend(at = c(70, 80, 90, 100), 
                     col_fun = col_fun, 
                     title_position = "topleft",
                     title = "identity(%)")
  lgd_list_vertical = packLegend(lgd_links)
  
  draw(lgd_list_vertical, x = unit(10, "mm"), 
       y = unit(10, "mm"), just = c("left", "bottom"))
}
NULL