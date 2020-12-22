## -----------------------------------------------------------------------------
library(RIdeogram)#
head(karyotype_dual_comparison)

## -----------------------------------------------------------------------------
head(synteny_dual_comparison)

## -----------------------------------------------------------------------------
blast_RIdeogram <- function(df1,df2){
library(RIdeogram)#
library(magrittr)
suppressPackageStartupMessages(library("dplyr"))
library("rsvg")
df3<-df2[,c(3,7,8,9,10)]
df3$fill<-ifelse(df3$V3>90,"0080cc",
                 ifelse(df3$V3<80,"0ab276","e64e60"))
df3$Species_1<-1
df3$Species_2<-1
df4 <- df3 %>%
  select(Species_1,V7,V8,Species_2,V9,V10,fill)
colnames(df4)<-colnames(synteny_dual_comparison)
#ideogram(karyotype =df1 ,
         #synteny=df4,output = "blast_RIdeogram.svg")
#rsvg_pdf("blast_RIdeogram.svg","blast_RIdeogram.pdf")
}

## -----------------------------------------------------------------------------
df1<-data.frame(Chr=c("I","I"),
               Start=c(1,1),
               End=c(131478,444567),
               fill=c("FF9D1E","FF9D1E"),
              species=c("chloroplast","mitochondrion"),
              size=c(12,12),
              color=c(252525,252525))
df2<-data.frame(V1=c(rep("NC_044701",13)),
               V2=c(rep("NC_044768",13)),
               V3=c(96.76,83.69,85.28,74.07,90.54,95.56,96.51,96.43,94.94,96.77,87.34,100.00,96.97),
               V4=c(1603,423,326,891,148,90,86,84,79,62,79,31,33),
               V5=c(43,37,38,172,11,3,2,3,4,2,6,0,0),
               V6=c(7,19,5,45,3,1,1,0,0,0,3,0,1),
               V7=c(54477,68192,66335,102938,36180,49,110947,31738,53913,105594,88210,10365,66299),
               V8=c(56074,68587,66654,103801,36324,138,111032,31821,53991,105655,88284,10395,66330),
               V9=c(375799,44003,360145,332888,64144,406657,420751,384266,2926,110375,240028,63403,360093),
               V10=c(374201,44420,360466,333746,64291,406745,420667,384349,3004,110314,240106,63373,360125),
               V11=c(0.0,2e-101,1e-88,1e-83,4e-48,4e-33,2e-32,5e-32,2e-27,2e-21,2e-16,2e-07,2e-06),
               V12=c(2663,370,327,311,193,143,141,139,124,104,87.9,58.4,54.7))
blast_RIdeogram(df1,df2) 

## -----------------------------------------------------------------------------
df<-data.frame(chr=c(rep("chloroplast",2),rep("mitochondrial",2)),
               x=c(1,131478,1,444567),
               y=c(0,1,0,1))

df

## -----------------------------------------------------------------------------
df3 <-df2
df3

## -----------------------------------------------------------------------------
blast_circlize <- function(df,df1){
suppressPackageStartupMessages(library(circlize))#RIdeogram,magrittr,dplyr,rsvg,circlize,RColorBrewer,ComplexHeatmap
library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))
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

## -----------------------------------------------------------------------------
blast_circlize(df,df3)

