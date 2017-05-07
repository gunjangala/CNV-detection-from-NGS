library(shape)

setwd("/home/ggg256/RD")
sample_1 <- read.delim("sample_01_RD.txt", header=FALSE)
sample_2 <- read.delim("sample_02_RD.txt", header=FALSE)
sample_3 <- read.delim("sample_03_RD.txt", header=FALSE)
sample_4 <- read.delim("sample_04_RD.txt", header=FALSE)
sample_5 <- read.delim("sample_05_RD.txt", header=FALSE)
sample_6 <- read.delim("sample_06_RD.txt", header=FALSE)
sample_7 <- read.delim("sample_07_RD.txt", header=FALSE)
sample_8 <- read.delim("sample_08_RD.txt", header=FALSE)
sample_9 <- read.delim("sample_09_RD.txt", header=FALSE)
sample_10 <- read.delim("sample_10_RD.txt", header=FALSE)

sample_11 <- read.delim("sample_11_RD.txt", header=FALSE)
sample_12 <- read.delim("sample_12_RD.txt", header=FALSE)
sample_13 <- read.delim("sample_13_RD.txt", header=FALSE)
sample_14 <- read.delim("sample_14_RD.txt", header=FALSE)
sample_15 <- read.delim("sample_15_RD.txt", header=FALSE)
sample_16 <- read.delim("sample_16_RD.txt", header=FALSE)
sample_17 <- read.delim("sample_17_RD.txt", header=FALSE)
sample_18 <- read.delim("sample_18_RD.txt", header=FALSE)
sample_19 <- read.delim("sample_19_RD.txt", header=FALSE)
sample_20 <- read.delim("sample_20_RD.txt", header=FALSE)

sample_21 <- read.delim("sample_21_RD.txt", header=FALSE)
sample_22 <- read.delim("sample_22_RD.txt", header=FALSE)
sample_23 <- read.delim("sample_23_RD.txt", header=FALSE)
sample_24 <- read.delim("sample_24_RD.txt", header=FALSE)
sample_25 <- read.delim("sample_25_RD.txt", header=FALSE)
sample_26 <- read.delim("sample_26_RD.txt", header=FALSE)
sample_27 <- read.delim("sample_27_RD.txt", header=FALSE)
sample_28 <- read.delim("sample_28_RD.txt", header=FALSE)
sample_29 <- read.delim("sample_29_RD.txt", header=FALSE)
sample_30 <- read.delim("sample_30_RD.txt", header=FALSE)

sample_31 <- read.delim("sample_31_RD.txt", header=FALSE)
sample_32 <- read.delim("sample_32_RD.txt", header=FALSE)
sample_33 <- read.delim("sample_33_RD.txt", header=FALSE)
sample_34 <- read.delim("sample_34_RD.txt", header=FALSE)
sample_35 <- read.delim("sample_35_RD.txt", header=FALSE)
sample_36 <- read.delim("sample_36_RD.txt", header=FALSE)
sample_37 <- read.delim("sample_37_RD.txt", header=FALSE)
sample_38 <- read.delim("sample_38_RD.txt", header=FALSE)
sample_39 <- read.delim("sample_39_RD.txt", header=FALSE)
# wild-type
sample_40 <- read.delim("sample_40_RD.txt", header=FALSE)


all.samples <- ls()
removeChrM <- function(sample){
  sample <- sample[sample$V1!="chrM",]
  return (sample)
}

sample_1 <- removeChrM(sample_1)
sample_2 <- removeChrM(sample_2)
sample_3 <- removeChrM(sample_3)
sample_4 <- removeChrM(sample_4)
sample_5 <- removeChrM(sample_5)
sample_6 <- removeChrM(sample_6)
sample_7 <- removeChrM(sample_7)
sample_8 <- removeChrM(sample_8)
sample_9 <- removeChrM(sample_9)
sample_10 <- removeChrM(sample_10)

sample_11 <- removeChrM(sample_11)
sample_12 <- removeChrM(sample_12)
sample_13 <- removeChrM(sample_13)
sample_14 <- removeChrM(sample_14)
sample_15 <- removeChrM(sample_15)
sample_16 <- removeChrM(sample_16)
sample_17 <- removeChrM(sample_17)
sample_18 <- removeChrM(sample_18)
sample_19 <- removeChrM(sample_19)
sample_20 <- removeChrM(sample_20)

sample_21 <- removeChrM(sample_21)
sample_22 <- removeChrM(sample_22)
sample_23 <- removeChrM(sample_23)
sample_24 <- removeChrM(sample_24)
sample_25 <- removeChrM(sample_25)
sample_26 <- removeChrM(sample_26)
sample_27 <- removeChrM(sample_27)
sample_28 <- removeChrM(sample_28)
sample_29 <- removeChrM(sample_29)
sample_30 <- removeChrM(sample_30)

sample_31 <- removeChrM(sample_31)
sample_32 <- removeChrM(sample_32)
sample_33 <- removeChrM(sample_33)
sample_34 <- removeChrM(sample_34)
sample_35 <- removeChrM(sample_35)
sample_36 <- removeChrM(sample_36)
sample_37 <- removeChrM(sample_37)
sample_38 <- removeChrM(sample_38)
sample_39 <- removeChrM(sample_39)
sample_40 <- removeChrM(sample_40)

end.position.of.each.chr <- as.data.frame(table(sample_40$V1))
end.position.of.each.chr$Var1 <- gsub( "chr", "", as.character(end.position.of.each.chr$Var1))
end.position.of.each.chr$Var1 <- as.numeric(as.roman(end.position.of.each.chr$Var1))


genome_plots <- function (s1,s2,s3,s4,s5,ancestral,name1,name2,name3,name4,name5,ancestral_name){
  
  par(mai = c(1.5,1.0,1.0,2.5)) 
  plot(seq(1:dim(s1)[1]),(s1$V3/mean(s1$V3)),type="n",pch=".",ylim=c(0,3),
       main="Whole genome plot", 
       xlab="Position across entire genome",ylab="Read depth relative to mean",
       xaxt="n",xaxs="i", yaxs="i",xlim=c(0,12000000))
  eaxis(1)
  
  rect(0,0.01,230218,2.99,col="lightcyan",border="transparent")
  rect(1043402,0.01,1360022,2.99,col="lightcyan",border="transparent")
  rect(2891908,0.01,3468767,2.99,col="lightcyan",border="transparent")
  rect(3738928,0.01,4829858,2.99,col="lightcyan",border="transparent")
  rect(5392501,0.01,5832388,2.99,col="lightcyan",border="transparent")
  rect(3738928,0.01,4829858,2.99,col="lightcyan",border="transparent")
  rect(6578139,0.01,7248330,2.99,col="lightcyan",border="transparent")
  rect(8326507,0.01,9250938,2.99,col="lightcyan",border="transparent")
  rect(10035271,0.01,11126562,2.99,col="lightcyan",border="transparent")
  
  
  lines(seq(1:dim(ancestral)[1]),runmed((sample_40$V3/mean(sample_40$V3)),140001), col="black",lwd=1)
  lines(seq(1:dim(s1)[1]),runmed((s1$V3/mean(s1$V3)),90001), col="violetred4",lwd=1)
  lines(seq(1:dim(s2)[1]),runmed((s2$V3/mean(s2$V3)),90001), col="green3",lwd=1)
  lines(seq(1:dim(s3)[1]),runmed((s3$V3/mean(s3$V3)),90001), col="blue",lwd=1)
  lines(seq(1:dim(s4)[1]),runmed((s4$V3/mean(s4$V3)),90001), col="red",lwd=1)
  lines(seq(1:dim(s5)[1]),runmed((s5$V3/mean(s5$V3)),90001), col="darkorchid1",lwd=1)
  
  
  legend(x = 12.7*(10^6), y = 3, xpd = TRUE,  merge = TRUE, xjust = 0,yjust = 1, lty=c(1,1), lwd=c(2.5,2.5),
         c(ancestral_name,name1,name2,name3,name4,name5),
         col=c("black","violetred4","green3","blue","red","darkorchid1"),
         bty = "n",inset=c(1,1), ncol = 1)
  
}



jpeg("/home/ggg256/RD/N0.12345.jpg")
genome_plots (sample_1,sample_2,sample_3,sample_4,sample_5,sample_40,
                 "sample1","sample2","sample3","sample4","sample5","ancestral")
dev.off()

jpeg("/home/ggg256/RD/N0.67890.jpg")
genome_plots (sample_6,sample_7,sample_8,sample_9,sample_10,sample_40,
              "sample6","sample7","sample8","sample9","sample10","ancestral")
dev.off()

#__________________________________
jpeg("/home/ggg256/RD/N1.12345.jpg")
genome_plots (sample_11,sample_12,sample_13,sample_14,sample_15,sample_40,
              "sample11","sample12","sample13","sample14","sample15","ancestral")
dev.off()

jpeg("/home/ggg256/RD/N1.67890.jpg")
genome_plots (sample_16,sample_17,sample_18,sample_19,sample_20,sample_40,
              "sample16","sample17","sample18","sample19","sample20","ancestral")
dev.off()
#__________________________________
jpeg("/home/ggg256/RD/N2.12345.jpg")
genome_plots (sample_21,sample_22,sample_23,sample_24,sample_25,sample_40,
              "sample21","sample22","sample23","sample24","sample25","ancestral")
dev.off()

jpeg("/home/ggg256/RD/N2.67890.jpg")
genome_plots (sample_26,sample_27,sample_28,sample_29,sample_30,sample_40,
              "sample26","sample27","sample28","sample29","sample30","ancestral")
dev.off()
#__________________________________
jpeg("/home/ggg256/RD/N3.12345.jpg")
genome_plots (sample_31,sample_32,sample_33,sample_34,sample_35,sample_40,
              "sample31","sample32","sample33","sample34","sample35","ancestral")
dev.off()

jpeg("/home/ggg256/RD/N3.67890.jpg")
genome_plots (sample_36,sample_37,sample_38,sample_39,sample_40,sample_40,
              "sample36","sample37","sample38","sample39","sample40","ancestral")
dev.off()



