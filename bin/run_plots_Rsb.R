#Rscript --vanilla run_plots_Rsb.R [filename] [times difference between pops] [number of pops] {name pop1} ... [name pop n]

args = commandArgs(trailingOnly=TRUE)

#if (!require("data.table")) install.packages("data.table")
library("data.table")
if (length(args) < 4 ){
  system("echo Error: Arguments are: [filename] [times difference between pops] [number of pops] [name pop1] ... [name pop n]")
  #filename <- "../example/Window_1.txt_Results_Tang.txt"
  #npops <- 6
  #times.diff <- log(100)
  #namepops <- c("BRG","BRM","CHCA","CHCU","CHFR","TXL")
  quit()
} else {
  filename <- args[1]
  times.diff <- log(as.numeric(args[2]))
  npops <- as.numeric(args[3])
  namepops <- NULL
  for(i in 1:npops) {
    namepops <- c(namepops,args[3+i])
  }
}

data_stats <- data.frame(fread(filename, header=T))
bp <- data_stats[,1]#/1000000

pdf(sprintf("%s_Results_Rsb.pdf",filename), width=12, height=6)
k <- 1
for(i in 1:(npops-1)) {
  for(j in (i+1):(npops)) {
    col_RsbN <- 1 + npops*2 + (npops*(npops-1)/2)
    lnRsb <- data.frame(bp,data_stats[,col_RsbN+k])
    
    filter_rows <- apply(lnRsb,c(1,2),is.na)
    filter_rows <- apply(filter_rows,1,sum)
    lnRsb <- lnRsb[!filter_rows,]
    
    plot (lnRsb[lnRsb[,2]>=0,1],lnRsb[lnRsb[,2]>=0,2], 
          col="green", main=sprintf("Rsb along chromosome: %s vs %s",namepops[i],namepops[j]), 
          xlab="Chromosome Position bp", ylab="ln(Rsb)",pch=10, type="p",cex=.2, 
          ylim=c(min(lnRsb[,2]),max(lnRsb[,2])), xlim=c(min(lnRsb[,1]),max(lnRsb[,1])))
    lines(lnRsb[lnRsb[,2]<0,1] ,lnRsb[lnRsb[,2]<0,2] , col="red"  ,pch=10, type="p",cex=.2)
    abline(h=0, col=1)
    abline(h=times.diff,col="black", lty=2)
    abline(h=-times.diff,col="black", lty=2)
    
    filter_rows_significance <- which(lnRsb[,2]>=times.diff | lnRsb[,2]<=-times.diff)
    write.table(x=lnRsb[filter_rows_significance,],file=sprintf("%s_Significant_Results_RsbN_%s-%s.txt",filename,namepops[i],namepops[j]), quote = F, sep = "\t",col.names = T,row.names = F,append = F, eol = '\n')

    k<- k + 1
  }
}
dev.off()
