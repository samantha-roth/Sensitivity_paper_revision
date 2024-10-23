# This script makes the grid plots (Figure 7 and 8 in the paper)
# Note: this script requires the full data. 

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

source("0_library.R")

# Load the required package for plotting
library(plot.matrix)
library(RColorBrewer)

# Tested dimension, method names, and evaluation time
tested_D_num <- c(2,5,10,15,20,30)
tested_D <- c("2D","5D","10D","15D","20D","30D")
tested_M <- c("Kriging","AKMCS","BASS")
tested_eval_time <- c(0.001,0.01,0.1,1,10,60,600,3600,3600*10)
# Label of evaluation time
eval_time_lab <- c("1ms","10ms","0.1s","1s","10s","1min","10min","1h","10h")

# Total time of each method & each scenario
Time_Sobol <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_Kriging <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_AKMCS <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_BASS <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))

# Load all related results for each test scenario
for (i in 1:length(tested_D)) {
  folder <- paste("./Ranking_Data/",tested_D[i],sep="")
  
  # General:
  load(paste(folder,"/avg_eval_time",sep=""))
  
  # Sobol:
  load(paste(folder,"/Sobol/T_Sobol",sep=""))
  load(paste(folder,"/Sobol/S_Sobol",sep=""))
  Sobol_convergesize <- S$C
  
  # Kriging:
  load(paste(folder,"/Kriging/T_Kriging",sep=""))
  load(paste(folder,"/Kriging/T_KrigingSobol",sep=""))
  load(paste(folder,"/Kriging/Kriging_size",sep=""))
  
  # AKMCS:
  load(paste(folder,"/AKMCS/T_AKMCS",sep=""))
  load(paste(folder,"/AKMCS/x",sep=""))
  AKMCS_size <- dim(x)[1]
  load(paste(folder,"/AKMCS/T_AKMCSSobol",sep=""))
  
  # BASS:
  load(paste(folder,"/BASS/T_BASS",sep=""))
  load(paste(folder,"/BASS/T_BASSSobol",sep=""))
  load(paste(folder,"/BASS/BASS_size",sep=""))
  BASS_size <- sample_size
  
  # Calculation of the computational time in each scenario
  # Sensitivity analysis time + model evaluation adjusted time + emulation time
  for (j in 1:length(tested_eval_time)) {
    Time_Sobol[i,j] <- T_Sobol + (tested_eval_time[j]-avg_eval_time)*Sobol_convergesize
    Time_Kriging[i,j] <- T_Kriging + (tested_eval_time[j]-avg_eval_time)*Kriging_size+ T_KrigingSobol
    Time_AKMCS[i,j] <- T_AKMCS + (tested_eval_time[j]-avg_eval_time)*AKMCS_size + T_AKMCSSobol
    Time_BASS[i,j] <- T_BASS + (tested_eval_time[j]-avg_eval_time)*BASS_size + T_BASSSobol
  }
}

# An adjustment for standard Sobol's time: for fast test models, it is likely the time becomes negative due
#     to random noises in the time recording.
Time_Sobol[Time_Sobol < 0] <- Time_Sobol[1,1]

PercentTime_Kriging <- Time_Sobol/Time_Kriging
PercentTime_AKMCS <- Time_Sobol/Time_AKMCS
PercentTime_BASS <- Time_Sobol/Time_BASS

#-------------------------------------------------------------
# Which method to choose? (Figure 7)
# Identify the fastest method for each scenario
Mat <- Time_Sobol
textMat <- matrix(NA,nrow=nrow(Mat),ncol=ncol(Mat))
for (i in 1:length(tested_D)) {
  for (j in 1:length(tested_eval_time)){
    Mat[i,j] <- min(Time_Sobol[i,j],Time_Kriging[i,j],
                    Time_AKMCS[i,j],Time_BASS[i,j],na.rm=TRUE)
    if (Mat[i,j]==Time_Sobol[i,j]){
      textMat[i,j] <- "Sobol"
    } 
    if (Mat[i,j]==Time_Kriging[i,j]){
      textMat[i,j] <- "Kriging"
    } 
    if (Mat[i,j]==Time_AKMCS[i,j]){
      textMat[i,j] <- "AKMCS"
    } 
    if (Mat[i,j]==Time_BASS[i,j]){
      textMat[i,j] <- "BASS"
    }
  }
}
# Row names and column names of the plot
rownames(textMat) <- tested_D_num
colnames(textMat) <- eval_time_lab
# Color palette
cols<-brewer.pal(n = 4, name = "Set3")


pdf(file = "./New_Figures/Figure_7.pdf",width = 12,height = 7)
par(mar=c(5,5,2.6,6))
plot(textMat[nrow(Mat):1, ],breaks = c("Sobol","Kriging","BASS","AKMCS"),
     xlab="Time of single run",ylab="Number of input parameters",col=cols,
     cex.lab=1.5,cex.axis=1,main="")
dev.off()
#----------------------------------------------------

# Absolute computational time of this method (Figure 8)
Mat<-Time_Sobol
textMat <- matrix(NA,nrow=nrow(Mat),ncol=ncol(Mat))
for (i in 1:length(tested_D)) {
  for (j in 1:length(tested_eval_time)){
    Mat[i,j] <- min(Time_Sobol[i,j],Time_Kriging[i,j],
                    Time_AKMCS[i,j],Time_BASS[i,j],na.rm=TRUE)
    if (Mat[i,j]==Time_Sobol[i,j]){
      textMat[i,j] <- "S"
    }
    if (Mat[i,j]==Time_Kriging[i,j]){
      textMat[i,j] <- "K"
    }
    if (Mat[i,j]==Time_AKMCS[i,j]){
      textMat[i,j] <- "A"
    }
    if (Mat[i,j]==Time_BASS[i,j]){
      textMat[i,j] <- "B"
    }
    Mat[i,j] <- Mat[i,j]/Time_Sobol[1,1]
  }
}
Mat[which(Mat<1)] <- 1
Mat <- round(log10(Mat))
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
Cols<-brewer.pal(n = 7, name = "Reds")
pdf(file = "./New_Figures/Figure_8.pdf",width = 12,height = 7)
par(mar=c(5,5,2.6,6))
time<-plot(Mat[nrow(Mat):1, ],breaks=c(0:7),digits=1,fmt.cell='%.0f',
           xlab="Time of single run",ylab="Number of input parameters",col=Cols,
           text.cell=list(pos=3, cex=1),max.col=100,
           cex.lab=1.5,cex.axis=1,main="",fmt.key = "%.0f")
for (i in 1:nrow(Mat)) {
  for (j in 1:ncol(Mat)) {
    args<-time$cell.text[[nrow(Mat)+1-i,j]]
    args$labels <- textMat[i,j]
    args$cex    <- 1.5
    args$pos    <- 1
    do.call(text, args)
  }
}
mtext("Orders of",side = 3, at = 15)
mtext("magnitude",side = 3, at = 15, line = -1)
dev.off()

#---------------------------------below is the script for table 3's data
rm(list=ls())
#Load Hymod and Sacsma's time
load("Ranking_Data/Hymod/avg_eval_time")

# Sobol:
load("Ranking_Data/Hymod/Sobol/T_Sobol")
load("Ranking_Data/Hymod/Sobol_convergencesize")
Time_Sobol <- T_Sobol + avg_eval_time*Sobol_convergesize

# Kriging:
load("Ranking_Data/Hymod/Kriging/T_Kriging")
load("Ranking_Data/Hymod/Kriging/T_KrigingSobol")
load("Ranking_Data/Hymod/Kriging/Kriging_size")
Time_Kriging <- T_Kriging + avg_eval_time*Kriging_size+ T_KrigingSobol

# AKMCS:
load(paste("Ranking_Data/Hymod/AKMCS/T_AKMCS",sep=""))
load(paste("Ranking_Data/Hymod/AKMCS/x",sep=""))
AKMCS_size <- dim(x)[1]
load(paste("Ranking_Data/Hymod/AKMCS/T_AKMCSSobol",sep=""))
Time_AKMCS <- T_AKMCS + T_AKMCSSobol

# BASS:
load(paste("Ranking_Data/Hymod/BASS/T_BASS",sep=""))
load(paste("Ranking_Data/Hymod/BASS/T_BASSSobol",sep=""))
load(paste("Ranking_Data/Hymod/BASS/BASS_size",sep=""))
Time_BASS <- T_BASS + avg_eval_time*sample_size + T_BASSSobol

#-----------sacsma

load("Ranking_Data/SacSma/avg_eval_time")

# Sobol:
load("Ranking_Data/SacSma/Sobol/T_Sobol")
load("Ranking_Data/SacSma/Sobol_convergencesize")
Time_Sobol <- T_Sobol + avg_eval_time*Sobol_convergesize

# Kriging:
load("Ranking_Data/SacSma/Kriging/T_Kriging")
load("Ranking_Data/SacSma/Kriging/T_KrigingSobol")
load("Ranking_Data/SacSma/Kriging/Kriging_size")
Time_Kriging <- T_Kriging + avg_eval_time*Kriging_size + T_KrigingSobol

# AKMCS:
load(paste("Ranking_Data/SacSma/AKMCS/T_AKMCS",sep=""))
load(paste("Ranking_Data/SacSma/AKMCS/x",sep=""))
AKMCS_size <- dim(x)[1]
load(paste("Ranking_Data/SacSma/AKMCS/T_AKMCSSobol",sep=""))
Time_AKMCS <- T_AKMCS + T_AKMCSSobol

# BASS:
load(paste("Ranking_Data/SacSma/BASS/T_BASS",sep=""))
load(paste("Ranking_Data/SacSma/BASS/T_BASSSobol",sep=""))
load(paste("Ranking_Data/SacSma/BASS/BASS_size",sep=""))
Time_BASS <- T_BASS + avg_eval_time*sample_size + T_BASSSobol


load(paste("Ranking_Data/SacSma/Sobol/S_Sobol",sep=""))
load(paste("Ranking_Data/SacSma/BASS/S_BASS",sep=""))
S_B <- S_BASS$S
