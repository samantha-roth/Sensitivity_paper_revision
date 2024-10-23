# This script plots the traceplot (Figure 6 of the paper)

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

source("0_library.R")

# Load the required package for plotting
library(RColorBrewer)

# Sizes used for standard Sobol traceplot
Size_S <- c(seq(40,1000,by=40),seq(1000,5000,by=1000))
# Sizes used for Kriging traceplot
Size_K <- c(5:10,seq(12,120,by=2))
# Sizes used for AKMCS traceplot
Size_A <- c(3:100)
# Sizes used for BASS traceplot
Size_B <- c(5:10,seq(12,100,by=2))

# Load the saved sensitivity indices
load("./Ranking_Data/Hymod/Traceplot_rank/T_S")
load("./Ranking_Data/Hymod/Traceplot_rank/T_K")
load("./Ranking_Data/Hymod/Traceplot_rank/T_B")
load("./Ranking_Data/Hymod/Traceplot_rank/T_A")

# Load the convergence size of the Hymod test model
load("./Ranking_Data/Hymod/Sobol_convergencesize")
C_S <- Sobol_convergesize
load("./Ranking_Data/Hymod/Kriging/Kriging_size")
C_K <- Kriging_size
load("./Ranking_Data/Hymod/AKMCS/x")
C_A <- dim(x)[1]
load("./Ranking_Data/Hymod/BASS/BASS_size")
C_B <- sample_size

# Create a folder to save figures
folder <- "./New_Figures"
if (!dir.exists(folder)){
  dir.create(folder, recursive = TRUE)
}

# 4 panels of trace plots + 1 line showing convergence locations
pdf(file = paste("./New_Figures/Figure_6.pdf",sep=""),width = 18,height = 12)
layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE))
par(mar=c(5,6,6,2.6))
Yrange <- c(min(T_S,T_A,T_B,T_K)-0.05,max(T_S,T_A,T_B,T_K)+0.05)
plot(Size_S,T_S[1, ],type="l",col="seagreen",xlab="Sample size",ylab="Sensitivity",main="Sobol",
     ylim=Yrange,cex.axis=2,cex.lab=2.5,cex.main=2)
for (i in 2:5){
  lines(Size_S,T_S[i, ],lty=i,col="seagreen")
}
legend("topright",lty = c(1:5), col = rep("seagreen",5), 
       legend = c("Parameter 1","Parameter 2","Parameter 3","Parameter 4","Parameter 5"), bty = "n", cex = 1.5)
mtext("a",side = 3, line = 1, at = 0, cex = 2)
arrows(C_S,1.1,C_S,0.9,length = 0.1, col = "seagreen")

par(mar=c(5,6,6,2.6))
plot(Size_K,T_K[1, ],type="l",col="purple",ylim = Yrange,main = "Kriging",
     xlab="Sample size",ylab="Sensitivity",cex.axis=2,cex.lab=2.5,cex.main=2)
for (i in 2:5){
  lines(Size_K,T_K[i, ],lty=i,col="purple")
}
legend("topright",lty = c(1:5), col = rep("purple",5), 
       legend = c("Parameter 1","Parameter 2","Parameter 3","Parameter 4","Parameter 5"), bty = "n", cex = 1.5)
mtext("b",side = 3, line = 1, at = min(Size_K), cex = 2)
arrows(C_K,0.35,C_K,0.55,length = 0.1,col="purple")

par(mar=c(5,6,2.6,2.6))
plot(Size_B,T_B[1, ],type="l",col="blue",ylim = Yrange,main = "BASS",
     xlab="Sample size",ylab="Sensitivity",cex.axis=2,cex.lab=2.5,cex.main=2)
for (i in 2:5){
  lines(Size_B,T_B[i, ],lty=i,col="blue")
}
legend("topright",lty = c(1:5), col = rep("blue",5), 
       legend = c("Parameter 1","Parameter 2","Parameter 3","Parameter 4","Parameter 5"), bty = "n", cex = 1.5)
mtext("c",side = 3, line = 1, at = min(Size_B), cex = 2)
arrows(C_B,1,C_B,0.8,length = 0.1,col="blue")

par(mar=c(5,6,2.6,2.6))
plot(Size_A,T_A[1, ],type="l",col="red",ylim = Yrange,main="AKMCS",
     xlab="Sample size",ylab="Sensitivity",cex.axis=2,cex.lab=2.5,cex.main=2)
for (i in 2:5){
  lines(Size_A,T_A[i, ],lty=i,col="red")
}
legend("topright",lty = c(1:5), col = rep("red",5), 
       legend = c("Parameter 1","Parameter 2","Parameter 3","Parameter 4","Parameter 5"), bty = "n", cex = 1.5)
mtext("d",side = 3, line = 1, at = min(Size_A), cex = 2)
arrows(C_A,1,C_A,0.8,length = 0.1,col="red")

plot(0,0,type = "n", xaxt = "n", yaxt = "n", bty="n", xlab = "", ylab="",
     xlim=c(10, 5000), ylim=c(0, 0.7),log="x")
axis(1, at = c(10,100,1000,10000),labels = c(10,100,1000,10000), pos = 0.5,
     cex.axis = 2.3, cex.lab = 2.5)
points(C_S, 0.5, col = "seagreen", pch = 20, cex = 2)
points(C_K, 0.5, col = "purple", pch = 20, cex = 2)
points(C_B, 0.5, col = "blue", pch = 20, cex = 2)
points(C_A, 0.5, col = "red", pch = 20, cex = 2)
text(C_S, 0.6, labels = "Sobol", col = "seagreen", cex = 2.5)
text(C_K, 0.6, labels = "Kriging", col = "purple", cex = 2.5)
text(C_B, 0.4, labels = "BASS", col = "blue", cex = 2.5)
text(C_A, 0.6, labels = "AKMCS", col = "red", cex = 2.5)
text(200, 0.2, labels = "Required sample size for convergence", col = "black", cex = 2.5)
dev.off()
