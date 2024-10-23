# This script records all the global variables, functions and libraries

#The working folder with scripts and data
setwd("/storage/group/pches/default/users/hxy46/Sensitivity")

# Load the required functions and packages
source("sobol_indices_boot.R")
library(GPfit)
library(lhs)
library(BASS)

# The dimensions we consider
D <- c(2,5,10,15,20,30)

# The list of sample sizes we consider
tot_size <- c(seq(500,1000,by=100),seq(1250,5000,by=250),seq(6000,50000,by=1000))

# The root folder that saves all the data during code running
folderpath <- "./Ranking_Data/"

# Choose the test model, the number can only be 1 (G function), 2 (Hymod) or 3 (Sacsma)
Testmodel_ind <- 1 

if (Testmodel_ind == 1){ # G function
  Testmodel<-function (X) {
    
    a = rep(NA,d)
    for (i in 1:d){
      a[i] <- i
    }
    
    up <- abs(4*X-2) + a
    down <- 1 + a
    prod <- prod(up/down)
    
    return(prod)
  }
}

if (Testmodel_ind == 2){ # Hymod
  source("Hymod.R")
  load("arnosubbiano.rda")
  
  Testmodel <- hymodr
  D <- 5
  
  ## Catchment area [km^2]
  area<-751  
  
  ## Precipitation [mm/day]
  precipit<-arnosubbiano[ ,2]
  
  ## Potential arnosubbiano [mm/day]
  evapo<-arnosubbiano[ ,3]
  
  ## Q observation
  Q_obs <- arnosubbiano[ ,4]
  
  Range <- matrix(data = c(70,400,0,0.95,0,0.9,0,0.1,0.1,1),nrow = 2,ncol = 5,byrow = FALSE)
  
}

if (Testmodel_ind == 3){ # SacSma
  source("sacSma.R")
  load("SacSma_dataframe")
  Testmodel <- sacSma
  D <- 16
  
  ## Precipitation [mm/day]
  precipit<-SacSma_dataframe$precipit
  
  ## Potential arnosubbiano [mm/day]
  evapo<-SacSma_dataframe$evapo
  
  ## Q observation
  Q_obs <- SacSma_dataframe$Q_obs
  
  Range <- matrix(data = c(5,300,5,150,10,700,100,1200,
                           5,500,0.1,0.75,0.001,0.05,0.01,0.6,
                           5,350,1,5,0,0.9,0,0.1,
                           0,0.4,0.3,0.3,0,0,0,0),nrow = 2,ncol = 16,byrow = FALSE)
  
}

# A function that maps [0,1] to any input range
Mapping <- function (X, Range){ 
  # Range is a matrix with 2 rows (min and max)
  if (dim(Range)[2]!=dim(X)[2]){
    return("The input dimensions must match")
  }
  for (i in 1:dim(X)[2]){
    X[ ,i] <- Range[1,i] + X[ ,i]*(Range[2,i]-Range[1,i])
  }
  return(X)
}


#Kriging function, used when calling the emulator in the Sobol analysis
Kriging <- function (X){
  a <- predict(GPmodel,X)
  a$Y_hat
}