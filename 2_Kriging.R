# Sobol based on the Kriging emulator
# Note: this script can take an extremely long running time. This is because high-dimensional test problems take
#       a long time to build the Kriging emulator, and we need to search for the required sample size 
#       (which means we need to repeat the emulation process). It would be more efficient to separate high
#       dimension tests by modifying the vector "D".

# For code check, keep only the first two terms in "D" (avoid high-dimensional runs) 
#       so that the total time is acceptable. (see the block with dashed lines below)

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
setwd("/Users/f007f8t/Documents/GitHub/Sensitivity_paper_revision")

# Load the required packages
source("0_library.R")

# Set a random seed
set.seed(3)

# Define the test model in each dimension, build the Kriging emulator and perform the Sobol analysis
for (k in 1:6){
  # model dimension
  d <- D[k]
  
  # Begin with 10*dimension as a base sample size
  Kriging_size <- 10*d
  
  # The training samples for emulator quality control
  # Select 20,000 samples by Latin Hypercube Sampling (LHS)
  x_test <- randomLHS(20000,d)
  
  # Folder for d dimension test scenario
  folder <- paste(folderpath,d,"D/Kriging",sep="")
  if (Testmodel_ind==2){
    folder <- paste(folderpath,"Hymod/Kriging",sep="") 
  }
  if (Testmodel_ind==3){
    folder <- paste(folderpath,"SacSma/Kriging",sep="") 
  }
  if (!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
  
  save(x_test,file = paste(folder,"/x_test",sep=""))

  # A while loop includes the stopping criterion
  while (1>0){
    # Begin with the base sample size (also sampled by LHS)
    X_GP <- randomLHS(Kriging_size,d)
    # Evaluate the model to get their true outputs
    if (Testmodel_ind >= 2){
      Y_GP <- apply(Mapping(X_GP,Range),1,Testmodel) 
    } else {
      Y_GP <- apply(X_GP,1,Testmodel)
    }
    
    # Record the Kriging model evaluation time
    start.time <- Sys.time()
    GPmodel <- GP_fit(X_GP,Y_GP)
    end.time <- Sys.time()
    
    # Emulator convergence check
    # First get the standard error of all the training samples
    a <- predict(GPmodel,x_test)
    std <- sqrt(a$MSE)
    print(paste("k =",k,"max(U) =",max(std),"range = ",(max(a$Y_hat)-min(a$Y_hat))/20,sep=" "))
    
    # If all the standard errors are less than or equal to 1, then convergence is reached, else add sample size
    if (max(std)<=(max(a$Y_hat)-min(a$Y_hat))/20){
      save(a,file = paste(folder,"/a",sep=""))
      break
    } else {
      Kriging_size <- Kriging_size + d
    }
  }
  
  # Time for building the emulator
  T_Kriging<-difftime(end.time,start.time,units = "secs")
  
  # Save this time and the required sample size
  save(T_Kriging,file = paste(folder,"/T_Kriging",sep=""))
  save(Kriging_size,file = paste(folder,"/Kriging_size",sep=""))
  
  # Then perform the sensitivity analysis
  # The sample size in sensitivity analysis is directly estimated by the required sample size of the original
  #     model for convenience.
  if (Testmodel_ind == 1){
    load(paste(folderpath,"Sobol_convergencesize",sep = ""))
  }
  if (Testmodel_ind == 2){
    load(paste(folderpath,"Hymod/Sobol_convergencesize",sep = ""))
  }
  if (Testmodel_ind == 3){
    load(paste(folderpath,"SacSma/Sobol_convergencesize",sep = ""))
  }
  N <- Sobol_convergesize[k]/(d+2+d*(d-1)/2)
  
  # Record the time for sensitivity analysis
  start.time<-Sys.time()
  mat <- sobol_matrices(N = N, params = as.character(c(1:d)), order = "second")
  Y_S <- Kriging(mat)
  
  # The sensitivity indices
  S_Kriging <- sobol_indices(Y=Y_S,N=N,params = as.character(c(1:d)),
                             boot=TRUE,R=100,order="second")
  end.time <- Sys.time()
  
  # Save the results
  T_KrigingSobol <- difftime(end.time,start.time,units = "secs")
  save(T_KrigingSobol,file = paste(folder,"/T_KrigingSobol",sep=""))
  save(S_Kriging,file=paste(folder,"/S_Kriging",sep=""))
}

