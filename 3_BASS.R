# Sobol based on the Bayesian Adaptive Spline Surface (BASS) method
# Note: the Sobol analysis in this script directly uses the function in "BASS" package instead of sensobol package

# This script also takes a time when dealing with high-dimensional models (but not extremely long).
#       Again change the dimension vector D for code replication check.
# Remove all existing environment and plots
rm(list = ls())
graphics.off()

source("0_library.R")

# Set a random seed
set.seed(1)

# Define the model in each dimension and apply BASS method
for (k in 1:6){
  # Model dimension
  d=D[k]
  S_total <- rep(NA,d)
  
  # Folder for d dimension test scenario
  folder <- paste(folderpath,d,"D/BASS",sep="")
  if (Testmodel_ind==2){
    folder <- paste(folderpath,"Hymod/BASS",sep="") 
  }
  if (Testmodel_ind==3){
    folder <- paste(folderpath,"SacSma/BASS",sep="") 
  }
  if (!dir.exists(folder)){
    dir.create(file.path(folder), showWarnings = FALSE)
  }
    
  # Use 20,000 LHS training data points to test emulator quality
  x_test <- randomLHS(20000,d)
  if (Testmodel_ind >= 2){
    x_test <- Mapping(x_test,Range)      
  }
  sample_size <- 10*d
  # Similar to Kriging method, we begin with 10 times model dimension samples
  
  while (1>0) {
    X <- randomLHS(sample_size, d)
    if (Testmodel_ind >= 2){
      X <- Mapping(X,Range)      
    }
    Y <- apply(X, 1, Testmodel)
    
    # Use a MCMC size of 500,000, burn-in period of 100,000, record the output every 1,000 steps
    mcmc_size <- 500000
    # Record the time of BASS emulation
    start.time <- Sys.time()
    mod <- bass(X, Y, nmcmc = mcmc_size, nburn = 100000, thin = 1000,verbose = FALSE) # fit BASS model
    end.time <- Sys.time()
    T_BASS <- difftime(end.time,start.time, units = "secs")
    
    y <- predict(mod,x_test)
    std <- sqrt(apply(y, 2, var))
    mean <- colMeans(y)
    print(paste("sample size = ",sample_size,sep=""))
    print(paste("max std = ",max(std), "thres value = ", ((max(mean)-min(mean))/20)^2,sep=""))
    
    # If still need to take more samples, add the sample size by d
    if (max(std) > ((max(mean)-min(mean))/20)){
      sample_size <- sample_size + d
    }
    # otherwise perform sensitivity analysis and record the time
    else{
      start.time <- Sys.time()
      S_BASS <- sobol(mod, verbose = FALSE)
      end.time <- Sys.time()
      T_BASSSobol <- difftime(end.time,start.time,units = "secs")
      
      Important_indices <- as.numeric(names(S_BASS$T[1, ]))

      for (j in 1:length(Important_indices)){
        S_total[Important_indices[j]] <- quantile(S_BASS$T[ ,j],probs = 0.975) - quantile(S_BASS$T[ ,j],probs = 0.025)
      }
      
      # Because we are not using sensobol package here, we need an additional convergence check for the Sobol
      #     analysis. Usually if the emulator's quality is good enough, the sensitivity analysis results should
      #     be stable. We add a convergence check here to make sure the results are reliable. If the results are 
      #     not converged, add the BASS sample size by d.
      
      # If converged, save the emulation time, sensitivity analysis time, sensitivity indices and the sample size
      Sens <- S_BASS$T[seq(4,400,by=4),c(1:length(Important_indices))]
      Rank <- t(apply(Sens, 1, rank))
      for (boot_ind1 in 1:99){
        T <- boot_ind1
        for (boot_ind2 in (T+1):100){
          Rho <- rep(NA,length(Important_indices))
          Weights <- rep(NA,length(Important_indices))
          for (para_ind in 1:length(Important_indices)){
            Weights[para_ind] <- (max(Sens[boot_ind1,para_ind],max(Sens[boot_ind2,para_ind])))^2
          }
          Weights_sum <- sum(Weights)
          for (para_ind in 1:length(Important_indices)){
            Rho[para_ind] <- abs(Rank[boot_ind1,para_ind]-Rank[boot_ind2,para_ind])*
              Weights[para_ind]/Weights_sum
          }
          if (boot_ind2 == 2){
            Rho_all <- Rho
          } else{
            Rho_all <- append(Rho_all,Rho)
          }
        }
      }
      Rho_all <- matrix(Rho_all, nrow = d)
      Rho_all <- apply(Rho_all, 2, sum)

      if (quantile(Rho_all,probs = 0.95, na.rm = TRUE) < 1){
        save(T_BASS, file = paste(folder, "/T_BASS", sep=""))
        save(T_BASSSobol,file = paste(folder,"/T_BASSSobol",sep=""))
        save(S_BASS,file = paste(folder,"/S_BASS",sep=""))
        save(sample_size,file = paste(folder,"/BASS_size",sep=""))
        break
      } else{
        sample_size <- sample_size + d
      }
    }
  }
}
