# This script prepares the required data of Figure 3 in the paper

# Remove all existing environment and plots
rm(list = ls())
graphics.off()
source("0_library.R")

# Plot the results in the example of 5D problem
d <- 5
folder<-paste("./Ranking_Data/",d,"D",sep="")
if (!dir.exists(paste(folder,"/Traceplot",sep=""))){
  dir.create(paste(folder,"/Traceplot",sep=""), recursive = TRUE)
}
# Test model in 5D
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



# Sizes used for Sobol trace plot
# Standard Sobol's results do not vary by random seeds, other methods' results do. 
Size_S <- c(seq(40,1000,by=40),seq(1000,5000,by=1000))
# Upper and lower bound of 95% CI of the total-order indices
T_S_high <-rep(NA, length(Size_S))
T_S_low <- rep(NA, length(Size_S))
T_S <- rep(NA, length(Size_S))

#Evaluate the results of the standard Sobol method
for (i in 1:length(Size_S)){
  N <- floor(Size_S[i]/(d+2+d*(d-1)/2))
  mat <- sobol_matrices(N = N, params = as.character(c(1:d)), order = "second")
  X <- split(t(mat), rep(1:dim(mat)[1], each = d))
  Y <- sapply(X, Testmodel)
  S <- sobol_indices(Y=Y,N=N,params = as.character(c(1:d)),
                     boot=TRUE,R=100,order="second")
  T_S_high[i] <- S$results$high.ci[6]
  T_S_low[i] <- S$results$low.ci[6]
  T_S[i] <- S$results$original[6]
}
save(T_S,file = paste(folder,"/Traceplot/T_S",sep=""))
save(T_S_high,file = paste(folder,"/Traceplot/T_S_high",sep=""))
save(T_S_low,file = paste(folder,"/Traceplot/T_S_low",sep=""))

# Sizes used for Kriging traceplot
# 5 different seeds for the emulation methods
Size_K <- c(5:10,seq(12,100,by=2))

# Total-order indices
T_K <- matrix(NA,nrow=5,ncol=length(Size_K))

# Sizes used for AKMCS traceplot
Size_A <- c(3:100)
T_A <- matrix(NA,nrow=5,ncol=length(Size_A))

# Sizes used for BASS traceplot
Size_B <- c(5:10,seq(12,100,by=2))
T_B <- matrix(NA,nrow=5,ncol=length(Size_B))

# Test for 5 seeds for the emulation methods
for (seed in 1:5){
  set.seed(seed*4)
  print(paste("seed=",seed," now Kriging",sep=""))
  # Kriging
  for (i in 1:length(Size_K)){
    # Take the corresponding sample size, fit a Kriging model and perform Sobol analysis
    X_GP <- randomLHS(Size_K[i],d)
    Y_GP <- apply(X_GP,1,Testmodel)
    GPmodel <- GP_fit(X_GP,Y_GP)
    mat <- sobol_matrices(N = 500, params = as.character(c(1:d)), order = "second")
    Y_K <- Kriging(mat)
    S <- sobol_indices(Y = Y_K,N = 500, params = as.character(c(1:d)),
                               boot=TRUE, R=100, order="second")
    T_K[seed,i] <- S$results$original[6]
  }
  # Save the results of the best estimates
  save(T_K,file = paste(folder,"/Traceplot/T_K",sep=""))

  set.seed(seed)
  print(paste("seed=",seed," now AKMCS",sep=""))
  # AKMCS records the results after each step
  # 20,000 training points
  candidate_size <- 20000
  X <- randomLHS(candidate_size,d)
  # Begin with 15 initial samples
  n_init <- 15
  indx <- sample(candidate_size,n_init)
  x <- X[indx, ]
  x_rest <- X[-indx, ]
  # Evaluate their outputs, fit a Kriging model
  y <- apply(x,1,Testmodel)
  GPmodel <- GP_fit(x, y)
  a<-predict(GPmodel,x_rest)
  # Learning function
  U<-sqrt(a$MSE)
  # Take the next sample adaptively, repeat these steps
  for (i in 1:length(Size_A)){
    k<-which(U==max(U))
   if (length(k)>1){
      k <- sample(k,1)
    }
    x_add<-x_rest[k, ]
    x_rest<-x_rest[-k, ]
    y_add<-Testmodel(x_add)
    y <- append(y,y_add)
    x <- rbind(x,x_add)
    GPmodel <- GP_fit(x, y)
    a<-predict(GPmodel,x_rest)
    U<-sqrt(a$MSE)
    
    # Save the sensitivity analysis results
    mat <- sobol_matrices(N = 500, params = as.character(c(1:d)), order = "second")
    Y_A <- Kriging(mat)
    S <- sobol_indices(Y = Y_A,N = 500, params = as.character(c(1:d)),
                               boot=TRUE, R=100, order="second")
    T_A[seed,i] <- S$results$original[6]
  }
  save(T_A,file = paste(folder,"/Traceplot/T_A",sep=""))

  set.seed(seed)
  print(paste("seed=",seed," now BASS",sep=""))
  # BASS method simply takes different initial sample sizes
  for (i in 1:length(Size_B)){
    X <- randomLHS(Size_B[i], d)
    Y <- apply(X, 1, Testmodel)
    mcmc_size <- 500000
    mod <- bass(X, Y, nmcmc = mcmc_size, nburn = 100000, thin = 1000,verbose = FALSE) # fit BASS model
    S_BASS <- sobol(mod, verbose = FALSE)
    # Take the mean of the ensemble as the best estimate
    T_B[seed,i] <- mean(S_BASS$T[ ,1])
  }
  save(T_B,file = paste(folder,"/Traceplot/T_B",sep=""))
}
