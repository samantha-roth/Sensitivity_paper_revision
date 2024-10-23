library(sensobol)
sobol_indices_boot <- function (matrices = c("A", "B", "AB"), Y, N, params, first = "saltelli", 
          total = "jansen", order = "first", boot = FALSE, R = NULL, 
          parallel = "no", ncpus = 1, conf = 0.95, type = "norm") 
{
  if (boot == FALSE & is.null(R) == FALSE | boot == TRUE & 
      is.null(R) == TRUE) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }
  sensitivity <- parameters <- NULL
  k <- length(params)
  d <- matrix(Y, nrow = N)
  if (boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, 
                      total = total, order = order, boot = FALSE, matrices = matrices)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  }
  else if (boot == TRUE) {
    tmp <- boot::boot(data = d, statistic = sobol_boot, 
                      R = R, N = N, params = params, first = first, total = total, 
                      order = order, matrices = matrices, parallel = parallel, 
                      ncpus = ncpus, boot = TRUE)
    out <- data.table::data.table(bootstats(tmp, conf = conf, 
                                            type = type))
  }
  else {
    stop("boot has to be TRUE or FALSE")
  }
  if (order == "first") {
    parameters <- c(rep(params, times = 2))
    sensitivity <- c(rep(c("Si", "Ti"), each = k))
  }
  else if (order == "second") {
    vector.second <- unlist(lapply(utils::combn(params, 
                                                2, simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    sensitivity <- c(rep(c("Si", "Ti"), each = length(params)), 
                     rep("Sij", times = length(vector.second)))
  }
  else if (order == "third") {
    vector.second <- unlist(lapply(utils::combn(params, 
                                                2, simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, 
                                               simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    sensitivity <- c(rep(c("Si", "Ti"), each = k), rep("Sij", 
                                                       times = length(vector.second)), rep("Sijl", times = length(vector.third)))
  }
  else if (order == "fourth") {
    vector.second <- unlist(lapply(utils::combn(params, 
                                                2, simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, 
                                               simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    vector.fourth <- unlist(lapply(utils::combn(params, 
                                                4, simplify = FALSE), function(x) paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.fourth)
    sensitivity <- c(rep(c("Si", "Ti"), each = k), rep("Sij", 
                                                       times = length(vector.second)), rep("Sijl", times = length(vector.third)), 
                     rep("Sijlm", times = length(vector.fourth)))
  }
  else {
    stop("order has to be first, second, third or fourth")
  }
  ind <- structure(list(), class = "sensobol")
  ind$boot <- tmp
  ind$results <- cbind(out, sensitivity, parameters)
  original <- NULL
  ind$si.sum <- ind$results[sensitivity == "Si", sum(original)]
  ind$first <- first
  ind$total <- total
  ind$C <- length(Y)
  return(ind)
}

sobol_boot <- function(d, i, N, params, matrices, R, first, total, order, boot) {
  
  # Stopping rule to check concordance between estimators and sample matrix
  # -------------------------------------------------------------------
  
  ms <- "Revise the correspondence between the matrices and the estimators"
  
  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    if (!first == "saltelli" & !first == "jansen" |
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen") {
      stop(ms)
    }
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    if (!first == "sobol"| !total == "saltelli") {
      stop(ms)
    }
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    
    if (!first == "azzini" | !total == "azzini" &
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen" & !total == "saltelli") {
      
      if (!total == "azzini" | !first == "saltelli" & !first == "jansen" &
          !first == "azzini" & !first == "sobol") {
        
        stop(ms)
      }
    }
  }
  
  # -------------------------------------
  
  k <- length(params)
  
  if (boot == TRUE) {
    m <- d[i, ]
    
  } else if (boot == FALSE) {
    m <- d
  }
  if (order == "second") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE))
    
  } else if (order == "third") {
    k <- length(params) +
      length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE))
    
  } else if (order == "fourth") {
    k <- length(params) +
      length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE)) +
      length(utils::combn(params, 4, simplify = FALSE))
  }
  
  # Define vectors based on sample design
  # ------------------------------------------------------------------
  
  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_BA <- m[, -c(1, 2)]
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(k + 2)]
    Y_BA <- m[, (k + 3):ncol(m)]
    
  } # A warning might be needed here
  if (isTRUE(all.equal(matrices, c("A", "B", "AB"))) |
      isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }
  
  # Define first-order estimators
  # --------------------------------------------------------------------
  
  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (first == "saltelli" | first == "jansen" | first == "sobol") {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }
  
  # ----------------------------------
  if (first == "sobol") {
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - f0^2)
    
  } else if (first == "saltelli") {
    Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))
    
  } else if (first == "jansen") {
    Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB)^2)
    
  } else if (first == "azzini") {
    VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
    Vi <- (2 * Rfast::colsums((Y_BA - Y_B) * (Y_A - Y_AB)))
    
  } else {
    stop("first should be sobol, saltelli, jansen or azzini")
  }
  
  if (first == "azzini") {
    Si <- Vi[1:length(params)] / VY[1:length(params)]
    
  } else {
    Si <- Vi[1:length(params)] / VY
  }
  
  # Define total-order estimators
  # --------------------------------------------------------------------
  
  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (total == "azzini" | total == "jansen" | total == "sobol" |
      total == "homma" | total == "janon" | total == "glen" | total == "saltelli") {
    
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }
  
  # ----------------------------------
  if (total == "jansen") {
    Ti <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB)^2)) / VY
    
  } else if (total == "sobol") {
    Ti <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY
    
  } else if (total == "homma") {
    Ti <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0^2) / VY
    
  } else if (total == "saltelli") {
    Ti <- 1 - ((1 / N * Rfast::colsums(Y_B * Y_BA - f0^2)) / VY)
    
  } else if (total == "janon") {
    Ti <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) -
                 (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2) /
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB^2) / 2) -
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2)
    
  } else if (total == "glen") {
    Ti <- 1 - (1 / (N - 1) *
                 Rfast::colsums(((Y_A - mean(Y_A)) * (Y_AB - Rfast::colmeans(Y_AB))) /
                                  sqrt(stats::var(Y_A) * Rfast::colVars(Y_AB))))
    
  } else if (total == "azzini") {
    Ti <- Rfast::colsums((Y_B - Y_BA)^2 + (Y_A - Y_AB)^2) /
      Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
    
  } else {
    stop("total should be jansen, sobol, homma saltelli, janon, glen or azzini")
  }
  Ti <- Ti[1:length(params)]
  
  # Define computation of second-order indices
  # ---------------------------------------------------------------------
  
  if (order == "second" | order == "third" | order == "fourth") {
    com2 <- utils::combn(1:length(params), 2, simplify = FALSE)
    tmp2 <- do.call(rbind, com2)
    Vi2 <- Vi[(length(params) + 1):(length(params) + length(com2))] # Second order
    Vi1 <- lapply(com2, function(x) Vi[x])
    final.pairwise <- t(mapply(c, Vi2, Vi1))
    Vij <- unname(apply(final.pairwise, 1, function(x) Reduce("-", x)))
    
    if (first == "azzini") {
      VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
      Sij <- Vij / VY[(length(params) + 1):(length(params) + ncol(utils::combn(1:length(params), 2)))]
      
    } else {
      Sij <- Vij / VY
    }
    
  } else {
    Sij <- NULL
  }
  
  # Define computation of third-order indices
  # ---------------------------------------------------------------------
  
  if (order == "third" | order == "fourth") {
    com3 <- utils::combn(1:length(params), 3, simplify = FALSE)
    tmp3 <- do.call(rbind, com3)
    Vi3 <- Vi[(length(params) + length(com2) + 1):(length(params) + length(com2) + length(com3))] # Third order
    Vi1 <- lapply(com3, function(x) Vi[x])
    
    # Pairs
    tmp <- do.call(rbind, com2)
    Vij.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], sep = ""))
    names(Vij) <- Vij.vec
    Vi2.1 <- unname(Vij[paste(tmp3[, 1], tmp3[, 2], sep  = "")])
    Vi2.2 <- unname(Vij[paste(tmp3[, 1], tmp3[, 3], sep  = "")])
    Vi2.3 <- unname(Vij[paste(tmp3[, 2], tmp3[, 3], sep  = "")])
    
    mat3 <- cbind(Vi3, Vi2.1, Vi2.2, Vi2.3, do.call(rbind, Vi1))
    Vijk <- unname(apply(mat3, 1, function(x) Reduce("-", x)))
    
    if (first == "azzini") {
      Sijl <- Vijk / VY[(length(params) + length(com2) + 1):(length(params) + length(com2) + length(com3))]
      
    } else {
      Sijl <- Vijk / VY
    }
    
  } else {
    Sijl <- NULL
  }
  
  # Define computation of fourth-order indices
  # ---------------------------------------------------------------------
  
  if(order == "fourth") {
    com4 <- utils::combn(1:length(params), 4, simplify = FALSE)
    tmp4 <- do.call(rbind, com4)
    Vi4 <- Vi[(length(params) + length(com2) + length(com3) + 1):
                (length(params) + length(com2) + length(com3) + length(com4))] # Fourth order
    Vi1 <- lapply(com4, function(x) Vi[x])
    
    # triplets
    tmp <- do.call(rbind, com3)
    Vijk.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], tmp[, 3], sep = ""))
    names(Vijk) <- Vijk.vec
    Vi3.1 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 2], tmp4[, 3], sep  = "")])
    Vi3.2 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 2], tmp4[, 4], sep  = "")])
    Vi3.3 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 3], tmp4[, 4], sep  = "")])
    Vi3.4 <- unname(Vijk[paste(tmp4[, 2], tmp4[, 3], tmp4[, 4], sep  = "")])
    
    # Pairs
    Vi2.1 <- unname(Vij[paste(tmp4[, 1], tmp4[, 2], sep  = "")])
    Vi2.2 <- unname(Vij[paste(tmp4[, 1], tmp4[, 3], sep  = "")])
    Vi2.3 <- unname(Vij[paste(tmp4[, 1], tmp4[, 4], sep  = "")])
    Vi2.4 <- unname(Vij[paste(tmp4[, 2], tmp4[, 3], sep  = "")])
    Vi2.5 <- unname(Vij[paste(tmp4[, 2], tmp4[, 4], sep  = "")])
    Vi2.6 <- unname(Vij[paste(tmp4[, 3], tmp4[, 4], sep  = "")])
    
    mat4 <- cbind(Vi4, Vi3.1, Vi3.2, Vi3.3,
                  Vi2.1, Vi2.2, Vi2.3, Vi2.4, Vi2.5,
                  Vi2.6, do.call(rbind, Vi1))
    
    Vijlm <- unname(apply(mat4, 1, function(x) Reduce("-", x)))
    
    if (first == "azzini") {
      Sijlm <- Vijlm / VY[(length(params) + length(com2) + length(com3) + 1):
                            (length(params) + length(com2) + length(com3) + length(com4))]
      
    } else {
      Sijlm <- Vijlm / VY
    }
    
  } else {
    Sijlm <- NULL
  }
  return(c(Si, Ti, Sij, Sijl, Sijlm))
}


bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])
    # confidence interval
    
    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }
      
    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }
      
    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }
      
    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(tmp)
}
