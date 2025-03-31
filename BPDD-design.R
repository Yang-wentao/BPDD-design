################################################################################
#  R code used in the manuscript entitled 'A Bayesian phase I/II platform design 
#  with data augmentation accounting for delayed outcomes' 
#  Author      : Wentao Yang, Zhangsheng Yu, Rongji Mu
#  R version 4.4.2 was used. R packages 'dlm', 'msm' and 'mvtnorm' are needed.
#
#  Contact     : Wentao Yang: jasper.yang@sjtu.edu.cn
#  Last update : 2025/3/31
################################################################################

rm(list = ls())
library(dlm);library(msm);library(mvtnorm)

##### Function Definition
#################

# norm2d function
norm2d <- function(v,rho,n) rmvnorm(n=n,mean=v,sigma=matrix(c(1,rho,rho,1),nrow=2))

# effdose function:
# standardizes dose levels and specifies priors for vbeta and vgamma
effdose <- function(prior.pt.mat, prior.pe.mat, ngroup, vd, type) {
  drange <- range(vd)
  coef.t <- c((drange[2] + drange[1]) / 2 / (drange[1] - drange[2]), 1 / (drange[2] - drange[1]))
  nd <- NULL
  u.mat <- array(list(), dim = c(1, 1, ngroup)) 
  vbeta <- vgamma <- matrix(0, ngroup, 2)
  for (k in 1:ngroup) {
    nd[k] <- length(prior.pt.mat[[k]])
    for (j in 1:nd[k]) {
      u.mat[[k]][j] <- round(coef.t[1] + coef.t[2] * vd[[k]][j], 2)
    }
    
    vbeta[k, ] <- lm(qnorm(prior.pt.mat[[k]]) ~ u.mat[[k]])$coefficients
    if (type == 'binary') {
      vgamma[k, ] <- lm(
      qnorm(prior.pe.mat[[k]]) ~ I(-(u.mat[[k]] - 0.1)^2)
      )$coefficients
    } else {
      vgamma[k, ] <- lm(
        prior.pe.mat[[k]] ~ I(-(u.mat[[k]] - 0.1)^2)
      )$coefficients
    }
  }

  list(
    umat = u.mat, 
    nd = nd, 
    vbeta = apply(vbeta, 2, mean), 
    vgamma = apply(vgamma, 2, mean) 
  )
  # The return value of effdose is a list containing the following elements:
  # 1. umat: a list of K numeric vectors; the k-th vector contains n_k standardized dose values for indication k.
  # 2. nd: an integer vector of length K; the k-th element represents the number of doses for indication k.
  # 3. vbeta: a 2-dimensional numeric vector containing the prior means beta0 and beta1.
  # 4. vgamma: a 2-dimensional numeric vector containing the prior means gamma0 and gamma1.
}

# tox.prob and eff.prob functions: given the matrix.uti, 
# compute the toxicity and efficacy probabilities for all dose levels under a fixed indication k.
tox.prob <- function(vp,u) {
  beta0 = vp[1]
  beta1 = vp[2]
  Z_T_mean <- beta0+beta1*u
  return(pnorm(Z_T_mean))
}

eff.prob <- function(vp,u,type) {
  c_Y = 2
  gamma0 = vp[3]
  gamma1 = vp[4]
  alpha = vp[5]
  Z_E_mean <- gamma0 - gamma1*(u-alpha)^2
  if(type=='binary'){
    pe = pnorm(Z_E_mean)
  }else{
    Z_E.sample <- rnorm(10000, Z_E_mean, sd=1)
    pe = mean(Z_E.sample>c_Y)
  }
  return(pe)
}

# get.utility function: computes the utility scores for all dose levels under a given indication k; 
# this is a helper function used within summary.mcmc.
get.utility <- function(vp,u,type) {
  beta0 = vp[1]  
  beta1 = vp[2]  
  gamma0 = vp[3]
  gamma1 = vp[4]
  alpha = vp[5]
  rho = vp[6]
  
  # V0
  Uti <- rbind(c(0, 50),
               c(25, 100))
  # V1
  # Uti <- rbind(c(0, 55),
  #              c(20, 100))
  
  # V2
  # Uti <- rbind(c(0, 66),
  #              c(33, 100))
  
  n.y <- 10000
  prob <- matrix(0, nrow = 2, ncol = 2)
  mu.T <- beta0+beta1*u 
  mu.E <- gamma0 - gamma1 * (u - alpha)^2
  sample2d <- norm2d(v=c(mu.T, mu.E), rho=rho, n=n.y)
  c_Y = ifelse(type=='binary', 0, y_E.cut)
  prob[1,1] <- sum((sample2d[,1]>0) & (sample2d[,2]<=c_Y))/n.y
  prob[1,2] <- sum((sample2d[,1]>0) & (sample2d[,2]>c_Y))/n.y
  prob[2,1] <- sum((sample2d[,1]<=0) & (sample2d[,2]<=c_Y))/n.y
  prob[2,2] <- sum((sample2d[,1]<=0) & (sample2d[,2]>c_Y))/n.y
  uti <- sum(Uti*prob)
  return(uti)
}

# summary.mcmc function: returns a cbind of (proportion of samples with toxicity probability below phi.T, 
# efficacy probability above phi.E, and utility score).
summary.mcmc <- function(matrix.uti,n.dose,u,type) {
  
  mean.matrix.uti <- apply(matrix.uti,2,mean)
  uti <- rep(0,n.dose)
  tox.mcmc <- rep(0,n.dose)
  eff.mcmc <- rep(0,n.dose)
  for (i in 1:n.dose) {
    tox.prob.mcmc <- apply(matrix.uti,1,tox.prob,u=u[i])
    eff.prob.mcmc <- apply(matrix.uti,1,eff.prob,u=u[i],type)
    tox.mcmc[i] <- sum(tox.prob.mcmc<phi.T)/nrow(matrix.uti)
    eff.mcmc[i] <- sum(eff.prob.mcmc>phi.E)/nrow(matrix.uti)
    uti[i] <- get.utility(mean.matrix.uti,u=u[i],type)
  }
  list(tox=tox.mcmc, eff=eff.mcmc, uti=uti)
}

# WT and WE are transformation functions, and invWT and invWE are their inverse functions 
# constructed via root-finding; they are used to simulate tT and tE for patients with toxicity or efficacy.
W_T_true <- function(tilde_t, u_i, cT=cT_true, eta_0T=eta_0T_true, eta_1T=eta_1T_true) {
  return(pnorm(cT * log(tilde_t) + eta_0T + eta_1T * u_i) / pnorm(eta_0T + eta_1T * u_i))
}
W_E_true <- function(tilde_t, u_i, cE=cE_true, eta_0E=eta_0E_true, eta_1E=eta_1E_true) {
  return(pnorm(cE * log(tilde_t) + eta_0E + eta_1E * u_i) / pnorm(eta_0E + eta_1E * u_i))
}
inv_W_T <- function(p, u_i) {
  uniroot(function(tilde_t) W_T_true(tilde_t, u_i) - p, interval = c(1e-6, 1))$root
}
inv_W_E <- function(p, u_i) {
  uniroot(function(tilde_t) W_E_true(tilde_t, u_i) - p, interval = c(1e-6, 1))$root
}

# outcome function: assigns yT, yE, tT, and tE for each patient 
# using the true parameter values and the W functions defined above.
outcome <- function(dose.ind, cohortsize, pt, vpara, u, type,day = 1) {
  
  gamma0.true = vpara[1]; 
  gamma1.true = vpara[2]; 
  alpha.true = vpara[3] 
  rho.true = vpara[4]    
  
  if(TRUE){Z_T_mean <- qnorm(pt[dose.ind])}
  Z_E_mean <- gamma0.true - gamma1.true * (u[dose.ind] - alpha.true)^2
  
  # Here, u represents the standardized doses for indication k, provided as u = umat[[k]].
  # In the outcome function, it is used as follows:
  # out.temp <- outcome(dose.ind[k], cohortsize = 1, pt.mat[[k]], vpt[k, ], umat[[k]], type, day = test_day)
  
  if(TRUE){
    Z <- rmvnorm(cohortsize, mean = c(Z_T_mean, Z_E_mean),
                 sigma = matrix(c(1,rho.true, rho.true,1),nrow=2)) 
    Y_T <- ifelse(Z[,1]>0, 1, 0)
    if(type=='binary'){
      Y_E <- ifelse(Z[,2]>0, 1, 0)
    }else{
      Y_E <- Z[,2]
    }
  }
  
  # assign toxicity time `t^T` and efficacy time `t^E` for each patient separately
  random_vals_T <- runif(cohortsize)
  random_vals_E <- runif(cohortsize)
  
  T_T <- round(max(5,ifelse(Y_T == 1, 
                sapply(random_vals_T, function(p) U_T * inv_W_T(p, u[dose.ind])), 
                Inf)),2)
  T_E <- round(max(15,ifelse(Y_E == 1, 
                sapply(random_vals_E, function(p) U_E * inv_W_E(p, u[dose.ind])), 
                Inf)),2)
  y_T <- rep(NA, cohortsize)
  y_E <- rep(NA, cohortsize)

  t_follow <- rep(0, cohortsize)
  t_entry <- rep(day, cohortsize)

  return(cbind(Y_T, Y_E, rep(u[dose.ind], cohortsize), T_T, T_E, y_T, y_E, t_entry, t_follow))
}

# Function of distribution of time to toxicity WT(d) and wT(d)
W_T <- function(tilde_t, u_i, c, eta_0, eta_1) {
  return(pnorm(c * log(tilde_t) + eta_0 + eta_1 * u_i) / pnorm(eta_0 + eta_1 * u_i))
}

w_T <- function(tilde_t, u_i, c, eta_0, eta_1) {
  num <- dnorm(c * log(tilde_t) + eta_0 + eta_1 * u_i) * c / tilde_t
  denom <- pnorm(eta_0 + eta_1 * u_i)
  num / denom
}

#  F_func=P(yT=1)
F_func <- function(d, beta_0, beta_1) {
  mu_k <- beta_0 + beta_1 * d  # 计算均值
  return(pnorm(mu_k))          # 返回标准正态分布的累积分布值
}

# Generate the patient enrollment schedule based on a Poisson process; 
# the first patient is enrolled on day 1.
generate_enrollment_times <- function(n.cohortsize, cohortsize, ngroup, mean_intervals) {
  if (length(mean_intervals) != ngroup) {
    stop("mean_intervals 的长度必须等于组数 ngroup")
  }
  total_patients <- n.cohortsize * cohortsize
  enrollment_times <- matrix(0, nrow = total_patients + 1, ncol = ngroup)
  
  for (group in 1:ngroup) {
    intervals <- round(rexp(total_patients, rate = 1 / mean_intervals[group]), 2)
    enrollment_times[1:total_patients, group] <- cumsum(intervals)
  }
  
  min_time <- min(enrollment_times[1:total_patients, ])
  enrollment_times[1:total_patients, ] <- enrollment_times[1:total_patients, ] - min_time + 1
  
  # The last row is marked as NA and reserved for storing the outcome time of all patients.
  enrollment_times[total_patients + 1, ] <- NA
  rownames(enrollment_times) <- c(1:total_patients, "finish_day")
  
  return(enrollment_times)
}

# Merge and sort the patient enrollment schedules.
merge_and_sort_enrollment_times <- function(enrollment_times) {
  ngroup <- ncol(enrollment_times)
  long_table <- data.frame(Time = numeric(0), Group = numeric(0), PatientID = numeric(0))
  
  for (group in 1:ngroup) {
    group_times <- enrollment_times[, group]
    total_patients <- nrow(enrollment_times) - 1
    patient_ids <- c(1:total_patients, n.cohortsize * cohortsize + 1)  # 将最后一个患者的 PatientID 设为总数 + 1
    
    group_table <- data.frame(
      Time = group_times,
      Group = group,
      PatientID = patient_ids
    )
    long_table <- rbind(long_table, group_table)
  }
  
  # Sort by time; if time values are identical, further sort by Group and PatientID.
  sorted_long_table <- long_table[order(long_table$Time, long_table$Group, long_table$PatientID, na.last = TRUE), ]
  rownames(sorted_long_table) <- NULL
  
  return(sorted_long_table)
}

# Define a function to reorder enrolltimes_long_table and reassign PatientIDs.
reorder_enrolltimes_table <- function(enrolltimes_long_table) {
  fake_patient_id <- n.cohortsize * cohortsize + 1
  enrolltimes_long_table <- enrolltimes_long_table[order(
    enrolltimes_long_table$Time, 
    enrolltimes_long_table$Group, 
    enrolltimes_long_table$PatientID, 
    na.last = TRUE
  ), ]
  
  enrolltimes_long_table$PatientID <- ave(
    seq_len(nrow(enrolltimes_long_table)), 
    enrolltimes_long_table$Group, 
    FUN = function(x) {
      ifelse(
        !is.na(enrolltimes_long_table$Time[x]) & enrolltimes_long_table$PatientID[x] != fake_patient_id, 
        seq_along(x),  
        fake_patient_id  
      )
    }
  )
  return(enrolltimes_long_table)
}

# Define a function to convert a long column vector back to a 4×n matrix with row and column names.
convert_to_named_original_format <- function(long_table, ngroup, total_patients) {
  original_matrix <- matrix(NA, nrow = total_patients + 1, ncol = ngroup)
  for (group in 1:ngroup) {
    group_data <- subset(long_table, Group == group)
    original_matrix[1:nrow(group_data), group] <- group_data$Time
  }
  
  rownames(original_matrix) <- c(1:total_patients, "finish_day")
  colnames(original_matrix) <- paste("Group", 1:ngroup)

  return(original_matrix)
}

# define mcmc_BPDD_1st 和 mcmc_BPDD_2nd
mcmc_BPDD_1st <- function(dat, vbeta=vbeta, vgamma=vgamma, ngroup, type) {
  
  # MCMC 1st with GIBBS Sampling
  # dat: the input dataset containing toxicity and efficacy outcomes for each group.
  #      It is a list of K elements, where each element dat[[k]] is a matrix representing patient data for the k-th indication.
  #      dat[[k]] = matrix(nrow = cohortsize, ncol = 3)
  #        - Column 1: toxicity outcome (Y_T)
  #        - Column 2: efficacy outcome (Y_E)
  #        - Column 3: dose level (u)
  
  N.post <- 1600
  N.burnin <- 800
  ngroup = ngroup  
  
  rho.hat <- 0.2
  beta0.hat <- vbeta[1] 
  beta1.hat <- vbeta[2] 
  gamma0.hat <- vgamma[1] 
  gamma1.hat <- vgamma[2] 
  alpha.hat <- 0.05
  
  # parameters appear in toxicity model
  beta0 <- rep(beta0.hat, ngroup) 
  beta1 <- rep(beta1.hat, ngroup) 
  beta1_share <- log(beta1.hat)   
  
  mu.beta0 <- beta0.hat 
  tau2.beta0 <- (4*mu.beta0)^2 
  # tau2.beta0 <- (6*mu.beta0)^2 # V3
  mu.beta1 <- log(beta1.hat) 
  tau2.beta1 <- ((log(3)-mu.beta1)/qnorm(0.9))^2  
  
  a.beta1 <- (2*mu.beta1)^2 
  # a.beta1 <- (3*mu.beta1)^2 # V3
  sigma2.beta1 <- 0.5*a.beta1 
  
  # parameters appear in efficacy model
  rho <- rep(rho.hat, ngroup) 
  gamma0 <- rep(gamma0.hat, ngroup) 
  gamma1 <- rep(gamma1.hat, ngroup) 
  alpha <- rep(alpha.hat, ngroup)   
  gamma1_share <- log(gamma1.hat)  
  
  mu.gamma0 <- gamma0.hat 
  tau2.gamma0 <- (4*mu.gamma0)^2 
  # tau2.gamma0 <- (6*mu.gamma0)^2 # V3
  mu.gamma1 <- log(gamma1.hat) 
  
  if(type=='binary'){
    tau2.gamma1 <- ((log(3)-mu.gamma1)/qnorm(0.9))^2
    a.gamma1 <- 4*mu.gamma1^2 
    # a.gamma1 <- 9*mu.gamma1^2 # V3
    sigma2.gamma1 <- 0.5*a.gamma1 
  }else{
    tau2.gamma1 <- 4*mu.gamma1^2
    a.gamma1 <- 4*mu.gamma1^2 
    sigma2.gamma1 <-  0.5*a.gamma1
  }
  
  n <- NULL
  for(k in 1:ngroup){
    n[k] = nrow(dat[[k]])
  }
  
  # N.post: used to construct Npost × K parameter matrices 
  # for all parameters to be updated during MCMC iterations.
  rho_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_share_t <- rep(0, N.post)  
  sigma2.beta1_t <- rep(0, N.post) 
  
  gamma0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  gamma1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  alpha_t <- matrix(0, nrow=N.post, ncol = ngroup) 
  gamma1_share_t <- rep(0, N.post)  
  sigma2.gamma1_t <- rep(0, N.post) 
  
  Z_T = array(list(),c(1,1,ngroup))
  Z_T_t <- array(list(), dim = c(1, 1, ngroup))
  Z_E = array(list(),c(1,1,ngroup))
  Z_E_t <- array(list(), dim = c(1, 1, ngroup)) # Y_E for normal
  
  # Assign initial values of -0.05, +0.5, or -0.5 
  # based on whether the observed y_T and y_E are missing or equal to 1, 0.
  for(k in 1:ngroup){
    Z_T[[k]] = ifelse(is.na(dat[[k]]$y_T), -0.05, ifelse(dat[[k]]$y_T==1, 0.5, -0.5)) 
    if(type=='binary'){
      Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, ifelse(dat[[k]]$y_E==1, 0.5, -0.5))
    }else{
      Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, dat[[k]]$y_E)  # Y_E
    }
  }
  
  for(ite in 1:N.post){
    xi_cut <- c(-Inf, 0, Inf)
    zeta_cut <- c(-Inf, 0, Inf)
    
    for(k in 1:ngroup){
      # rho
      log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
        Z_T.C <- Z_E.C <- NULL
        mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
        mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][,3] - alpha[k])^2
        Z_T.C <- c(Z_T.C, Z_T[[k]]-mean.Z_T)
        Z_E.C <- c(Z_E.C, Z_E[[k]]-mean.Z_E)
        return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
      }
      
      logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
      rho[k] <- arms(rho[k], logden, function(x) ((x>-0.8)*(x<0.8)), 1)
      rho_t[ite,k] <- rho[k]
      var.z <- 1-rho[k]^2 
      
      mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho[k] * (
        Z_E[[k]] - gamma0[k] + gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
      )
      
      Z_T[[k]] <- rtnorm(
        n = nrow(dat[[k]]),
        mean = mean.Z.T12,
        sd = sqrt(var.z),
        lower = ifelse(!is.na(dat[[k]]$y_T), xi_cut[dat[[k]]$y_T + 1], xi_cut[1]),
        upper = ifelse(!is.na(dat[[k]]$y_T), xi_cut[dat[[k]]$y_T + 2], xi_cut[3])
      )

      mean.Z.E21 <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2 + 
        rho[k] * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][, 3])
      
      Z_E[[k]] <- rtnorm(
        n = nrow(dat[[k]]),
        mean = mean.Z.E21,
        sd = sqrt(var.z),
        lower = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 1], xi_cut[1]),
        upper = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 2], xi_cut[3])
      )
      
      # beta0k 
      mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
      mu.n.beta0 <- mean(Z_T[[k]]-beta1[k]*dat[[k]][,3]-rho[k]*(Z_E[[k]]-mean.Z_E)) 
      sigma2.n.beta0 <- (1-rho[k]^2)/n[k] 
      
      mean.beta0 <- (tau2.beta0*mu.n.beta0+sigma2.n.beta0*mu.beta0)/(tau2.beta0+sigma2.n.beta0) 
      var.beta0 <- tau2.beta0*sigma2.n.beta0/(tau2.beta0+sigma2.n.beta0) 
      beta0[k]  <- rnorm(1, mean = mean.beta0, sd=sqrt(var.beta0))
      beta0_t[ite,k] <- beta0[k] 
      
      # beta1k
      log.likelihood.beta1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.beta1, beta1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3] 
        mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
        Z_T.C <- Z_T-mean.Z_T 
        Z_E.C <- Z_E-mean.Z_E 
        return(-1/2/(1-rho^2)*sum(Z_T.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.beta1*(log(beta1)-beta1_share)^2)
      }
      
      logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho[k], sigma2.beta1, beta1_share,
                                                                 dat[[k]], Z_T[[k]], Z_E[[k]])-log(x) 
      beta1[k]  <- arms(beta1[k], logden_beta1k, function(x) ((x>0)*(x<3)), 1)
      beta1_t[ite,k] <- beta1[k]
      
      # gamma0k
      mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
      vstar <- -(dat[[k]][,3]-alpha[k])^2
      mu.n.gamma0 <- mean(Z_E[[k]]-gamma1[k]*vstar-rho[k]*(Z_T[[k]]-mean.Z_T))
      sigma2.n.gamma0 <-(1-rho[k]^2)/n[k]
      
      mean.gamma0 <- (tau2.gamma0*mu.n.gamma0+sigma2.n.gamma0*mu.gamma0)/(tau2.gamma0+sigma2.n.gamma0)
      var.gamma0 <- tau2.gamma0*sigma2.n.gamma0 /(tau2.gamma0+sigma2.n.gamma0)
      gamma0[k]  <- rnorm(1, mean = mean.gamma0, sd=sqrt(var.gamma0))
      gamma0_t[ite,k] <- gamma0[k]
      
      # gamma1k
      log.likelihood.gamma1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.gamma1, gamma1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3]
        mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
        Z_T.C <- Z_T-mean.Z_T
        Z_E.C <- Z_E-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.gamma1*(log(gamma1)-gamma1_share)^2)
      }
      
      logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho[k], sigma2.gamma1, gamma1_share,
                                                                   dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)
      gamma1.Star <- ifelse(type=='binary', 5, 5) 
      gamma1[k]  <- arms(gamma1[k], logden_gamma1k, function(x) ((x>0)*(x<gamma1.Star)), 1)
      gamma1_t[ite,k] <- gamma1[k]
      
      
      # alphak
      log.likelihood.alpha <- function(alphak, rho, beta0k, beta1k, gamma0k, gamma1k, Z_Tk, Z_Ek, datk) {
        u = datk[,3] 
        mean.Z_T <- beta0k + beta1k*u
        mean.Z_E <- gamma0k - gamma1k * (u - alphak)^2
        Z_T.C <- Z_Tk-mean.Z_T
        Z_E.C <- Z_Ek-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C))
      }
      logden_alphak <- function(x) log.likelihood.alpha(x, rho[k], beta0[k], beta1[k], gamma0[k], gamma1[k], Z_T[[k]], Z_E[[k]], dat[[k]])
      
      alpha[k] <- arms(alpha[k],logden_alphak, function(x) ((x>-0.5)*(x<0.5)), 1)
      alpha_t[ite,k] <- alpha[k]
  
    }
    
    # Note: In each iteration, the shared information layer should be updated only after completing updates 
    # for the apparent layers of all K groups (k = 1 to K).
    
    # beta1_share
    mu.n.beta1_share <- mean(log(beta1))
    sigma2.n.beta1_share <- sigma2.beta1/ngroup  
    mean.beta1_share <- (tau2.beta1*mu.n.beta1_share+sigma2.n.beta1_share*mu.beta1)/(tau2.beta1+sigma2.n.beta1_share)
    var.beta1_share <- (tau2.beta1*sigma2.n.beta1_share)/(tau2.beta1+sigma2.n.beta1_share)
    
    beta1_share <- rnorm(1, mean=mean.beta1_share, sd = sqrt(var.beta1_share))
    beta1_share_t[ite] <- beta1_share
    
    # sigma2.beta1
    logden_sigma2.beta1 <- function(x) -ngroup/2*log(x)- sum((log(beta1)-beta1_share)^2)/(2*x)
    sigma2.beta1 <- arms(sigma2.beta1, logden_sigma2.beta1, function(x) (x>0)*(x<a.beta1), 1)
    sigma2.beta1_t[ite] <- sigma2.beta1
    
    # gamma1_share
    mu.n.gamma1_share <- mean(log(gamma1))
    sigma2.n.gamma1_share <- sigma2.gamma1/ngroup
    mean.gamma1_share <- (tau2.gamma1*mu.n.gamma1_share + sigma2.n.gamma1_share*mu.gamma1)/(tau2.gamma1+sigma2.n.gamma1_share)
    var.gamma1_share <- (tau2.gamma1*sigma2.n.gamma1_share)/(tau2.gamma1+sigma2.n.gamma1_share)
    
    gamma1_share <- rnorm(1, mean = mean.gamma1_share, sd=sqrt(var.gamma1_share))
    gamma1_share_t[ite] <- gamma1_share
    
    # sigma2.gamma1
    logden_sigma2.gamma1 <- function(x) -ngroup/2*log(x) - sum((log(gamma1)-gamma1_share)^2)/(2*x)
    sigma2.gamma1 <- arms(sigma2.gamma1, logden_sigma2.gamma1, function(x) (x>0)*(x<a.gamma1), 1)
    sigma2.gamma1_t[ite] <- sigma2.gamma1
    
    # This marks the completion of one full iteration; 
    # the process is repeated for a total of N.post iterations.
  }
  
  # Define an index sequence to select stable posterior samples after the burn-in period.
  ind <- seq((N.burnin+1),N.post)
  matrix.uti <- array(list(),dim=c(1,1,ngroup))
  
  for(k in 1:ngroup){
    matrix.uti[[k]] = cbind(beta0_t[,k], beta1_t[,k], gamma0_t[,k],
                            gamma1_t[,k], alpha_t[,k], rho_t[,k])[ind, ]
    
    # matrix.uti is a list, where each element matrix.uti[[k]] is a parameter matrix for the k-th indication.
    # Each matrix.uti[[k]] stores the MCMC samples of all apparent layer parameters for indication k across N.post iterations.
    # For example, beta0_t is an N.post × K matrix that records the values of beta0 for each indication k over all iterations.
  }  
  return(matrix.uti) 
}

mcmc_BPDD_2nd <- function(dat, vbeta=vbeta, vgamma=vgamma, ngroup, type, c_hat, eta_0_hat, eta_1_hat) {
  
  # MCMC 2nd with GIBBS Sampling
  # dat: the input dataset containing toxicity and efficacy outcomes for each group.
  #      It is a list of K elements, where each element dat[[k]] is a matrix representing patient data for the k-th indication.
  #      dat[[k]] = matrix(nrow = cohortsize, ncol = 3)
  #        - Column 1: toxicity outcome (Y_T)
  #        - Column 2: efficacy outcome (Y_E)
  #        - Column 3: dose level (u)
  
  N.post <- 3000
  N.burnin <- 1000
  ngroup = ngroup   

  rho.hat <- 0.2
  beta0.hat <- vbeta[1] 
  beta1.hat <- vbeta[2] 
  gamma0.hat <- vgamma[1] 
  gamma1.hat <- vgamma[2] 
  
  alpha.hat <- 0.05
  
  # parameters appear in toxicity model
  beta0 <- rep(beta0.hat, ngroup) 
  beta1 <- rep(beta1.hat, ngroup) 
  beta1_share <- log(beta1.hat)   
  
  mu.beta0 <- beta0.hat 
  tau2.beta0 <- (4*mu.beta0)^2 
  # tau2.beta0 <- (6*mu.beta0)^2 # V3
  mu.beta1 <- log(beta1.hat) 
  tau2.beta1 <- ((log(3)-mu.beta1)/qnorm(0.9))^2  
  
  a.beta1 <- (2*mu.beta1)^2 
  # a.beta1 <- (3*mu.beta1)^2 # V3
  sigma2.beta1 <- 0.5*a.beta1
  
  
  # parameters appear in efficacy model
  rho <- rep(rho.hat, ngroup) 
  gamma0 <- rep(gamma0.hat, ngroup) 
  gamma1 <- rep(gamma1.hat, ngroup) 
  alpha <- rep(alpha.hat, ngroup)   
  gamma1_share <- log(gamma1.hat)  
  
  mu.gamma0 <- gamma0.hat 
  tau2.gamma0 <- (4*mu.gamma0)^2 
  # tau2.gamma0 <- (6*mu.gamma0)^2 # V3
  mu.gamma1 <- log(gamma1.hat)  
  
  if(type=='binary'){
    tau2.gamma1 <- ((log(3)-mu.gamma1)/qnorm(0.9))^2 
    a.gamma1 <- 4*mu.gamma1^2 
    # a.gamma1 <- 9*mu.gamma1^2 # V3
    sigma2.gamma1 <- 0.5*a.gamma1 
  }else{
    tau2.gamma1 <- 4*mu.gamma1^2
    a.gamma1 <- 4*mu.gamma1^2 
    sigma2.gamma1 <-  0.5*a.gamma1
  }

  n <- NULL
  for(k in 1:ngroup){
    n[k] = nrow(dat[[k]])
  }
  
  rho_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_share_t <- rep(0, N.post)  
  sigma2.beta1_t <- rep(0, N.post) 
  
  gamma0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  gamma1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  alpha_t <- matrix(0, nrow=N.post, ncol = ngroup) 
  gamma1_share_t <- rep(0, N.post)  
  sigma2.gamma1_t <- rep(0, N.post) 
  
  Z_T = array(list(),c(1,1,ngroup))
  Z_T_t <- array(list(), dim = c(1, 1, ngroup))
  
  
  Z_E = array(list(),c(1,1,ngroup))
  Z_E_t <- array(list(), dim = c(1, 1, ngroup)) # Y_E for normal
  
  for(k in 1:ngroup){
    Z_T[[k]] = ifelse(is.na(dat[[k]]$y_T), -0.05, ifelse(dat[[k]]$y_T==1, 0.5, -0.5)) 
    if(type=='binary'){
      Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, ifelse(dat[[k]]$y_E==1, 0.5, -0.5))
    }else{
      Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, dat[[k]]$y_E)  # Y_E
    }
  }
  
  for(ite in 1:N.post){
    xi_cut <- c(-Inf, 0, Inf)
    zeta_cut <- c(-Inf, 0, Inf)
    for(k in 1:ngroup){
      # rho
      log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
        Z_T.C <- Z_E.C <- NULL
        
        mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
        mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][,3] - alpha[k])^2
        Z_T.C <- c(Z_T.C, Z_T[[k]]-mean.Z_T)
        Z_E.C <- c(Z_E.C, Z_E[[k]]-mean.Z_E)
        
        return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
      }
      
      logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
      rho[k] <- arms(rho[k], logden, function(x) ((x>-0.8)*(x<0.8)), 1)
      rho_t[ite,k] <- rho[k]
      
      var.z <- 1-rho[k]^2 
      mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho[k] * (
        Z_E[[k]] - gamma0[k] + gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
      )
      
      # The most critical step in the second-stage MCMC!!!
      # Generate Z_T[[k]] using sapply element-wise — updated to match the generation method for the second-stage MCMC.
      Z_T[[k]] <- sapply(1:nrow(dat[[k]]), function(i) {
        if (!is.na(dat[[k]]$y_T[i])) {
          rtnorm(
            n = 1,
            mean = mean.Z.T12[i], 
            sd = sqrt(var.z),
            lower = xi_cut[dat[[k]]$y_T[i] + 1], 
            upper = xi_cut[dat[[k]]$y_T[i] + 2]
          )
        } else {
          # imputation method !
          F_val <- pnorm(beta0[k] + beta1[k] * dat[[k]][i, 3])  
          W_T_val <- W_T(dat[[k]]$t_follow[i] / U_T, dat[[k]][i, 3], c_hat, eta_0_hat, eta_1_hat) 
          p_ki <- F_val * (1 - W_T_val) / (1 - F_val * W_T_val)  
          u <- runif(1)
          if (u <= p_ki) {
            rtnorm(n = 1, mean = mean.Z.T12[i], sd = sqrt(var.z), lower = 0, upper = Inf)
          } else {
            rtnorm(n = 1, mean = mean.Z.T12[i], sd = sqrt(var.z), lower = -Inf, upper = 0)
          }
        }
      })
      
      mean.Z.E21 <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2 + 
        rho[k] * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][, 3])
      
      if (type == 'binary') {
        Z_E[[k]] <- rtnorm(
          n = nrow(dat[[k]]),
          mean = mean.Z.E21,
          sd = sqrt(var.z),
          lower = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 1], xi_cut[1]),
          upper = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 2], xi_cut[3])
        )
      } else {
        Z_E[[k]] <- mean.Z.E21
      }

      # beta0k 
      mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
      mu.n.beta0 <- mean(Z_T[[k]]-beta1[k]*dat[[k]][,3]-rho[k]*(Z_E[[k]]-mean.Z_E)) 
      sigma2.n.beta0 <- (1-rho[k]^2)/n[k]
      
      mean.beta0 <- (tau2.beta0*mu.n.beta0+sigma2.n.beta0*mu.beta0)/(tau2.beta0+sigma2.n.beta0) 
      var.beta0 <- tau2.beta0*sigma2.n.beta0/(tau2.beta0+sigma2.n.beta0) 
      beta0[k]  <- rnorm(1, mean = mean.beta0, sd=sqrt(var.beta0)) 
      beta0_t[ite,k] <- beta0[k] 
      
      # beta1k
      log.likelihood.beta1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.beta1, beta1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3] 
        mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
        Z_T.C <- Z_T-mean.Z_T 
        Z_E.C <- Z_E-mean.Z_E 
        return(-1/2/(1-rho^2)*sum(Z_T.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.beta1*(log(beta1)-beta1_share)^2)
      }
      
      logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho[k], sigma2.beta1, beta1_share,
                                                                 dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)
      
      beta1[k]  <- arms(beta1[k], logden_beta1k, function(x) ((x>0)*(x<3)), 1)
      beta1_t[ite,k] <- beta1[k]
      
      # gamma0k
      mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
      vstar <- -(dat[[k]][,3]-alpha[k])^2
      mu.n.gamma0 <- mean(Z_E[[k]]-gamma1[k]*vstar-rho[k]*(Z_T[[k]]-mean.Z_T))
      sigma2.n.gamma0 <-(1-rho[k]^2)/n[k]
      
      mean.gamma0 <- (tau2.gamma0*mu.n.gamma0+sigma2.n.gamma0*mu.gamma0)/(tau2.gamma0+sigma2.n.gamma0)
      var.gamma0 <- tau2.gamma0*sigma2.n.gamma0 /(tau2.gamma0+sigma2.n.gamma0)
      gamma0[k]  <- rnorm(1, mean = mean.gamma0, sd=sqrt(var.gamma0))
      gamma0_t[ite,k] <- gamma0[k]
      
      # gamma1k
      log.likelihood.gamma1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.gamma1, gamma1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3]
        mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
        Z_T.C <- Z_T-mean.Z_T
        Z_E.C <- Z_E-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.gamma1*(log(gamma1)-gamma1_share)^2)
      }
      
      logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho[k], sigma2.gamma1, gamma1_share,
                                                                   dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)
      
      gamma1.Star <- ifelse(type=='binary', 5, 5) 
      gamma1[k]  <- arms(gamma1[k], logden_gamma1k, function(x) ((x>0)*(x<gamma1.Star)), 1)
      gamma1_t[ite,k] <- gamma1[k]
      
      
      # alphak
      log.likelihood.alpha <- function(alphak, rho, beta0k, beta1k, gamma0k, gamma1k, Z_Tk, Z_Ek, datk) {
        u = datk[,3] 
        mean.Z_T <- beta0k + beta1k*u
        mean.Z_E <- gamma0k - gamma1k * (u - alphak)^2
        Z_T.C <- Z_Tk-mean.Z_T
        Z_E.C <- Z_Ek-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C))
      }
      logden_alphak <- function(x) log.likelihood.alpha(x, rho[k], beta0[k], beta1[k], gamma0[k], gamma1[k], Z_T[[k]], Z_E[[k]], dat[[k]])
      
      alpha[k] <- arms(alpha[k],logden_alphak, function(x) ((x>-0.5)*(x<0.5)), 1)
      alpha_t[ite,k] <- alpha[k]
      
    }
    
    # beta1_share
    mu.n.beta1_share <- mean(log(beta1))
    sigma2.n.beta1_share <- sigma2.beta1/ngroup  
    mean.beta1_share <- (tau2.beta1*mu.n.beta1_share+sigma2.n.beta1_share*mu.beta1)/(tau2.beta1+sigma2.n.beta1_share)
    var.beta1_share <- (tau2.beta1*sigma2.n.beta1_share)/(tau2.beta1+sigma2.n.beta1_share)
    
    beta1_share <- rnorm(1, mean=mean.beta1_share, sd = sqrt(var.beta1_share))
    beta1_share_t[ite] <- beta1_share
    
    # sigma2.beta1
    logden_sigma2.beta1 <- function(x) -ngroup/2*log(x)- sum((log(beta1)-beta1_share)^2)/(2*x)
    sigma2.beta1 <- arms(sigma2.beta1, logden_sigma2.beta1, function(x) (x>0)*(x<a.beta1), 1)
    sigma2.beta1_t[ite] <- sigma2.beta1
    
    # gamma1_share
    mu.n.gamma1_share <- mean(log(gamma1))
    sigma2.n.gamma1_share <- sigma2.gamma1/ngroup
    mean.gamma1_share <- (tau2.gamma1*mu.n.gamma1_share + sigma2.n.gamma1_share*mu.gamma1)/(tau2.gamma1+sigma2.n.gamma1_share)
    var.gamma1_share <- (tau2.gamma1*sigma2.n.gamma1_share)/(tau2.gamma1+sigma2.n.gamma1_share)
    
    gamma1_share <- rnorm(1, mean = mean.gamma1_share, sd=sqrt(var.gamma1_share))
    gamma1_share_t[ite] <- gamma1_share
    
    # sigma2.gamma1
    logden_sigma2.gamma1 <- function(x) -ngroup/2*log(x) - sum((log(gamma1)-gamma1_share)^2)/(2*x)
    sigma2.gamma1 <- arms(sigma2.gamma1, logden_sigma2.gamma1, function(x) (x>0)*(x<a.gamma1), 1)
    sigma2.gamma1_t[ite] <- sigma2.gamma1
    
  }
  
  ind <- seq((N.burnin+1),N.post)
  matrix.uti <- array(list(),dim=c(1,1,ngroup))
  
  for(k in 1:ngroup){
    matrix.uti[[k]] = cbind(beta0_t[,k], beta1_t[,k], gamma0_t[,k],
                            gamma1_t[,k], alpha_t[,k], rho_t[,k])[ind, ]
  }  
  return(matrix.uti) 
}

#################



###### True value assignment section
#################

cT_true=1; eta_0T_true=0; eta_1T_true=2;
cE_true=2; eta_0E_true=0; eta_1E_true=2;

# V4; delay time to events
# cT_true=1; eta_0T_true=-0.5; eta_1T_true=2;
# cE_true=2; eta_0E_true=-0.5; eta_1E_true=2;

U_T=30; U_E=90;
C_T = 0.05; C_E = 0.05;

ngroup = 4 # number of indications
pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialization

cohortsize = 3 
phi.T = 0.3 
phi.E = 0.3 
type = 'binary' 

######## scenario1234
# scenario1 OBD:1234
n.cohortsize = 14
vd <- list(
  c(0.35, 0.7, 0.80, 0.90),
  c(0.20, 0.49, 0.73, 0.81),
  c(0.10, 0.24, 0.50, 0.7),
  c(0.11, 0.18, 0.28, 0.52)
)

pt.mat[[1]] = c(0.18, 0.42, 0.50, 0.58)
pt.mat[[2]] = c(0.14, 0.27, 0.41, 0.46)
pt.mat[[3]] = c(0.12, 0.17, 0.30, 0.44)
pt.mat[[4]] = c(0.10, 0.12, 0.15, 0.25)


vpt <- rbind(
  c(-0.22, 1, 0.1, 0.5),
  c(-0.25, 3, 0.1, 0.5),
  c(-0.1, 2.4, 0.1, 0.5),
  c(-0.1, 4, 0.1, 0.5)
)

# # scenario2 OBD:22(23)3
# n.cohortsize = 14
# vd <- list(
#   c(0.20, 0.45, 0.72, 0.85),
#   c(0.20, 0.53, 0.73, 0.90),
#   c(0.19, 0.45, 0.61, 0.8),
#   c(0.21, 0.35, 0.5, 0.8)
# )
# 
# pt.mat[[1]] = c(0.08, 0.21, 0.42, 0.54)
# pt.mat[[2]] = c(0.10, 0.27, 0.41, 0.54)
# pt.mat[[3]] = c(0.11, 0.22, 0.30, 0.45)
# pt.mat[[4]] = c(0.10, 0.15, 0.22, 0.39)
# 
# 
# vpt <- rbind(
#   c(-0.1, 3, 0.1, 0.5),
#   c(-0.25, 3, 0.1, 0.5),
#   c(-0.2, 3, 0.1, 0.5),
#   c(-0.15, 4, 0.1, 0.5)
# )


# # scenario3 OBD:0120
# n.cohortsize = 14
# vd <- list(
#   c(0.20, 0.45, 0.65),
#   c(0.3, 0.65, 0.85),
#   c(0.10, 0.4, 0.75, 0.9),
#   c(0.1, 0.3, 0.6, 0.85)
# )
# 
# pt.mat[[1]] = c(0.13, 0.25, 0.37)
# pt.mat[[2]] = c(0.16, 0.38, 0.54)
# pt.mat[[3]] = c(0.11, 0.25, 0.49, 0.60)
# pt.mat[[4]] = c(0.18, 0.27, 0.44, 0.59)
# 
# vpt <- rbind(
#   c(-1.1, 1, 0.1, 0.5),
#   c(-0.4, 0.4, 0.1, 0.5),
#   c(-0.3, 2, 0.1, 0.5),
#   c(-1.18, 0.8, 0.1, 0.5)
# )


# # scenario4 OBD:0(12)23
# n.cohortsize = 12
# vd <- list(
#   c(0.3, 0.5, 0.75),
#   c(0.42, 0.58, 0.9),
#   c(0.20, 0.45, 0.75),
#   c(0.1, 0.25, 0.54)
# )
# 
# pt.mat[[1]] = c(0.21, 0.34, 0.54)
# pt.mat[[2]] = c(0.22, 0.27, 0.38)
# pt.mat[[3]] = c(0.13, 0.25, 0.43)
# pt.mat[[4]] = c(0.14, 0.18, 0.29)
# 
# vpt <- rbind(
#   c(-0.9, 1, 0.1, 0.5),
#   c(-0.32, 3, 0.1, 0.5),
#   c(-0.3, 3, 0.1, 0.5),
#   c(-0.25,3.6, 0.1, 0.5)
# )

########

# Initialize prior.pt.mat and prior.pe.mat
prior.pt.mat <- vector("list", length(vd))
prior.pe.mat <- vector("list", length(vd))

# # Prior values for toxicity and efficacy
for (i in 1:ngroup){
  if(length(vd[[i]])==4){
    prior.pt.mat[[i]] = c(0.10, 0.20, 0.30, 0.40)
    prior.pe.mat[[i]] = c(0.15, 0.30, 0.45, 0.40)
  } else if (length(vd[[i]])==3) {
    prior.pt.mat[[i]] = c(0.10, 0.20, 0.30)
    prior.pe.mat[[i]] = c(0.15, 0.30, 0.45)
  }
}

#################

# Start the main process via the main_BPDD() function
main_BPDD <- function(n.cohortsize, cohortsize, phi.T, phi.E, pt.mat, prior.pt.mat,
                              prior.pe.mat, vpt, vd, type, C_T = 0.05, C_E = 0.05) {

  # Preprocessing of true values and initialization of related variables
  #####################
  ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
  umat <- ed$umat
  nd <- ed$nd
  vbeta <- ed$vbeta  
  vgamma <- ed$vgamma 
  
  k.status <- rep(1, ngroup)    # Status indicator for each group: 1 = active, 0 = stopped
  dose.ind <- rep(1, ngroup)    # Current dose index for each group
  h.dose <- rep(1, ngroup)      # Highest dose level tried in each group
  dat <- vector("list", ngroup) # A list to store patient-level data for each group
  finish_day <- rep(1, ngroup)  # The trial completion day for each group
  
  #####################
  
  # Generate the enrollment time table: first create 4 × Nk entries per group, 
  # then merge them into a long-format table where patients are enrolled and assessed sequentially.
  # Patient enrollment follows a Poisson process,
  # with average inter-arrival times of 10, 12, 14, and 16 days, respectively.
  enrollment_times <- generate_enrollment_times(n.cohortsize, cohortsize, ngroup, mean_intervals=c(10, 12, 14, 16))
  enrolltimes_long_table <- merge_and_sort_enrollment_times(enrollment_times)
  total_patients <- n.cohortsize * cohortsize
  test_day = 0
  
  ### trials start
  row <- 1  
  while (row <= nrow(enrolltimes_long_table)){
    test_day <- enrolltimes_long_table[row, 1]
    j <- enrolltimes_long_table[row, 2]  
    patient_id <- enrolltimes_long_table[row, 3] 
    
    if (k.status[j]==0){ row = row + 1 ; next }
    
    if (patient_id == 4) { # deal with the second cohort
      
      # Update the toxicity (y_T) and efficacy (y_E) outcomes for the first three patients in group j
      if (!is.null(dat[[j]]) && nrow(dat[[j]]) > 0) {
        dat[[j]]$t_follow <- test_day - dat[[j]]$t_entry
        dat[[j]]$y_T <- ifelse(
          dat[[j]]$Y_T == 1 & dat[[j]]$T_T <= dat[[j]]$t_follow + 1e-5 , 1,
          ifelse(dat[[j]]$Y_T == 0 & dat[[j]]$t_follow + 1e-5 >= U_T, 0, NA)
        )
        dat[[j]]$y_E <- ifelse(
          dat[[j]]$Y_E == 1 & dat[[j]]$T_E <= dat[[j]]$t_follow + 1e-5 , 1,
          ifelse(dat[[j]]$Y_E == 0 & dat[[j]]$t_follow + 1e-5 >= U_E, 0, NA)
        )
      }

      previous_patients <- dat[[j]][1:3, ]
      toxicity_2patients_no_occurred <- sum(previous_patients$y_T == 0, na.rm = TRUE) >= 2
      toxicity_occurred <- sum(previous_patients$y_T == 1, na.rm = TRUE) >= 1
      
      if (!toxicity_2patients_no_occurred && !toxicity_occurred) {
        completed_times <- ifelse(
          previous_patients$Y_T == 1,
          previous_patients$t_entry + previous_patients$T_T,  
          previous_patients$t_entry + U_T                     
        )
        
        # Determine delay times based on different toxicity patterns: (000), (001), (011), (111)
        if (sum(previous_patients$Y_T == 0) == 3) {
          sorted_times <- sort(completed_times, decreasing = TRUE)
          delay_time <- sorted_times[2]  
        } else if (sum(previous_patients$Y_T == 0) == 2 && sum(previous_patients$Y_T == 1) == 1) {
          max_no_toxicity_time <- max(completed_times[previous_patients$Y_T == 0])
          single_toxicity_time <- min(completed_times[previous_patients$Y_T == 1])
          delay_time <- min(max_no_toxicity_time, single_toxicity_time)
        } else if (sum(previous_patients$Y_T == 1) == 2) {
          delay_time <- min(completed_times[previous_patients$Y_T == 1])
        } else if (sum(previous_patients$Y_T == 1) == 3) {
          delay_time <- min(completed_times[previous_patients$Y_T == 1])
        } else {
          stop("Unexpected case: Toxicity data does not fit expected patterns.")
        }
        
        enrolltimes_long_table[row,1] <- delay_time + 0.05
        
        # After delaying enrollment times, reorder the table and reassign PatientIDs
        enrolltimes_long_table <- reorder_enrolltimes_table(enrolltimes_long_table)
        next  
      } else { 
        if (toxicity_occurred) {
          dose.ind[j] <- 1
        } else if (toxicity_2patients_no_occurred) {
          dose.ind[j] <- min(h.dose[j] + 1, nd[j])  
          h.dose[j] <- dose.ind[j]
        } 
      }  
      
      } else if (patient_id != 1 && patient_id %% cohortsize == 1) {
        data_no_yT_NA <- vector("list", ngroup) 
        dat_simulation <- vector("list", ngroup)
        for (kk in 1:ngroup){
          if(!is.null(dat[[kk]]) && nrow(dat[[kk]]) > 0) {
            dat[[kk]]$t_follow <- test_day - dat[[kk]]$t_entry
            dat[[kk]]$y_T <- ifelse(
              dat[[kk]]$Y_T == 1 & dat[[kk]]$T_T <= dat[[kk]]$t_follow + 1e-5, 1,
              ifelse(dat[[kk]]$Y_T == 0 & dat[[kk]]$t_follow + 1e-5 >= U_T, 0, NA)
            )
            dat[[kk]]$y_E <- ifelse(
              dat[[kk]]$Y_E == 1 & dat[[kk]]$T_E <= dat[[kk]]$t_follow + 1e-5, 1,
              ifelse(dat[[kk]]$Y_E == 0 & dat[[kk]]$t_follow + 1e-5 >= U_E, 0, NA)
            )
            data_no_yT_NA[[kk]] <- dat[[kk]][!is.na(dat[[kk]]$y_T), ]
            dat_simulation[[kk]] <- dat[[kk]]
          }
        }  
        
        for (kk in 1:ngroup){
          if(is.null(dat[[kk]]) || nrow(dat[[kk]]) == 0){
            u <- umat[[kk]]
            out.temp <- outcome(dose.ind[kk], cohortsize = 1, pt.mat[[kk]], vpt[kk, ], u, type, day = test_day)
            if (is.null(out.temp) || nrow(out.temp) == 0) {
              stop("out.temp is empty or unavailable")
            }
            out.temp <- as.data.frame(out.temp)
            colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")
            
            out.temp$group <- kk
            out.temp$Y_T = 0                      
            out.temp$Y_E = 0                      
            out.temp$dose = dose.ind[kk]          
            out.temp$T_T = Inf                    
            out.temp$T_E = Inf                    
            out.temp$y_T = 0                      
            out.temp$y_E = 0                      
            out.temp$t_entry = test_day-U_T       
            out.temp$t_follow = U_T               
            
            dat_simulation[[kk]] <- out.temp
          }
        }  
        
        all_elements_not_null_nrow_more_than0 <- all(sapply(data_no_yT_NA, function(x) !is.null(x) && nrow(x) > 0))
        
        for (kk in 1:ngroup){
          if(is.null(data_no_yT_NA[[kk]]) || nrow(data_no_yT_NA[[kk]]) == 0){
            u <- umat[[kk]]
            out.temp <- outcome(dose.ind[kk], cohortsize = 1, pt.mat[[kk]], vpt[kk, ], u, type, day = test_day)
            if (is.null(out.temp) || nrow(out.temp) == 0) {
              stop("out.temp is empty or unavailable")
            }
            out.temp <- as.data.frame(out.temp)
            colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")
            
            out.temp$group <- kk
            out.temp$Y_T = 0                      
            out.temp$Y_E = 0                     
            out.temp$dose = dose.ind[kk]          
            out.temp$T_T = Inf                    
            out.temp$T_E = Inf                   
            out.temp$y_T = 0                      
            out.temp$y_E = 0                      
            out.temp$t_entry = test_day-U_T      
            out.temp$t_follow = U_T               
            
            data_no_yT_NA[[kk]] <- out.temp
          }
        }
        
        all_dat <- do.call(rbind, dat_simulation)

        # Step 1: Perform MCMC on data_no_yT_NA to obtain mcmc.temp_1st
        mcmc.temp_1st <- mcmc_BPDD_1st(data_no_yT_NA, vbeta=vbeta, vgamma=vgamma, ngroup=ngroup, type)
        matrix.uti <- mcmc.temp_1st  
        beta_0_1st <- numeric(ngroup) 
        beta_1_1st <- numeric(ngroup)
        
        for (kk in 1:ngroup) { 
          if (!is.null(matrix.uti[[kk]]) && nrow(matrix.uti[[kk]]) > 0) {
            beta_0_1st[kk] <- mean(matrix.uti[[kk]][, 1])  
            beta_1_1st[kk] <- mean(matrix.uti[[kk]][, 2])  
          } else {
            warning(paste("mcmc.temp_1st ", kk, "indication matrix.uti[[kk]] is empty or unavailable"))
          }
        }
        
        # Step 2: Specify and solve the maximum likelihood function
        ################
        # Likelihood function based on all toxicity-related data from enrolled patients
        
        log_L <- function(params, all_dat, beta_0_1st, beta_1_1st) {
          if (length(params) != 3) stop("params must have 3 elements: c, eta_0, eta_1")
          if (!all(c("y_T", "T_T", "t_follow", "dose", "group") %in% names(all_dat))) 
            stop("all_dat must contain columns: y_T, T_T, t_follow, dose, group")
          
          c <- params[1]
          eta_0 <- params[2]
          eta_1 <- params[3]
          
          y_T <- all_dat$y_T
          tilde_T_T <- all_dat$T_T / U_T
          tilde_t_follow <- all_dat$t_follow / U_T 
          d <- all_dat$dose
          group <- all_dat$group
          
          F_vals <- pnorm(beta_0_1st[group] + beta_1_1st[group] * d)
          w_T_vals <- w_T(tilde_T_T, d, c, eta_0, eta_1)
          W_T_vals <- W_T(tilde_t_follow, d, c, eta_0, eta_1)
          
          if (length(w_T_vals) != length(y_T) || length(W_T_vals) != length(y_T) || length(F_vals) != length(y_T)) {
            warning(paste("Vector lengths of w_T_vals, W_T_vals, F_vals, and y_T must match"))
          }
          
          # Case (1,1): A toxicity response is observed; compute g = wT * F
          log_g <- rep(0, length(y_T))
          valid_idx <- !is.na(y_T) & y_T == 1 
          log_g[valid_idx] <- log(w_T_vals[valid_idx] * F_vals[valid_idx])
    
          # Case (0,0): Toxicity unknown → 1 - W^T * F
          log_survival <- rep(0, length(y_T))
          log_survival[is.na(y_T)] <- log(1 - W_T_vals[is.na(y_T)] * F_vals[is.na(y_T)]) 
          
          log_L_total <- sum(log_g+log_survival, na.rm = TRUE) 
          penalty = 0
          return(log_L_total-penalty)
        }
        
        # Maximize log_L to estimate c, eta_0, and eta_1
        optim_res <- optim(
          par = c(1.5, 0, 1.5),  
          fn = function(params) -log_L(params, all_dat, beta_0_1st, beta_1_1st),  
          method = "L-BFGS-B", 
          lower = c(0.1, -1, 0.5), 
          upper = c(4, 2, 4),
        )
        c_hat <- optim_res$par[1]
        eta_0_hat <- optim_res$par[2]
        eta_1_hat <- optim_res$par[3]
        
        ################
        
        # Step 3: Perform another MCMC run using the mcmc_BPDD_2nd function to obtain mcmc.temp_2nd, 
        # which will be used in Step 4 for dose decision-making
        mcmc.temp_2nd <- mcmc_BPDD_2nd(dat_simulation, vbeta=vbeta, vgamma=vgamma, ngroup=ngroup, type,
                                                    c_hat, eta_0_hat, eta_1_hat)
        
        
        # Step 4: Based on the latest MCMC results (mcmc.temp_2nd), 
        # update the dose assignment for the next patient in group k
        u=umat[[j]]
        summ.temp<-summary.mcmc(matrix.uti=mcmc.temp_2nd[[j]],n.dose=nd[j], u=u, type)
        tox.summ <- summ.temp$tox
        
        if (nrow(dat[[j]]) == n.cohortsize * cohortsize){h.dose[j]=nd[j]} 
        if ((tox.summ[h.dose[j]]>.5) & (h.dose[j]!=nd[j])){
          dose.ind[j] <- h.dose[j]+1
          h.dose[j] <- dose.ind[j]
        }else{
          
          eff.summ <- summ.temp$eff
          capA.ind <- sort(which((tox.summ>C_T)&(eff.summ>C_E)))
          if(length(capA.ind)==0){ # admissible dose set is empty
            dose.ind[j] = 0
            k.status[j] = 0
            finish_day[j] = test_day
          }else{# admissible dose set empty is not empty
            
            # Extract the utility scores of all admissible dose levels for use in randomized dose assignment,
            # where the selection probabilities are based on these scores.
            uti <- summ.temp$uti[capA.ind] 
            uti.normalize <- uti/sum(uti)
            cumu.uti <- uti
            
            for (i in 1:length(uti)) cumu.uti[i] <- sum(uti.normalize[1:i])
            r <- runif(1,0,1)
            dose.ind[j] <- min(h.dose[j]+1,capA.ind[min(which(r<cumu.uti))])
            # During randomized assignment, doses beyond the two highest previously tried levels are not allowed.
            
            # update h.dose[j]
            if(dose.ind[j]>h.dose[j]){
              dose.ind[j] <- h.dose[j]+1
              h.dose[j] <- dose.ind[j]
            }
            
            if (nrow(dat[[j]]) == n.cohortsize * cohortsize){
              dose.ind[j] = capA.ind[which.max(uti)]
              k.status[j] = 0
              finish_day[j] = test_day
            }
          }
        }  
      }  # the end for the dose-decision of indication j
          
      if ( k.status[j] == 0) { row = row + 1; next }
      u <- umat[[j]]
      out.temp <- outcome(dose.ind[j], cohortsize = 1, pt.mat[[j]], vpt[j, ], u, type, day = test_day)
      
      if (is.null(out.temp) || nrow(out.temp) == 0) {
        stop("out.temp is empty or unavailable")
      }
      
      out.temp <- as.data.frame(out.temp)
      colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")
      out.temp$group <- j
      if (is.null(dat[[j]]) || nrow(dat[[j]]) == 0) {
        dat[[j]] <- out.temp
      } else {
        dat[[j]] <- rbind(dat[[j]], out.temp)
      }
      
      if (patient_id==n.cohortsize * cohortsize){
      j_group_efficacy_times <- ifelse(
        dat[[j]]$Y_E == 1,
        dat[[j]]$t_entry + dat[[j]]$T_E,  
        dat[[j]]$t_entry + U_E            
      )
      
      # Locate the position of the pseudo-patient in group j within enrolltimes_long_table
      
      fake_patient_row <- which(
        enrolltimes_long_table$Group == j & 
          enrolltimes_long_table$PatientID == ( n.cohortsize * cohortsize + 1 )
      )
      
      if (length(fake_patient_row) == 1) {
        enrolltimes_long_table[fake_patient_row,1] <- max(j_group_efficacy_times, na.rm = TRUE)
        enrolltimes_long_table <- reorder_enrolltimes_table(enrolltimes_long_table)
      
      } else {
        stop(paste("no fake patient record"))
      }
      
      }
      
      # Current row's patient has completed enrollment; proceed to the next row for enrollment evaluation
      row = row + 1
      
    } # The decision-making and enrollment steps for each patient in enrolltimes_long_table are now complete

  
  # The entire trial process is complete; proceed to data processing
  enrolltimes_final = convert_to_named_original_format(enrolltimes_long_table, ngroup, total_patients)
  matrix_means <- lapply(mcmc.temp_2nd, function(mat) {
    if (!is.null(mat) && nrow(mat) > 0) {
      colMeans(mat) 
    } else {
      NA 
    }
  })
  
  matrix_means_df <- do.call(rbind, lapply(seq_along(matrix_means), function(k) {
    data.frame(
      Indicator = k,
      t(matrix_means[[k]])
    )
  }))
  
  colnames(matrix_means_df) <- c(
    "Indicator",
    "beta0", "beta1", "gamma0", "gamma1", "alpha", "rho"  # 参数列名
  )
  
  
  # Output the selected OBD dose for each group k 
  # and the dose assignment sequence across all iterations
  res <- vector("list", ngroup)
  for(k in 1:ngroup){
    dtemp <- NULL
    for(id in 1:length(dat[[k]][,3])){
      dtemp[id] = vd[[k]][which(umat[[k]]==dat[[k]][,3][id])] 
    }
    res[[k]] = list(d.select = dose.ind[k], d.alloc=dtemp)
  }
  return(list(dose_pick=res,parameter_mean=matrix_means_df,finish_day=finish_day,enroll=enrolltimes_final,dat=dat))
}



##### Numerical simulation and result saving

library("tictoc")
set.seed(1213)
tic('start trials')
result <- main_BPDD(n.cohortsize, cohortsize, phi.T, phi.E, pt.mat, prior.pt.mat,
            prior.pe.mat, vpt, vd, type, C_T = 0.05, C_E = 0.05)
toc()

# Check the structure of the returned result
str(result)

# Extract and inspect each component
print(result$parameter_mean)    # View parameter means
print(result$dose_pick)         # View dose selection and allocation
print(result$finish_day)        # View trial duration for each group
print(result$enroll)            # View the enrollment schedule of all patients
print(result$dat)               # View all patient data at trial completion; each patient has 9 variables


