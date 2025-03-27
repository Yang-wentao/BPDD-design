#' Main function for BPDD: Bayesian Platform Data-augmentation for Delayed outcomes
#'
#' This function performs the main analysis for Bayesian data augmentation
#' in the context of delayed outcomes in Phase I/II clinical trials. The
#' function handles cohort sizes, treatment effects, prior distributions, and
#' other parameters to model delayed clinical outcomes. It applies Bayesian
#' methods to simulate and augment data, supporting robust analysis of clinical trial data.
#'
#' @param n.cohortsize The total number of cohorts to simulate, default is 14 or 12.
#' @param cohortsize The size of each cohort, default is 3.
#' @param phi.T The upper bound of the acceptable toxicity level in a drug trial.
#' @param phi.E The lower bound of the acceptable efficacy level in a drug trial.
#' @param pt.mat A matrix representing true toxicity values.
#' @param prior.pt.mat A matrix representing the priors of toxicity.
#' @param prior.pe.mat A matrix representing the priors of efficacy.
#' @param vpt A matrix representing the true parameters in setting models.
#' @param vd A matrix representing the initial doses from 0 to 1.
#' @param type The type outcomes, default is "binary".
#' @param C_T The minimum acceptable probability that the drug's toxicity is below the upper bound in the dose-finding process, with a default value of 0.05.
#' @param C_E The minimum acceptable probability that the drug's efficacy is above the lower bound in the dose-finding process, with a default value of 0.05.
#' @return A list containing the results of the Bayesian analysis, including list of OBD_dose_pick, parameter_estimate, finish_day, Patient Enrollment Schedule
#' @export
main_BPDD<- function(n.cohortsize, cohortsize, phi.T, phi.E, pt.mat, prior.pt.mat,
                     prior.pe.mat, vpt, vd, type="binary", C_T = 0.05, C_E = 0.05) {
  library(dlm);library(msm);library(mvtnorm)
  # 子函数的定义
  #################
  # norm2d函数
  norm2d <- function(v,rho,n) rmvnorm(n=n,mean=v,sigma=matrix(c(1,rho,rho,1),nrow=2))

  # effdose函数：标准化剂量，并给出vbeta，vgamma先验
  effdose <- function(prior.pt.mat, prior.pe.mat, ngroup, vd, type) {
    # return the standardized doses and the prior estimates for beta and gamma
    # ngroup：整数，表示指标K的数量（即不同临床试验的子组数量）

    drange <- range(vd)
    coef.t <- c((drange[2] + drange[1]) / 2 / (drange[1] - drange[2]), 1 / (drange[2] - drange[1]))
    # 抓取剂量的最大最小值，并把所有剂量标准化至-0.5～0.5

    nd <- NULL
    u.mat <- array(list(), dim = c(1, 1, ngroup)) # 我不知道为啥不写dim=ngroup
    vbeta <- vgamma <- matrix(0, ngroup, 2)
    # vbeta、vgamma是一个ngroup行，2列的窄矩阵

    for (k in 1:ngroup) {
      nd[k] <- length(prior.pt.mat[[k]])
      for (j in 1:nd[k]) {
        u.mat[[k]][j] <- round(coef.t[1] + coef.t[2] * vd[[k]][j], 2)
        # 剂量标准化，随后放入u.mat的第k个列表元素的第j元
      }

      # 线性拟合毒性模型的beta参数
      vbeta[k, ] <- lm(qnorm(prior.pt.mat[[k]]) ~ u.mat[[k]])$coefficients
      # 返回的 coefficients 顺序为：
      # 1. vbeta[k, 1]：截距（即 \beta_0k）
      # 2. vbeta[k, 2]：斜率（即 \beta_1k）

      # 二次拟合疗效模型的gamma参数
      if (type == 'binary') {
        # 将疗效的潜变量与二次函数形式拟合qnorm(prior.pe)=gamma0k+gamma1k(-(d-0.1)^2)
        vgamma[k, ] <- lm(
          qnorm(prior.pe.mat[[k]]) ~ I(-(u.mat[[k]] - 0.1)^2)
        )$coefficients
      } else {
        vgamma[k, ] <- lm(
          prior.pe.mat[[k]] ~ I(-(u.mat[[k]] - 0.1)^2)
        )$coefficients
      }
    }

    # effdose返回值为一个列表
    list(
      umat = u.mat, # 标准化剂量
      nd = nd, # 每组剂量数
      vbeta = apply(vbeta, 2, mean), # beta先验均值
      vgamma = apply(vgamma, 2, mean) # gamma先验均值
    )
    # effdose返回值为一个列表
    # 列表第一个元素umat: 有K个数组类型元素、第k个数组元素是n_k个标准化剂量值
    # 列表第二个元素nd：有K个正整数元素、第k个元素是第k个指标的剂量数
    # 列表第三个元素vbeta：2维数组，均值beta0，均值beta1，用于先验
    # 列表第四个元素vgamma：2维数组，均值gamma0，均值gamma1，用于先验
  }

  # tox.prob,eff.prob函数：已知matrix.uti的情况下，fit某一指标k下，所有剂量的毒性率和疗效率;
  # 在summary.mcmc函数中，使用到这两个函数时, vp=matrix.uti[[k]]，并且是sapply逐行运算
  tox.prob <- function(vp,u) {
    # return fit toxicity rate
    # u=umat[[k]]
    beta0 = vp[1]
    beta1 = vp[2]
    Z_T_mean <- beta0+beta1*u
    return(pnorm(Z_T_mean))
    # pnorm是概率，标准正态分布小于等于x的概率；q是分位数；d是密度值
  }
  eff.prob <- function(vp,u,type) {
    # return fit efficacy rate
    # u=umat[[k]]
    c_Y = 2
    gamma0 = vp[3]
    gamma1 = vp[4]
    alpha = vp[5]
    # gamma2 = vp[6]
    # Z_E_mean <- gamma0+gamma1*ifelse(u<=alpha, u, alpha) 12.4更改为下面，12.6更改为下下面
    # Z_E_mean <- gamma0+gamma1*(ifelse(u<=alpha, u, alpha)-gamma2*(u-alpha)*ifelse(u>alpha, 1, 0))
    Z_E_mean <- gamma0 - gamma1*(u-alpha)^2
    if(type=='binary'){
      pe = pnorm(Z_E_mean)
    }else{
      Z_E.sample <- rnorm(10000, Z_E_mean, sd=1)
      pe = mean(Z_E.sample>c_Y)
    }
    return(pe)
  }

  # get.utility函数：在已知获得某个指标k下，所有剂量的得分；这个函数是使用在summary.mcmc中的辅助函数
  get.utility <- function(vp,u,type) {
    # return the utility
    # 通过抽取10000次样本，获得某一个指标k下，所有剂量对应的得分
    # 调用的语句：uti[i] <- get.utility(mean.matrix.uti,u=u[i],type)

    if(FALSE){y_E.cut = 2}# 这句话在YE为二元变量时应该用不到
    #注意下面的这些参数，应该都是向量的形式（单行数组）
    beta0 = vp[1]
    beta1 = vp[2]
    gamma0 = vp[3]
    gamma1 = vp[4]
    alpha = vp[5]
    rho = vp[6]
    # gamma2 = vp[7]

    # 第1种标准情况
    Uti <- rbind(c(0, 50),
                 c(25, 100))
    # 第2种情况V1
    # Uti <- rbind(c(0, 55),
    #              c(20, 100))

    # 第3种情况V2
    # Uti <- rbind(c(0, 66),
    #              c(33, 100))

    n.y <- 10000
    prob <- matrix(0, nrow = 2, ncol = 2)
    mu.T <- beta0+beta1*u
    # mu.E <- gamma0+gamma1*(ifelse(u<=alpha, u, alpha)-gamma2*(u-alpha)*ifelse(u>alpha, 1, 0))
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

  # summary.mcmc函数：输出cbind(毒性概率低于phi.T、疗效高于phi.E的样本比例、得分）
  summary.mcmc <- function(matrix.uti,n.dose,u,type) {
    # 这段代码已经基本看懂，目的是在得到稳定的mcmc链的参数矩阵列表matrix.uti之后
    # 用后验抽样均值代替后验分布均值，对毒性率、疗效率、得分进行计算，维度为ndose
    # 这里的n.dose应该是指标k的剂量数，即nd[k]

    mean.matrix.uti <- apply(matrix.uti,2,mean)
    uti <- rep(0,n.dose)
    tox.mcmc <- rep(0,n.dose)
    eff.mcmc <- rep(0,n.dose)
    for (i in 1:n.dose) {
      tox.prob.mcmc <- apply(matrix.uti,1,tox.prob,u=u[i])
      eff.prob.mcmc <- apply(matrix.uti,1,eff.prob,u=u[i],type)
      # matrix.uti其实在使用过程中，就是之后提到的matrix.uti[[k]];
      # 每一行都代表烧录期之后的一次迭代数据
      # 这里是对matrix.uti的每一行（全参数的某一次迭代数据）操作，计算拟合毒性率和拟合疗效率
      # 随后就能计算出毒性率/疗效率低于或高于阈值的比例
      # 毒性率低于阈值的概率（即后验样本比例）>0.5,则后续用于算法里的剂量升级

      tox.mcmc[i] <- sum(tox.prob.mcmc<phi.T)/nrow(matrix.uti)
      #计算了剂量d_ki下，毒性概率低于phi.T的样本比例！

      eff.mcmc[i] <- sum(eff.prob.mcmc>phi.E)/nrow(matrix.uti)
      #计算了剂量d_ki下，疗效概率高于phi.E的样本比例！

      uti[i] <- get.utility(mean.matrix.uti,u=u[i],type)
    }
    list(tox=tox.mcmc, eff=eff.mcmc, uti=uti)
    # ～～summary.mcmc子函数感觉有点apply上的问题
    # ～～gpt建议对列表类型数据matrix.uti使用lapply,但是源代码本身是能跑的
    # ～～在加上我基本能理解这几行代码的含义，我决定先暂时不去动这个子函数
  }


  ### 下面是自己写的新函数or有更新的函数

  # WT,WE和invWT,invWE函数构造,即将用于模拟有毒性/有疗效患者的tT,tE
  # Inverse functions for W^T and W^E using root finding，用于生成模拟数据tT/tE
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

  # outcome函数：使用三个参数真值和上面的W函数，给患者赋值yTyEtTtE
  # 输出cbind(Y_T, Y_E, rep(u[dose.ind], cohortsize), T_T, T_E,y_T,y_E,t_entry,t_follow)一个患者对应一个9维信息向量
  outcome <- function(dose.ind, cohortsize, pt, vpara, u, type,day = 1) {
    # return the toxicity and efficacy outcomes and tT,tE

    # outcome函数在使用时是这样赋值的：
    # k是指标，out.temp <- outcome(dose.ind[k], cohortsize, pt.mat[[k]], vpt[k,], u=umat[[k]], type)
    # 函数使用时，vpara = vpt[k,]，vpt为所有K组（K行）的下面四个参数的真值矩阵
    # dose.ind是指标，表示（第k组的）第几个剂量

    gamma0.true = vpara[1];
    gamma1.true = vpara[2];
    alpha.true = vpara[3] # efficacy model (4)
    rho.true = vpara[4]   # used in the joint bivariate model
    # gamma2.true = vpara[5]

    if(TRUE){Z_T_mean <- qnorm(pt[dose.ind])}

    # 毒性真实值是直接已知的（未必很好的符合给出的“潜变量+线性模型”，之后可以验证一下线性模型是否合理，是否有必要变为2次函数）
    # beta是用来通过线性拟合ZT从而拟合毒性真实值的
    # 疗效的真实值是通过gamma的真实值带入不同剂量，然后计算得到的

    # Z_E_mean <- gamma0.true+
    #   gamma1.true*(ifelse(u[dose.ind]<= alpha.true, u[dose.ind], alpha.true)-gamma2.true*(u[dose.ind]-alpha.true)*ifelse(u[dose.ind]>alpha.true, 1, 0))
    Z_E_mean <- gamma0.true - gamma1.true * (u[dose.ind] - alpha.true)^2

    # 可见这里的u是（第k组）标准化后的剂量,使用时u=umat[[k]]，第k组的标准化剂量数据
    # 在outcome函数使用时，使用方式为：
    # out.temp <- outcome(dose.ind[k], cohortsize=1, pt.mat[[k]], vpt[k,], umat[[k]], type, day=test\_day)

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


    # 分配 toxicity time `t^T` and efficacy time `t^E` for each patient separately
    random_vals_T <- runif(cohortsize)
    random_vals_E <- runif(cohortsize)

    T_T <- round(max(5,ifelse(Y_T == 1,
                              sapply(random_vals_T, function(p) U_T * inv_W_T(p, u[dose.ind])),
                              Inf)),2)
    T_E <- round(max(15,ifelse(Y_E == 1,
                               sapply(random_vals_E, function(p) U_E * inv_W_E(p, u[dose.ind])),
                               Inf)),2)
    # 把随机生成的概率p，传入inv_W_E函数中，得到相应的tilde_t

    y_T <- rep(NA, cohortsize)

    y_E <- rep(NA, cohortsize)

    #避免引入字符型数据，从而让整个数据矩阵的数据都变成字符型
    # y_T <- ifelse(Y_T != 2, "mis",NA)
    #
    # y_E <- ifelse(Y_E != 2, "mis",NA)

    t_follow <- rep(0, cohortsize)

    t_entry <- rep(day, cohortsize)


    # Combine results
    return(cbind(Y_T, Y_E, rep(u[dose.ind], cohortsize), T_T, T_E, y_T, y_E, t_entry, t_follow))
    # 模拟输出这一组患者的真实Y_T和Y_E.
    # 还要在Y_E列右边新增两列t^T和t^E
  }

  # 定义W_T 这个函数可以接受向量输入！也可以接受tilde_t=Inf的输入！输出为1.188（不会报错）
  W_T <- function(tilde_t, u_i, c, eta_0, eta_1) {
    return(pnorm(c * log(tilde_t) + eta_0 + eta_1 * u_i) / pnorm(eta_0 + eta_1 * u_i))
  }

  # 定义w_T 这个函数可以接受向量输入！也可以接受tilde_t=Inf的输入！输出为0（不会报错）
  w_T <- function(tilde_t, u_i, c, eta_0, eta_1) {
    num <- dnorm(c * log(tilde_t) + eta_0 + eta_1 * u_i) * c / tilde_t
    denom <- pnorm(eta_0 + eta_1 * u_i)
    num / denom
  }

  # 定义 F_func
  F_func <- function(d, beta_0, beta_1) {
    mu_k <- beta_0 + beta_1 * d  # 计算均值
    return(pnorm(mu_k))          # 返回标准正态分布的累积分布值
  }


  # 生成患者poisson过程的入组时间表格,第一个患者在第一天入组
  generate_enrollment_times <- function(n.cohortsize, cohortsize, ngroup, mean_intervals) {
    # 参数检查
    if (length(mean_intervals) != ngroup) {
      stop("mean_intervals 的长度必须等于组数 ngroup")
    }

    # 总患者数
    total_patients <- n.cohortsize * cohortsize

    # 初始化表格
    enrollment_times <- matrix(0, nrow = total_patients + 1, ncol = ngroup)

    for (group in 1:ngroup) {
      # 使用不同的平均间隔生成时间间隔并累积生成时间
      intervals <- round(rexp(total_patients, rate = 1 / mean_intervals[group]), 2)
      enrollment_times[1:total_patients, group] <- cumsum(intervals)
    }

    # 找到所有生成时间的最小值
    min_time <- min(enrollment_times[1:total_patients, ])

    # 调整使最小时间为 1
    enrollment_times[1:total_patients, ] <- enrollment_times[1:total_patients, ] - min_time + 1

    # 最后一行标记为 NA：后续放所有患者出结果的时间
    enrollment_times[total_patients + 1, ] <- NA
    rownames(enrollment_times) <- c(1:total_patients, "finish_day")

    return(enrollment_times)
  }

  # 合并并排序患者入组时间表
  merge_and_sort_enrollment_times <- function(enrollment_times) {
    # 提取列数（组数）
    ngroup <- ncol(enrollment_times)

    # 初始化长表格
    long_table <- data.frame(Time = numeric(0), Group = numeric(0), PatientID = numeric(0))

    for (group in 1:ngroup) {
      # 提取组的时间列和患者编号，包括最后一行的 finish_day
      group_times <- enrollment_times[, group]
      total_patients <- nrow(enrollment_times) - 1
      patient_ids <- c(1:total_patients, n.cohortsize * cohortsize + 1)  # 将最后一个患者的 PatientID 设为总数 + 1

      # 创建当前组的表格
      group_table <- data.frame(
        Time = group_times,
        Group = group,
        PatientID = patient_ids
      )

      # 合并到长表格
      long_table <- rbind(long_table, group_table)
    }

    # 按时间排序，如果时间相同，则按 Group 和 PatientID 排序
    sorted_long_table <- long_table[order(long_table$Time, long_table$Group, long_table$PatientID, na.last = TRUE), ]

    # 删除行名
    rownames(sorted_long_table) <- NULL

    return(sorted_long_table)
  }

  # 定义函数：重新排序 enrolltimes_long_table 并重新分配 PatientID
  reorder_enrolltimes_table <- function(enrolltimes_long_table) {
    # 获取假患者的 ID
    fake_patient_id <- n.cohortsize * cohortsize + 1

    # 按 Time, Group, 原始 PatientID 排序，NA 排到最后
    enrolltimes_long_table <- enrolltimes_long_table[order(
      enrolltimes_long_table$Time,
      enrolltimes_long_table$Group,
      enrolltimes_long_table$PatientID,
      na.last = TRUE
    ), ]

    # 对每组内的 PatientID 重新编号
    enrolltimes_long_table$PatientID <- ave(
      seq_len(nrow(enrolltimes_long_table)),
      enrolltimes_long_table$Group,
      FUN = function(x) {
        ifelse(
          !is.na(enrolltimes_long_table$Time[x]) & enrolltimes_long_table$PatientID[x] != fake_patient_id,
          seq_along(x),  # 为非假患者重新分配连续编号
          fake_patient_id  # 保持假患者的 PatientID 不变
        )
      }
    )

    return(enrolltimes_long_table)
  }

  # 定义函数将长列转换回 4*n 的矩阵形式，带行列名称
  convert_to_named_original_format <- function(long_table, ngroup, total_patients) {
    # 初始化空矩阵
    original_matrix <- matrix(NA, nrow = total_patients + 1, ncol = ngroup)

    # 遍历每个组
    for (group in 1:ngroup) {
      # 筛选属于该组的患者
      group_data <- subset(long_table, Group == group)

      # 将患者的入组时间填入矩阵中
      original_matrix[1:nrow(group_data), group] <- group_data$Time
    }

    # 添加行名称
    rownames(original_matrix) <- c(1:total_patients, "finish_day")

    # 添加列名称
    colnames(original_matrix) <- paste("Group", 1:ngroup)

    # 返回结果
    return(original_matrix)
  }



  # # 函数使用示例
  # # 平均间隔设置：前两组为5天，后两组为10天
  # mean_intervals <- c(5, 5, 10, 10)
  #
  # # 生成入组时间表格
  # enrollment_times <- generate_enrollment_times(n.cohortsize, cohortsize, ngroup, mean_intervals)
  #
  # # 打印结果
  # print(enrollment_times)



  # 定义 mcmc_Delayed_BPCC_1st 和 mcmc_Delayed_BPCC_3rd
  mcmc.arms_Delayed_BPCC_1st <- function(dat, vbeta=vbeta, vgamma=vgamma, ngroup, type) {

    # MCMC with GIBBS Sampling
    #	含义：dat 是输入的数据集，每个组的毒性和疗效结果数据。
    #	类型：K个元素的列表，第k个列表元素是一个二维矩阵，对应第k个指标的患者数据
    # dat[[k]] = matrix(nrow = cohortsize, ncol = 3)
    # 第一列：毒性结果 Y_T，第二列：疗效结果 Y_E，第三列：剂量水平 u

    # 局部丢弃与全局 Burn-in 的对比
    #	局部丢弃：你提到的局部丢弃是每次在 Gibbs 更新某个参数时，通过多次采样选择“靠后”的一个值。这相当于对每个参数进行更深入的采样，类似于提高每次采样的质量，可能加快局部的收敛。
    #	这种做法在参数之间相互依赖性较强时，可能会有帮助，因为它可以使参数更新的结果更“稳定”。
    #	全局 Burn-in：全局 Burn-in 是从整个链的角度出发，丢弃初始非平稳的部分，确保从稳定的样本中进行推断。它通常已经足够，因为 MCMC 在经过足够多次迭代后，整个链都进入了平稳状态，采样值来自目标分布。

    N.post <- 1600
    N.burnin <- 800
    ngroup = ngroup   # ngroup: number of indications or trials
    # set initial values
    # parameter appears in model (2)

    # 11.20改进，rho初值改为0.2
    rho.hat <- 0.2

    # 这里的5个hat均为先验均值，在effdose函数里已经hat_beta_{gk}求平均，
    # 生成两个两维数组vbeta，vgamma（g=0,1,beta可换成gamma）
    beta0.hat <- vbeta[1] # 0.01, 0.3
    beta1.hat <- vbeta[2] # bingo
    gamma0.hat <- vgamma[1] # bingo
    gamma1.hat <- vgamma[2] # bingo

    # 11.20改进,alpha初值改为0.05
    alpha.hat <- 0.05

    # parameters appear in toxicity model (3)
    beta0 <- rep(beta0.hat, ngroup) # K维数组，各指标k均给出beta0.hat作为先验初值
    beta1 <- rep(beta1.hat, ngroup) # K维数组，各指标k均给出beta1.hat作为先验初值
    beta1_share <- log(beta1.hat)   # ～～？没看懂这里 shared among ngroup groups # 看懂了

    mu.beta0 <- beta0.hat # 第三层参数指定
    tau2.beta0 <- (4*mu.beta0)^2 # 第三层参数指定
    # tau2.beta0 <- (6*mu.beta0)^2 # 第三层参数指定 # 第4种情况更宽的先验
    mu.beta1 <- log(beta1.hat) # 第三层参数指定
    tau2.beta1 <- ((log(3)-mu.beta1)/qnorm(0.9))^2  # 第三层参数指定

    a.beta1 <- (2*mu.beta1)^2 # 看懂了，这个a.beta1应该就是delta.beta1
    # a.beta1 <- (3*mu.beta1)^2 # 看懂了，这个a.beta1应该就是delta.beta1 # 第四种情况 更宽的先验
    sigma2.beta1 <- 0.5*a.beta1 # sigma2.beta1的先验值取为均匀分布中间值
    # \sigma_{\beta_1}^2 in manuscript


    # parameters appear in efficacy model (4)
    rho <- rep(rho.hat, ngroup) # 第五种情况ind修改
    gamma0 <- rep(gamma0.hat, ngroup) # K维数组，各指标k均给出gamma0.hat作为先验初值
    gamma1 <- rep(gamma1.hat, ngroup) # K维数组，各指标k均给出gamma1.hat作为先验初值
    alpha <- rep(alpha.hat, ngroup)   # K维数组，各指标k均给出alpha.hat作为先验初值
    # gamma2 <- rep(0.5, ngroup)   # 12.4更新 K维数组，各指标k均给出0.5作为先验初值，这部分区域用于存放各indication的参数在当前轮次的值！！！
    gamma1_share <- log(gamma1.hat)  # ～～？同上，这里也没怎么看懂 shared among ngroup groups # 看懂了

    mu.gamma0 <- gamma0.hat # 第三层参数指定
    tau2.gamma0 <- (4*mu.gamma0)^2 # 第三层参数指定
    # tau2.gamma0 <- (6*mu.gamma0)^2 # 第4种情况更宽的先验
    mu.gamma1 <- log(gamma1.hat)  # 第三层参数指定

    if(type=='binary'){
      tau2.gamma1 <- ((log(3)-mu.gamma1)/qnorm(0.9))^2 # 第三层参数指定
      a.gamma1 <- 4*mu.gamma1^2 # 看懂了，这个a.gamma1应该就是delta.gamma1
      # a.gamma1 <- 9*mu.gamma1^2 # 看懂了，这个a.gamma1应该就是delta.gamma1 # 第四种情况 更宽的先验
      sigma2.gamma1 <- 0.5*a.gamma1 # sigma2.gamma1的先验值取为均匀分布中间值
      # \sigma_{\gamma_1}^2 in manuscript
    }else{
      tau2.gamma1 <- 4*mu.gamma1^2
      a.gamma1 <- 4*mu.gamma1^2
      sigma2.gamma1 <-  0.5*a.gamma1# \sigma_{\gamma_1}^2 in manuscript
    }

    #下面这两行注释不是ywt写的
    #sigma2.gamma1 <- 0.5*a.gamma1 # \sigma_{\gamma_1}^2 in manuscript
    #tau2.gamma1 <- 0.5*a.gamma1


    # 存储每个指标k里的,已知试验结果的患者数量n[k]
    # dat[[k]]是第k指标里的患者毒性/疗效/剂量（一共三列）数据，每行表示一个患者
    # 故n[k] = nrow(dat[[k]])
    n <- NULL
    for(k in 1:ngroup){
      n[k] = nrow(dat[[k]])
    }

    # N.post:给所有将要迭代的参数构造Npost*K的参数矩阵

    # alpha初值0.05
    # N.post:给所有将要迭代的参数构造Npost*K的参数矩阵
    rho_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta0_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta1_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta1_share_t <- rep(0, N.post)
    sigma2.beta1_t <- rep(0, N.post)

    gamma0_t <- matrix(0, nrow=N.post, ncol = ngroup)
    gamma1_t <- matrix(0, nrow=N.post, ncol = ngroup)
    alpha_t <- matrix(0, nrow=N.post, ncol = ngroup)
    # gamma2_t <- matrix(0, nrow=N.post, ncol = ngroup) #12.4更新，gamma2_t行为应该基本和alpha_t一致
    gamma1_share_t <- rep(0, N.post)
    sigma2.gamma1_t <- rep(0, N.post)


    # latent variables
    # Z_T[[k]] 是一个向量，长度与 dat[[k]]$y_T 的长度一致。
    Z_T = array(list(),c(1,1,ngroup))
    Z_T_t <- array(list(), dim = c(1, 1, ngroup))


    Z_E = array(list(),c(1,1,ngroup))
    Z_E_t <- array(list(), dim = c(1, 1, ngroup)) # Y_E for normal


    # 11.20 已做相应改动，根据观察数据y_T,y_E是否缺失、是否为1，赋值-0.05,+0.5,-0.5
    # 这里相当于给Z_T,Z_E赋初值了
    for(k in 1:ngroup){

      # 检测 dat[[k]]$y_T 是否包含 NA 值
      if (!is.data.frame(dat[[k]])) {
        # print(k)
        # print(dat[[k]])
        # print(dat[[1]])
        # print(dat[[2]])
        # print(dat[[3]])
        # print(dat[[4]])
        # print(class(dat[[k]]))
        stop(paste("dat[[", k, "]] is not a data frame. Please check the input data."))
      }

      if (any(is.na(dat[[k]]$y_T))) {
        warning(paste("第一次mcmc循环里，第", k, "组的 y_T 中有 NA 值，请仔细检查！"))
      }

      Z_T[[k]] = ifelse(is.na(dat[[k]]$y_T), -0.05, ifelse(dat[[k]]$y_T==1, 0.5, -0.5))
      if(type=='binary'){
        Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, ifelse(dat[[k]]$y_E==1, 0.5, -0.5))
      }else{
        Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, dat[[k]]$y_E)  # Y_E
      }
    }

    # ite是iteration迭代的缩写，代表当前的迭代次数
    # 为了减少资源浪费，可以将所有不需要在每次迭代中重新定义的函数
    # 移动到 for(ite in 1:N.post) 循环外，并通过参数传递适应每次迭代的不同输入。

    for(ite in 1:N.post){
      # generate the latent variables
      xi_cut <- c(-Inf, 0, Inf)
      zeta_cut <- c(-Inf, 0, Inf)

      # rho
      # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
      #	dat 是输入的数据集列表，第k个元素是第k组的毒性和疗效、剂量数据，nk行*3列的矩阵。
      #	dat 是 K 个元素的列表，第k个列表元素是一个二维矩阵，对应第k个指标的患者数据
      # 第一列：毒性结果 Y_T，第二列：疗效结果 Y_E，第三列：剂量水平 u

      # 用到了全局变量dat,Z_T,Z_E
      # log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
      #   Z_T.C <- Z_E.C <- NULL
      #   for(j in 1:ngroup){
      #     mean.Z_T <- beta0[j] + beta1[j]*dat[[j]][,3]
      #     # mean.Z_E <- gamma0[j]
      #     # + gamma1[j]*(ifelse(dat[[j]][,3]<=alpha[j], dat[[j]][,3], alpha[j])-gamma2[j]*(dat[[j]][,3]-alpha[j])*ifelse(dat[[j]][,3]>alpha[j], 1, 0))
      #     mean.Z_E <- gamma0[j] - gamma1[j] * (dat[[j]][,3] - alpha[j])^2
      #     Z_T.C <- c(Z_T.C, Z_T[[j]]-mean.Z_T)
      #     Z_E.C <- c(Z_E.C, Z_E[[j]]-mean.Z_E)
      #     # 妙啊，包含了所有sum_k nk的已治疗的患者数据；之后用sum求和即可
      #   }
      #   return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
      #   # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
      # }

      # # rho
      # logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
      # rho <- arms(rho, logden, function(x) ((x>-0.8)*(x<0.8)), 1)
      # rho_t[ite] <- rho
      # 关于arms第三个自变量位置的作用：
      # 所以，这个支持函数的作用就是筛选有效的采样值
      # 如果当前抽样的值不符合支持区间的条件（即返回FALSE）
      # 那么ARMS算法会重新进行抽样，直到找到一个满足条件的值。

      # 虽然ARMS从一个初始值开始（例如当前的rho）
      # 但这只是为了开始采样过程。在生成新的rho值时
      # ARMS会完全基于当前的对数密度函数
      # 因此rho的初始值并不会影响到最终采样的结果。
      # 采样的结果仅依赖于当前其他参数的取值和目标分布的形状。


      # rho: sigma_{12}^2 本次循环里将用var.z省略表示1-rho_t[ite]^2

      for(k in 1:ngroup){

        # rho_k   # 第五种情况ind修改
        log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
          Z_T.C <- Z_E.C <- NULL

          mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
          # mean.Z_E <- gamma0[k]
          # + gamma1[k]*(ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))
          mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][,3] - alpha[k])^2
          Z_T.C <- c(Z_T.C, Z_T[[k]]-mean.Z_T)
          Z_E.C <- c(Z_E.C, Z_E[[k]]-mean.Z_E)

          return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
          # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
        }

        logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
        rho[k] <- arms(rho[k], logden, function(x) ((x>-0.8)*(x<0.8)), 1)
        rho_t[ite,k] <- rho[k]

        var.z <- 1-rho[k]^2

        # 计算 mean.Z.T12 (这也是一个和dat[[k]][, 3]维度相同的向量)
        # mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho * (
        #   Z_E[[k]] - gamma0[k] - gamma1[k] * (ifelse(dat[[k]][, 3] <= alpha[k], dat[[k]][, 3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))
        # )
        mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho[k] * (
          Z_E[[k]] - gamma0[k] + gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
        )

        # 生成 Z_T[[k]]（逐元素生成）～～～～这里的步骤需要大幅度更改，或者说我重新写一个mcmc3，在那里再更改？
        # 因为 y_T 为NA的情况，说明已经是到了 第三步的那一次mcmc了～～～～～～～～～～～～～！！！
        Z_T[[k]] <- rtnorm(
          n = nrow(dat[[k]]),
          mean = mean.Z.T12,
          sd = sqrt(var.z),
          lower = ifelse(!is.na(dat[[k]]$y_T), xi_cut[dat[[k]]$y_T + 1], xi_cut[1]),
          upper = ifelse(!is.na(dat[[k]]$y_T), xi_cut[dat[[k]]$y_T + 2], xi_cut[3])
        )

        # 计算 mean.Z.E21
        # mean.Z.E21 <- gamma0[k] + gamma1[k] * (ifelse(dat[[k]][, 3] <= alpha[k], dat[[k]][, 3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0)) +
        #   rho * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][,3])
        mean.Z.E21 <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2 +
          rho[k] * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][, 3])

        # 根据 type 赋值 Z_E[[k]] 这个Z_E的赋值代码似乎已经没问题了，第一步和第三步mcmc都是这样赋值Z_E的
        if (type == 'binary') {
          Z_E[[k]] <- rtnorm(
            n = nrow(dat[[k]]),
            mean = mean.Z.E21,
            sd = sqrt(var.z),
            lower = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 1], xi_cut[1]),
            upper = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 2], xi_cut[3])
          )
        } else {
          Z_E[[k]] <- mean.Z.E21  # 非二元情况下直接赋值均值，这个暂时无所谓
        }

        # beta0k
        # 依然注意，下面提到的很多变量都是1维向量的形式，长度为nk
        # 是的，在 R 中，当你将一个向量和一个标量相加时
        # R会将该标量自动广播（broadcasting）到向量的每个元素上。
        # 也就是说，向量的每个元素都会加上该标量。
        # mean.Z_E <- gamma0[k]
        # + gamma1[k]*(ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))# bingo!
        mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
        mu.n.beta0 <- mean(Z_T[[k]]-beta1[k]*dat[[k]][,3]-rho[k]*(Z_E[[k]]-mean.Z_E)) # bingo! # 11.20 帅！
        sigma2.n.beta0 <- (1-rho[k]^2)/n[k] # bingo!

        # 检查并计算 beta0 的所有可能异常
        if (is.na(tau2.beta0)) {
          warning("mcmc1 tau2.beta0 is NA. It must be a positive number.")
        }
        if (tau2.beta0 <= 0) {
          warning("mcmc1 tau2.beta0 is invalid. It must be greater than 0.")
        }
        if (is.na(mu.n.beta0)) {
          warning("mcmc1 mu.n.beta0 is NA. Please check the input data or calculation.")
        }
        if (is.na(sigma2.n.beta0)) {
          warning("mcmc1 sigma2.n.beta0 is NA. It must be a valid positive number.")
        }
        if (sigma2.n.beta0 <= 0) {
          warning("mcmc1 sigma2.n.beta0 is invalid. It must be greater than 0.")
        }
        if (n[k] <= 0) {
          warning("mcmc1 n[k] is invalid. It must be a positive integer.")
        }

        mean.beta0 <- (tau2.beta0*mu.n.beta0+sigma2.n.beta0*mu.beta0)/(tau2.beta0+sigma2.n.beta0) # bingo!
        var.beta0 <- tau2.beta0*sigma2.n.beta0/(tau2.beta0+sigma2.n.beta0) # bingo!
        beta0[k]  <- rnorm(1, mean = mean.beta0, sd=sqrt(var.beta0)) # bingo!
        beta0_t[ite,k] <- beta0[k] # bingo!

        # 从这里上面最后一行和上面的两个for循环我们可以知道，ite是大循环，k是小循环；
        # 把第k组的每一个参数都迭代过之后，才会到第k+1组开启下一组的全参数迭代；
        # 每个第k组的参数都循环之后，然后才是ite进入下一轮

        # beta1k
        log.likelihood.beta1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.beta1, beta1_share, dat, Z_T, Z_E){
          mean.Z_T <- beta0 + beta1*dat[,3] # bingo！
          # mean.Z_E <- gamma0 + gamma1*(ifelse(dat[,3]<=alpha, dat[,3], alpha)-gamma2*(dat[,3]-alpha)*ifelse(dat[,3]>alpha, 1, 0)) # bingo！
          mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
          Z_T.C <- Z_T-mean.Z_T # bingo！
          Z_E.C <- Z_E-mean.Z_E # bingo！
          return(-1/2/(1-rho^2)*sum(Z_T.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.beta1*(log(beta1)-beta1_share)^2)
        }
        # 将logden的函数名修改为了logden_beta1k

        # 正确的代码
        if(TRUE){logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho[k], sigma2.beta1, beta1_share,
                                                                   dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)} #正确的！
        # # 去掉jacobi项之后的代码
        # if(TRUE){logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho, sigma2.beta1, beta1_share,
        #                                                            dat[[k]], Z_T[[k]], Z_E[[k]])}

        beta1[k]  <- arms(beta1[k], logden_beta1k, function(x) ((x>0)*(x<3)), 1)

        # beta1[k]记录当前迭代轮次和即将到来的迭代轮次的取值，beta1_t记录下所有迭代轮次的取值
        beta1_t[ite,k] <- beta1[k]


        # gamma0k
        # 依然注意，下面提到的很多变量都是1维向量的形式，长度为nk
        # 这一段代码检查完毕，检查无误

        mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
        # vstar <- ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0)
        vstar <- -(dat[[k]][,3]-alpha[k])^2
        mu.n.gamma0 <- mean(Z_E[[k]]-gamma1[k]*vstar-rho[k]*(Z_T[[k]]-mean.Z_T))
        sigma2.n.gamma0 <-(1-rho[k]^2)/n[k]

        # 检查并计算 gamma0 的所有可能异常
        if (is.na(tau2.gamma0)) {
          warning("mcmc1 tau2.gamma0 is NA. It must be a positive number.")
        }
        if (tau2.gamma0 <= 0) {
          warning("mcmc1 tau2.gamma0 is invalid. It must be greater than 0.")
        }
        if (is.na(mu.n.gamma0)) {
          warning("mcmc1 mu.n.gamma0 is NA. Please check the input data or calculation.")
        }
        if (is.na(sigma2.n.gamma0)) {
          warning("mcmc1 sigma2.n.gamma0 is NA. It must be a valid positive number.")
        }
        if (sigma2.n.gamma0 <= 0) {
          warning("mcmc1 sigma2.n.gamma0 is invalid. It must be greater than 0.")
        }
        if (n[k] <= 0) {
          warning("mcmc1 n[k] is invalid. It must be a positive integer.")
        }

        mean.gamma0 <- (tau2.gamma0*mu.n.gamma0+sigma2.n.gamma0*mu.gamma0)/(tau2.gamma0+sigma2.n.gamma0)
        # ～～下面这一行代码应该是多写了，现注释掉
        # var.n.gamma0 <- tau2.gamma0*sigma2.n.gamma0/(tau2.gamma0+sigma2.n.gamma0)
        var.gamma0 <- tau2.gamma0*sigma2.n.gamma0 /(tau2.gamma0+sigma2.n.gamma0)
        gamma0[k]  <- rnorm(1, mean = mean.gamma0, sd=sqrt(var.gamma0))
        gamma0_t[ite,k] <- gamma0[k]


        # gamma1k
        log.likelihood.gamma1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.gamma1, gamma1_share, dat, Z_T, Z_E){
          mean.Z_T <- beta0 + beta1*dat[,3]
          # mean.Z_E <- gamma0 + gamma1*(ifelse(dat[,3]<=alpha, dat[,3], alpha)-gamma2*(dat[,3]-alpha)*ifelse(dat[,3]>alpha, 1, 0))
          mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
          Z_T.C <- Z_T-mean.Z_T
          Z_E.C <- Z_E-mean.Z_E
          return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.gamma1*(log(gamma1)-gamma1_share)^2)
        }
        #正确的代码
        if(TRUE){logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho[k], sigma2.gamma1, gamma1_share,
                                                                     dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)}

        # #去掉jacobi项之后的代码
        # if(TRUE){logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho, sigma2.gamma1, gamma1_share,
        #                                                              dat[[k]], Z_T[[k]], Z_E[[k]])}

        # 限制斜率大小，gamma1k大于0、但不超过5
        # 12.6把上面的上界4更新为5
        gamma1.Star <- ifelse(type=='binary', 5, 5)
        gamma1[k]  <- arms(gamma1[k], logden_gamma1k, function(x) ((x>0)*(x<gamma1.Star)), 1)
        gamma1_t[ite,k] <- gamma1[k]


        # alphak
        # 这一段代码已检查完毕

        log.likelihood.alpha <- function(alphak, rho, beta0k, beta1k, gamma0k, gamma1k, Z_Tk, Z_Ek, datk) {
          u = datk[,3] # 第k组标准化后的剂量
          mean.Z_T <- beta0k + beta1k*u
          # mean.Z_E <- gamma0k + gamma1k*(ifelse(u<=alphak, u, alphak)-gamma2k*(u-alphak)*(ifelse(u>alphak,1,0)))
          mean.Z_E <- gamma0k - gamma1k * (u - alphak)^2
          Z_T.C <- Z_Tk-mean.Z_T
          Z_E.C <- Z_Ek-mean.Z_E
          return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C))
        }
        logden_alphak <- function(x) log.likelihood.alpha(x, rho[k], beta0[k], beta1[k], gamma0[k], gamma1[k], Z_T[[k]], Z_E[[k]], dat[[k]])

        # 12.4改为-0.25~0.5之间；剂量被标准化为-0.5～0.5之间，所以感觉应该不是-1～1，这里FALSE区是原代码
        if(FALSE){alpha[k] <- arms(alpha[k],logden, function(x) ((x>-1)*(x<1)), 1)}
        if(TRUE){alpha[k] <- arms(alpha[k],logden_alphak, function(x) ((x>-0.5)*(x<0.5)), 1)}
        alpha_t[ite,k] <- alpha[k]

      }


      # 注意：在每一轮下，在所有第k=1～K组的表观信息层都迭代完之后，再来迭代共享信息层！

      # beta1_share
      # 这一段代码已检查完毕
      mu.n.beta1_share <- mean(log(beta1))
      sigma2.n.beta1_share <- sigma2.beta1/ngroup

      mean.beta1_share <- (tau2.beta1*mu.n.beta1_share+sigma2.n.beta1_share*mu.beta1)/(tau2.beta1+sigma2.n.beta1_share)
      var.beta1_share <- (tau2.beta1*sigma2.n.beta1_share)/(tau2.beta1+sigma2.n.beta1_share)

      beta1_share <- rnorm(1, mean=mean.beta1_share, sd = sqrt(var.beta1_share))
      beta1_share_t[ite] <- beta1_share


      # sigma2.beta1
      # 这一段代码已检查完毕
      logden_sigma2.beta1 <- function(x) -ngroup/2*log(x)- sum((log(beta1)-beta1_share)^2)/(2*x)
      sigma2.beta1 <- arms(sigma2.beta1, logden_sigma2.beta1, function(x) (x>0)*(x<a.beta1), 1)
      sigma2.beta1_t[ite] <- sigma2.beta1


      # gamma1_share
      # 这一段代码已检查完毕
      mu.n.gamma1_share <- mean(log(gamma1))
      sigma2.n.gamma1_share <- sigma2.gamma1/ngroup

      mean.gamma1_share <- (tau2.gamma1*mu.n.gamma1_share + sigma2.n.gamma1_share*mu.gamma1)/(tau2.gamma1+sigma2.n.gamma1_share)
      var.gamma1_share <- (tau2.gamma1*sigma2.n.gamma1_share)/(tau2.gamma1+sigma2.n.gamma1_share)
      gamma1_share <- rnorm(1, mean = mean.gamma1_share, sd=sqrt(var.gamma1_share))
      gamma1_share_t[ite] <- gamma1_share


      # sigma2.gamma1
      # 这一段代码已检查完毕
      logden_sigma2.gamma1 <- function(x) -ngroup/2*log(x) - sum((log(gamma1)-gamma1_share)^2)/(2*x)
      sigma2.gamma1 <- arms(sigma2.gamma1, logden_sigma2.gamma1, function(x) (x>0)*(x<a.gamma1), 1)
      sigma2.gamma1_t[ite] <- sigma2.gamma1

      #迭代过程至此完成完整的一轮，一共循环完成N.post轮
    }

    #取一个序列（index）来表明，只取烧录期之后的稳定样本
    ind <- seq((N.burnin+1),N.post)

    # 下面这行注释不是ywt写的
    # matrix.post <- cbind()[ind, ]
    matrix.uti <- array(list(),dim=c(1,1,ngroup))

    for(k in 1:ngroup){
      matrix.uti[[k]] = cbind(beta0_t[,k], beta1_t[,k], gamma0_t[,k],
                              gamma1_t[,k], alpha_t[,k], rho_t[,k])[ind, ]

      # matrix.uti 是一个列表，其中的每一个元素（matrix.uti[[k]]）是对应某个指标的参数矩阵。
      # matrix.uti[[k]]记录了k指标下，所有**表观信息层参数**的所有迭代轮次的取值
      # 以beta0_t为例子，这是一个Npost行*K列的数组，记录了每个指标k下，在Npost次迭代过程中的所有值
      # 所以这里等号右边是选取了第k个指标下的全体参数的Npost次的迭代值
    }
    return(matrix.uti)
  }

  mcmc.arms_Delayed_BPCC_3rd <- function(dat, vbeta=vbeta, vgamma=vgamma, ngroup, type, c_hat, eta_0_hat, eta_1_hat) {

    # MCMC with GIBBS Sampling
    #	含义：dat 是输入的数据集，每个组的毒性和疗效结果数据。
    #	类型：K个元素的列表，第k个列表元素是一个二维矩阵，对应第k个指标的患者数据
    # dat[[k]] = matrix(nrow = cohortsize, ncol = 3)
    # 第一列：毒性结果 Y_T，第二列：疗效结果 Y_E，第三列：剂量水平 u

    # 局部丢弃与全局 Burn-in 的对比
    #	局部丢弃：你提到的局部丢弃是每次在 Gibbs 更新某个参数时，通过多次采样选择“靠后”的一个值。这相当于对每个参数进行更深入的采样，类似于提高每次采样的质量，可能加快局部的收敛。
    #	这种做法在参数之间相互依赖性较强时，可能会有帮助，因为它可以使参数更新的结果更“稳定”。
    #	全局 Burn-in：全局 Burn-in 是从整个链的角度出发，丢弃初始非平稳的部分，确保从稳定的样本中进行推断。它通常已经足够，因为 MCMC 在经过足够多次迭代后，整个链都进入了平稳状态，采样值来自目标分布。

    N.post <- 3000
    N.burnin <- 1000
    ngroup = ngroup   # ngroup: number of indications or trials
    # set initial values
    # parameter appears in model (2)

    # 11.20改进，rho初值改为0.2
    rho.hat <- 0.2

    # 这里的5个hat均为先验均值，在effdose函数里已经hat_beta_{gk}求平均，
    # 生成两个两维数组vbeta，vgamma（g=0,1,beta可换成gamma）
    beta0.hat <- vbeta[1] # 0.01, 0.3
    beta1.hat <- vbeta[2] # bingo
    gamma0.hat <- vgamma[1] # bingo
    gamma1.hat <- vgamma[2] # bingo

    # 11.20改进,alpha初值改为0.05
    alpha.hat <- 0.05

    # parameters appear in toxicity model (3)
    beta0 <- rep(beta0.hat, ngroup) # K维数组，各指标k均给出beta0.hat作为先验初值
    beta1 <- rep(beta1.hat, ngroup) # K维数组，各指标k均给出beta1.hat作为先验初值
    beta1_share <- log(beta1.hat)   # ～～？没看懂这里 shared among ngroup groups # 看懂了

    mu.beta0 <- beta0.hat # 第三层参数指定
    tau2.beta0 <- (4*mu.beta0)^2 # 第三层参数指定
    # tau2.beta0 <- (6*mu.beta0)^2 # 第三层参数指定 # 第4种情况更宽的先验
    mu.beta1 <- log(beta1.hat) # 第三层参数指定
    tau2.beta1 <- ((log(3)-mu.beta1)/qnorm(0.9))^2  # 第三层参数指定

    a.beta1 <- (2*mu.beta1)^2 # 看懂了，这个a.beta1应该就是delta.beta1
    # a.beta1 <- (3*mu.beta1)^2 # 看懂了，这个a.beta1应该就是delta.beta1 # 第四种情况 更宽的先验
    sigma2.beta1 <- 0.5*a.beta1 # sigma2.beta1的先验值取为均匀分布中间值
    # \sigma_{\beta_1}^2 in manuscript


    # parameters appear in efficacy model (4)
    rho <- rep(rho.hat, ngroup) # 第五种情况ind修改
    gamma0 <- rep(gamma0.hat, ngroup) # K维数组，各指标k均给出gamma0.hat作为先验初值
    gamma1 <- rep(gamma1.hat, ngroup) # K维数组，各指标k均给出gamma1.hat作为先验初值
    alpha <- rep(alpha.hat, ngroup)   # K维数组，各指标k均给出alpha.hat作为先验初值
    # gamma2 <- rep(0.5, ngroup)   # 12.4更新 K维数组，各指标k均给出0.5作为先验初值，这部分区域用于存放各indication的参数在当前轮次的值！！！
    gamma1_share <- log(gamma1.hat)  # ～～？同上，这里也没怎么看懂 shared among ngroup groups # 看懂了

    mu.gamma0 <- gamma0.hat # 第三层参数指定
    tau2.gamma0 <- (4*mu.gamma0)^2 # 第三层参数指定
    # tau2.gamma0 <- (6*mu.gamma0)^2 # 第4种情况更宽的先验
    mu.gamma1 <- log(gamma1.hat)  # 第三层参数指定

    if(type=='binary'){
      tau2.gamma1 <- ((log(3)-mu.gamma1)/qnorm(0.9))^2 # 第三层参数指定
      a.gamma1 <- 4*mu.gamma1^2 # 看懂了，这个a.gamma1应该就是delta.gamma1
      # a.gamma1 <- 9*mu.gamma1^2 # 看懂了，这个a.gamma1应该就是delta.gamma1 # 第四种情况 更宽的先验
      sigma2.gamma1 <- 0.5*a.gamma1 # sigma2.gamma1的先验值取为均匀分布中间值
      # \sigma_{\gamma_1}^2 in manuscript
    }else{
      tau2.gamma1 <- 4*mu.gamma1^2
      a.gamma1 <- 4*mu.gamma1^2
      sigma2.gamma1 <-  0.5*a.gamma1# \sigma_{\gamma_1}^2 in manuscript
    }

    #下面这两行注释不是ywt写的
    #sigma2.gamma1 <- 0.5*a.gamma1 # \sigma_{\gamma_1}^2 in manuscript
    #tau2.gamma1 <- 0.5*a.gamma1


    # 存储每个指标k里的,已知试验结果的患者数量n[k]
    # dat[[k]]是第k指标里的患者毒性/疗效/剂量（一共三列）数据，每行表示一个患者
    # 故n[k] = nrow(dat[[k]])
    n <- NULL
    for(k in 1:ngroup){
      n[k] = nrow(dat[[k]])
    }

    # N.post:给所有将要迭代的参数构造Npost*K的参数矩阵
    # alpha初值0.05
    # N.post:给所有将要迭代的参数构造Npost*K的参数矩阵
    rho_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta0_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta1_t <- matrix(0, nrow=N.post, ncol = ngroup)
    beta1_share_t <- rep(0, N.post)
    sigma2.beta1_t <- rep(0, N.post)

    gamma0_t <- matrix(0, nrow=N.post, ncol = ngroup)
    gamma1_t <- matrix(0, nrow=N.post, ncol = ngroup)
    alpha_t <- matrix(0, nrow=N.post, ncol = ngroup)
    # gamma2_t <- matrix(0, nrow=N.post, ncol = ngroup) #12.4更新，gamma2_t行为应该基本和alpha_t一致
    gamma1_share_t <- rep(0, N.post)
    sigma2.gamma1_t <- rep(0, N.post)


    # latent variables
    # ～～～？？？这里四行没怎么看懂？？？
    # Z_T[[k]] 是一个向量，长度与 dat[[k]]$y_T 的长度一致。
    Z_T = array(list(),c(1,1,ngroup))
    Z_T_t <- array(list(), dim = c(1, 1, ngroup))


    Z_E = array(list(),c(1,1,ngroup))
    Z_E_t <- array(list(), dim = c(1, 1, ngroup)) # Y_E for normal


    # 11.20 已做相应改动，根据观察数据y_T,y_E是否缺失、是否为1，赋值-0.05,+0.5,-0.5
    # 这里相当于给Z_T,Z_E赋初值了
    for(k in 1:ngroup){

      # 第三轮不用再检测 dat[[k]]$y_T 是否包含 NA 值
      # if (any(is.na(dat[[k]]$y_T))) {
      #   warning(paste("第一次mcmc循环里，第", k, "组的 y_T 中有 NA 值，请仔细检查！"))
      # }

      Z_T[[k]] = ifelse(is.na(dat[[k]]$y_T), -0.05, ifelse(dat[[k]]$y_T==1, 0.5, -0.5))
      if(type=='binary'){
        Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, ifelse(dat[[k]]$y_E==1, 0.5, -0.5))
      }else{
        Z_E[[k]] = ifelse(is.na(dat[[k]]$y_E), -0.05, dat[[k]]$y_E)  # Y_E
      }
    }

    # ite是iteration迭代的缩写，代表当前的迭代次数
    # 为了减少资源浪费，可以将所有不需要在每次迭代中重新定义的函数
    # 移动到 for(ite in 1:N.post) 循环外，并通过参数传递适应每次迭代的不同输入。

    for(ite in 1:N.post){
      # generate the latent variables
      xi_cut <- c(-Inf, 0, Inf)
      zeta_cut <- c(-Inf, 0, Inf)

      # # rho
      # # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
      # #	dat 是输入的数据集列表，第k个元素是第k组的毒性和疗效、剂量数据，nk行*3列的矩阵。
      # #	dat 是 K 个元素的列表，第k个列表元素是一个二维矩阵，对应第k个指标的患者数据
      # # 第一列：毒性结果 Y_T，第二列：疗效结果 Y_E，第三列：剂量水平 u
      #
      # # 用到了全局变量dat,Z_T,Z_E
      # log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
      #   Z_T.C <- Z_E.C <- NULL
      #   for(j in 1:ngroup){
      #     mean.Z_T <- beta0[j] + beta1[j]*dat[[j]][,3]
      #     # mean.Z_E <- gamma0[j]
      #     # + gamma1[j]*(ifelse(dat[[j]][,3]<=alpha[j], dat[[j]][,3], alpha[j])-gamma2[j]*(dat[[j]][,3]-alpha[j])*ifelse(dat[[j]][,3]>alpha[j], 1, 0))
      #     mean.Z_E <- gamma0[j] - gamma1[j] * (dat[[j]][,3] - alpha[j])^2
      #     Z_T.C <- c(Z_T.C, Z_T[[j]]-mean.Z_T)
      #     Z_E.C <- c(Z_E.C, Z_E[[j]]-mean.Z_E)
      #     # 妙啊，包含了所有sum_k nk的已治疗的患者数据；之后用sum求和即可
      #   }
      #   return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
      #   # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
      # }

      # # rho
      # logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
      # rho <- arms(rho, logden, function(x) ((x>-0.8)*(x<0.8)), 1)
      # rho_t[ite] <- rho
      # 关于arms第三个自变量位置的作用：
      # 所以，这个支持函数的作用就是筛选有效的采样值
      # 如果当前抽样的值不符合支持区间的条件（即返回FALSE）
      # 那么ARMS算法会重新进行抽样，直到找到一个满足条件的值。

      # 虽然ARMS从一个初始值开始（例如当前的rho）
      # 但这只是为了开始采样过程。在生成新的rho值时
      # ARMS会完全基于当前的对数密度函数
      # 因此rho的初始值并不会影响到最终采样的结果。
      # 采样的结果仅依赖于当前其他参数的取值和目标分布的形状。
      # rho: sigma_{12}^2 本次循环里将用var.z省略表示1-rho_t[ite]^2

      for(k in 1:ngroup){

        # rho_k   # 第五种情况ind修改
        log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
          Z_T.C <- Z_E.C <- NULL

          mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
          # mean.Z_E <- gamma0[k]
          # + gamma1[k]*(ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))
          mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][,3] - alpha[k])^2
          Z_T.C <- c(Z_T.C, Z_T[[k]]-mean.Z_T)
          Z_E.C <- c(Z_E.C, Z_E[[k]]-mean.Z_E)

          return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
          # 这里就是所有基于已知全部参数数据的对数似然函数log L(zT,zE,\theta)
        }

        logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
        rho[k] <- arms(rho[k], logden, function(x) ((x>-0.8)*(x<0.8)), 1)
        rho_t[ite,k] <- rho[k]

        var.z <- 1-rho[k]^2
        # 这里的dat[[k]][,3]是一个nk长度的列向量数据，表明第k组所有患者剂量数据
        # 同理这里提到的许多变量都是nk长度的
        # ～～?!出大问题，我小导写的代码疑似把条件均值求错了，+号写成了-号

        # if(FALSE){ #～～？这一块代码是我小导写的，疑似条件均值求错了
        #   mean.Z.E <- gamma0[k] + gamma1[k]*ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k]) - rho*(Z_T[[k]]-beta0[k] - beta1[k]*dat[[k]][,3])
        #   if(type=='binary'){
        #    Z_E[[k]] <- rtnorm(n[k], mean=mean.Z.E, sd = sqrt(var.z), lower = xi_cut[dat[[k]][,2]+1], upper = xi_cut[dat[[k]][,2]+2])
        #   }
        #   #
        #   mean.Z.T <- beta0[k] + beta1[k]*dat[[k]][,3] - rho*(Z_E[[k]]-gamma0[k] - gamma1[k]*ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k]))
        #   Z_T[[k]] <- rtnorm(n[k], mean=mean.Z.T, sd = sqrt(var.z), lower = xi_cut[dat[[k]][,1]+1], upper = xi_cut[dat[[k]][,1]+2])
        # }

        # 根据y_T,y_E和已有的前一轮参数，迭代赋值Z_T,Z_E
        # 11.20 更改，感觉Z_T最好比Z_E先赋值，因为主要是Z_T的结果去影响Z_E,后者的已知信息更少
        # 11.20 更改，重新书写逻辑，加上了缺失数据y_T,y_E时，Z_T,Z_E的分支赋值

        # 计算 mean.Z.T12 (这也是一个和dat[[k]][, 3]维度相同的向量)
        # mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho * (
        #   Z_E[[k]] - gamma0[k] - gamma1[k] * (ifelse(dat[[k]][, 3] <= alpha[k], dat[[k]][, 3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))
        # )
        mean.Z.T12 <- beta0[k] + beta1[k] * dat[[k]][, 3] + rho[k] * (
          Z_E[[k]] - gamma0[k] + gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
        )

        ########## 第三步的mcmc中最关键的一步！！！！！！
        # 生成 Z_T[[k]]（使用sapply逐元素生成) 已修改为适应第三步mcmc的生成方式！！！
        # 第三步的 MCMC 中 Z_T 的生成逻辑
        if (anyNA(dat[[k]][, 3])){
          print(k)
          print(nrow(dat[[k]]))
          print((dat[[k]]))
          stop("NA detected in dat[[k]][, 3]")
        }

        Z_T[[k]] <- sapply(1:nrow(dat[[k]]), function(i) {
          if (!is.na(dat[[k]]$y_T[i])) {
            # 情况 (y_T=1/0)：直接根据 y_T 的值生成
            rtnorm(
              n = 1,
              mean = mean.Z.T12[i],
              sd = sqrt(var.z),
              lower = xi_cut[dat[[k]]$y_T[i] + 1],
              upper = xi_cut[dat[[k]]$y_T[i] + 2]
            )
          } else {
            # 情况 (y_T = NA)：需要计算 p_{k,i}
            F_val <- pnorm(beta0[k] + beta1[k] * dat[[k]][i, 3])  # 计算 F 值
            W_T_val <- W_T(dat[[k]]$t_follow[i] / U_T, dat[[k]][i, 3], c_hat, eta_0_hat, eta_1_hat)  # 计算 W^T
            p_ki <- F_val * (1 - W_T_val) / (1 - F_val * W_T_val)  # 计算 p_{k,i}
            # 采样随机数 u 并根据 u 决定截断范围
            u <- runif(1)
            if (u <= p_ki) {
              rtnorm(n = 1, mean = mean.Z.T12[i], sd = sqrt(var.z), lower = 0, upper = Inf)
            } else {
              rtnorm(n = 1, mean = mean.Z.T12[i], sd = sqrt(var.z), lower = -Inf, upper = 0)
            }
          }
        })

        # 计算 mean.Z.E21
        # mean.Z.E21 <- gamma0[k] + gamma1[k] * (ifelse(dat[[k]][, 3] <= alpha[k], dat[[k]][, 3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0)) +
        #   rho * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][,3])
        mean.Z.E21 <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2 +
          rho[k] * (Z_T[[k]] - beta0[k] - beta1[k] * dat[[k]][, 3])

        # 根据 type 赋值 Z_E[[k]] 这个Z_E的赋值代码似乎已经没问题了，第一步和第三步mcmc都是这样赋值Z_E的
        if (type == 'binary') {
          Z_E[[k]] <- rtnorm(
            n = nrow(dat[[k]]),
            mean = mean.Z.E21,
            sd = sqrt(var.z),
            lower = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 1], xi_cut[1]),
            upper = ifelse(!is.na(dat[[k]]$y_E), xi_cut[dat[[k]]$y_E + 2], xi_cut[3])
          )
        } else {
          Z_E[[k]] <- mean.Z.E21  # 非二元情况下直接赋值均值，这个暂时无所谓
        }


        # beta0k
        # 依然注意，下面提到的很多变量都是1维向量的形式，长度为nk
        # 是的，在 R 中，当你将一个向量和一个标量相加时
        # R会将该标量自动广播（broadcasting）到向量的每个元素上。
        # 也就是说，向量的每个元素都会加上该标量。
        # mean.Z_E <- gamma0[k]
        # + gamma1[k]*(ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0))# bingo!
        mean.Z_E <- gamma0[k] - gamma1[k] * (dat[[k]][, 3] - alpha[k])^2
        mu.n.beta0 <- mean(Z_T[[k]]-beta1[k]*dat[[k]][,3]-rho[k]*(Z_E[[k]]-mean.Z_E)) # bingo! # 11.20 帅！
        sigma2.n.beta0 <- (1-rho[k]^2)/n[k] # bingo!

        # 检查并计算 beta0 的所有可能异常
        if (is.na(tau2.beta0)) {
          warning("mcmc1 tau2.beta0 is NA. It must be a positive number.")
        }
        if (tau2.beta0 <= 0) {
          warning("mcmc1 tau2.beta0 is invalid. It must be greater than 0.")
        }
        if (is.na(mu.n.beta0)) {
          warning("mcmc1 mu.n.beta0 is NA. Please check the input data or calculation.")
        }
        if (is.na(sigma2.n.beta0)) {
          warning("mcmc1 sigma2.n.beta0 is NA. It must be a valid positive number.")
        }
        if (sigma2.n.beta0 <= 0) {
          warning("mcmc1 sigma2.n.beta0 is invalid. It must be greater than 0.")
        }
        if (n[k] <= 0) {
          warning("mcmc1 n[k] is invalid. It must be a positive integer.")
        }
        # var.beta0 <- tau2.beta0 * sigma2.n.beta0 / (tau2.beta0 + sigma2.n.beta0)
        # mean.beta0 <- (tau2.beta0 * mu.n.beta0 + sigma2.n.beta0 * mu.beta0) / (tau2.beta0 + sigma2.n.beta0)
        # beta0[k] <- rnorm(1, mean = mean.beta0, sd = sqrt(var.beta0))
        # 上面这三行需要注释掉，莫名其妙重复了
        mean.beta0 <- (tau2.beta0*mu.n.beta0+sigma2.n.beta0*mu.beta0)/(tau2.beta0+sigma2.n.beta0) # bingo!
        var.beta0 <- tau2.beta0*sigma2.n.beta0/(tau2.beta0+sigma2.n.beta0) # bingo!
        beta0[k]  <- rnorm(1, mean = mean.beta0, sd=sqrt(var.beta0)) # bingo!
        beta0_t[ite,k] <- beta0[k] # bingo!

        # 从这里上面最后一行和上面的两个for循环我们可以知道，ite是大循环，k是小循环；
        # 把第k组的每一个参数都迭代过之后，才会到第k+1组开启下一组的全参数迭代；
        # 每个第k组的参数都循环之后，然后才是ite进入下一轮

        # beta1k
        log.likelihood.beta1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.beta1, beta1_share, dat, Z_T, Z_E){
          mean.Z_T <- beta0 + beta1*dat[,3] # bingo！
          # mean.Z_E <- gamma0 + gamma1*(ifelse(dat[,3]<=alpha, dat[,3], alpha)-gamma2*(dat[,3]-alpha)*ifelse(dat[,3]>alpha, 1, 0)) # bingo！
          mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
          Z_T.C <- Z_T-mean.Z_T # bingo！
          Z_E.C <- Z_E-mean.Z_E # bingo！
          return(-1/2/(1-rho^2)*sum(Z_T.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.beta1*(log(beta1)-beta1_share)^2)
        }
        # 将logden的函数名修改为了logden_beta1k

        # 正确的代码
        if(TRUE){logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho[k], sigma2.beta1, beta1_share,
                                                                   dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)} #正确的！
        # # 去掉jacobi项之后的代码
        # if(TRUE){logden_beta1k <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho, sigma2.beta1, beta1_share,
        #                                                            dat[[k]], Z_T[[k]], Z_E[[k]])}

        beta1[k]  <- arms(beta1[k], logden_beta1k, function(x) ((x>0)*(x<3)), 1)

        # beta1[k]记录当前迭代轮次和即将到来的迭代轮次的取值，beta1_t记录下所有迭代轮次的取值
        beta1_t[ite,k] <- beta1[k]


        # gamma0k
        # 依然注意，下面提到的很多变量都是1维向量的形式，长度为nk
        # 这一段代码检查完毕，检查无误

        mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
        # vstar <- ifelse(dat[[k]][,3]<=alpha[k], dat[[k]][,3], alpha[k])-gamma2[k]*(dat[[k]][,3]-alpha[k])*ifelse(dat[[k]][,3]>alpha[k], 1, 0)
        vstar <- -(dat[[k]][,3]-alpha[k])^2
        mu.n.gamma0 <- mean(Z_E[[k]]-gamma1[k]*vstar-rho[k]*(Z_T[[k]]-mean.Z_T))
        sigma2.n.gamma0 <-(1-rho[k]^2)/n[k]

        # 检查并计算 gamma0 的所有可能异常
        if (is.na(tau2.gamma0)) {
          warning("mcmc1 tau2.gamma0 is NA. It must be a positive number.")
        }
        if (tau2.gamma0 <= 0) {
          warning("mcmc1 tau2.gamma0 is invalid. It must be greater than 0.")
        }
        if (is.na(mu.n.gamma0)) {
          warning("mcmc1 mu.n.gamma0 is NA. Please check the input data or calculation.")
        }
        if (is.na(sigma2.n.gamma0)) {
          warning("mcmc1 sigma2.n.gamma0 is NA. It must be a valid positive number.")
        }
        if (sigma2.n.gamma0 <= 0) {
          warning("mcmc1 sigma2.n.gamma0 is invalid. It must be greater than 0.")
        }
        if (n[k] <= 0) {
          warning("mcmc1 n[k] is invalid. It must be a positive integer.")
        }

        mean.gamma0 <- (tau2.gamma0*mu.n.gamma0+sigma2.n.gamma0*mu.gamma0)/(tau2.gamma0+sigma2.n.gamma0)
        # ～～下面这一行代码应该是多写了，现注释掉
        # var.n.gamma0 <- tau2.gamma0*sigma2.n.gamma0/(tau2.gamma0+sigma2.n.gamma0)
        var.gamma0 <- tau2.gamma0*sigma2.n.gamma0 /(tau2.gamma0+sigma2.n.gamma0)
        gamma0[k]  <- rnorm(1, mean = mean.gamma0, sd=sqrt(var.gamma0))
        gamma0_t[ite,k] <- gamma0[k]


        # gamma1k
        log.likelihood.gamma1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.gamma1, gamma1_share, dat, Z_T, Z_E){
          mean.Z_T <- beta0 + beta1*dat[,3]
          # mean.Z_E <- gamma0 + gamma1*(ifelse(dat[,3]<=alpha, dat[,3], alpha)-gamma2*(dat[,3]-alpha)*ifelse(dat[,3]>alpha, 1, 0))
          mean.Z_E <- gamma0 - gamma1 * (dat[, 3] - alpha)^2
          Z_T.C <- Z_T-mean.Z_T
          Z_E.C <- Z_E-mean.Z_E
          return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.gamma1*(log(gamma1)-gamma1_share)^2)
        }
        #正确的代码
        if(TRUE){logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho[k], sigma2.gamma1, gamma1_share,
                                                                     dat[[k]], Z_T[[k]], Z_E[[k]])-log(x)}

        # #去掉jacobi项之后的代码
        # if(TRUE){logden_gamma1k <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho, sigma2.gamma1, gamma1_share,
        #                                                              dat[[k]], Z_T[[k]], Z_E[[k]])}

        # 限制斜率大小，gamma1k大于0、但不超过5
        # 12.6把上面的上界4更新为5
        gamma1.Star <- ifelse(type=='binary', 5, 5)
        gamma1[k]  <- arms(gamma1[k], logden_gamma1k, function(x) ((x>0)*(x<gamma1.Star)), 1)
        gamma1_t[ite,k] <- gamma1[k]


        # alphak
        # 这一段代码已检查完毕

        log.likelihood.alpha <- function(alphak, rho, beta0k, beta1k, gamma0k, gamma1k, Z_Tk, Z_Ek, datk) {
          u = datk[,3] # 第k组标准化后的剂量
          mean.Z_T <- beta0k + beta1k*u
          # mean.Z_E <- gamma0k + gamma1k*(ifelse(u<=alphak, u, alphak)-gamma2k*(u-alphak)*(ifelse(u>alphak,1,0)))
          mean.Z_E <- gamma0k - gamma1k * (u - alphak)^2
          Z_T.C <- Z_Tk-mean.Z_T
          Z_E.C <- Z_Ek-mean.Z_E
          return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C))
        }
        logden_alphak <- function(x) log.likelihood.alpha(x, rho[k], beta0[k], beta1[k], gamma0[k], gamma1[k], Z_T[[k]], Z_E[[k]], dat[[k]])

        # 12.4改为-0.25~0.5之间；剂量被标准化为-0.5～0.5之间，所以感觉应该不是-1～1，这里FALSE区是原代码
        if(FALSE){alpha[k] <- arms(alpha[k],logden, function(x) ((x>-1)*(x<1)), 1)}
        if(TRUE){alpha[k] <- arms(alpha[k],logden_alphak, function(x) ((x>-0.5)*(x<0.5)), 1)}
        alpha_t[ite,k] <- alpha[k]

      }


      # 注意：在每一轮下，在所有第k=1～K组的表观信息层都迭代完之后，再来迭代共享信息层！

      # beta1_share
      # 这一段代码已检查完毕
      mu.n.beta1_share <- mean(log(beta1))
      sigma2.n.beta1_share <- sigma2.beta1/ngroup

      mean.beta1_share <- (tau2.beta1*mu.n.beta1_share+sigma2.n.beta1_share*mu.beta1)/(tau2.beta1+sigma2.n.beta1_share)
      var.beta1_share <- (tau2.beta1*sigma2.n.beta1_share)/(tau2.beta1+sigma2.n.beta1_share)

      beta1_share <- rnorm(1, mean=mean.beta1_share, sd = sqrt(var.beta1_share))
      beta1_share_t[ite] <- beta1_share


      # sigma2.beta1
      # 这一段代码已检查完毕
      logden_sigma2.beta1 <- function(x) -ngroup/2*log(x)- sum((log(beta1)-beta1_share)^2)/(2*x)
      sigma2.beta1 <- arms(sigma2.beta1, logden_sigma2.beta1, function(x) (x>0)*(x<a.beta1), 1)
      sigma2.beta1_t[ite] <- sigma2.beta1


      # gamma1_share
      # 这一段代码已检查完毕
      mu.n.gamma1_share <- mean(log(gamma1))
      sigma2.n.gamma1_share <- sigma2.gamma1/ngroup

      mean.gamma1_share <- (tau2.gamma1*mu.n.gamma1_share + sigma2.n.gamma1_share*mu.gamma1)/(tau2.gamma1+sigma2.n.gamma1_share)
      var.gamma1_share <- (tau2.gamma1*sigma2.n.gamma1_share)/(tau2.gamma1+sigma2.n.gamma1_share)
      gamma1_share <- rnorm(1, mean = mean.gamma1_share, sd=sqrt(var.gamma1_share))
      gamma1_share_t[ite] <- gamma1_share


      # sigma2.gamma1
      # 这一段代码已检查完毕
      logden_sigma2.gamma1 <- function(x) -ngroup/2*log(x) - sum((log(gamma1)-gamma1_share)^2)/(2*x)
      sigma2.gamma1 <- arms(sigma2.gamma1, logden_sigma2.gamma1, function(x) (x>0)*(x<a.gamma1), 1)
      sigma2.gamma1_t[ite] <- sigma2.gamma1

      #迭代过程至此完成完整的一轮，一共循环完成N.post轮
    }

    #取一个序列（index）来表明，只取烧录期之后的稳定样本
    ind <- seq((N.burnin+1),N.post)

    # 下面这行注释不是ywt写的
    # matrix.post <- cbind()[ind, ]
    matrix.uti <- array(list(),dim=c(1,1,ngroup))

    for(k in 1:ngroup){
      matrix.uti[[k]] = cbind(beta0_t[,k], beta1_t[,k], gamma0_t[,k],
                              gamma1_t[,k], alpha_t[,k], rho_t[,k])[ind, ]

      # matrix.uti 是一个列表，其中的每一个元素（matrix.uti[[k]]）是对应某个指标的参数矩阵。
      # matrix.uti[[k]]记录了k指标下，所有**表观信息层参数**的所有迭代轮次的取值
      # 以beta0_t为例子，这是一个Npost行*K列的数组，记录了每个指标k下，在Npost次迭代过程中的所有值
      # 所以这里等号右边是选取了第k个指标下的全体参数的Npost次的迭代值
    }
    return(matrix.uti)
  }

  #################

  ### 真值的预处理和部分变量的初始化
  #####################
  # 用 effdose 函数处理真实/已知数据
  ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
  umat <- ed$umat
  nd <- ed$nd
  vbeta <- ed$vbeta  # beta0.hat, beta1.hat
  vgamma <- ed$vgamma # gamma0.hat, gamma1.hat

  k.status <- rep(1, ngroup)  # 每组的状态，1 表示继续进行，0 表示停止
  dose.ind <- rep(1, ngroup)  # 每组当前的剂量索引
  h.dose <- rep(1, ngroup)    # 每组已尝试的最高剂量
  dat <- vector("list", ngroup)  # 存储各组的患者数据
  finish_day <- rep(1, ngroup) # 存储各组的结束试验日期

  #####################

  # 生成入组时间表格，先生成4*Nk的表格数据，再合并为一个长表格，按长表格里的患者依次入组/判定
  # 平均间隔设置：患者入组是poisson过程，平均入组天数分别为10,12,14,16天
  enrollment_times <- generate_enrollment_times(n.cohortsize, cohortsize, ngroup, mean_intervals=c(10, 12, 14, 16))
  enrolltimes_long_table <- merge_and_sort_enrollment_times(enrollment_times)
  # print(enrollment_times)
  # print(enrolltimes_long_table)
  total_patients <- n.cohortsize * cohortsize
  test_day = 0
  ### 开始实验
  row <- 1  # 初始化索引
  while (row <= nrow(enrolltimes_long_table)){
    test_day <- enrolltimes_long_table[row, 1]
    j <- enrolltimes_long_table[row, 2]  # j是group_id,当前患者是第j小组的
    patient_id <- enrolltimes_long_table[row, 3] # 当前患者是小组里第patient_id个入组治疗的

    if (k.status[j]==0){ row = row + 1 ; next } # 对于提前试验结束的小组 跳过之后的流程

    if (patient_id == 4) { # 第四人特殊处理 本if代码在12.14重写完善

      # message(paste("第", j, "组第 4 位患者在第", test_day, "天即将到来，进入特殊处理逻辑"))

      # 更新第j小组前三位患者的毒性结果 y_T 和疗效结果 y_E
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

      # 提取前三位患者的毒性数据
      previous_patients <- dat[[j]][1:3, ]

      # 判定是否至少2个患者已观察完无毒性，或至少1个患者已观察到出现毒性
      toxicity_2patients_no_occurred <- sum(previous_patients$y_T == 0, na.rm = TRUE) >= 2
      toxicity_occurred <- sum(previous_patients$y_T == 1, na.rm = TRUE) >= 1

      # 如果前三位患者的毒性结果不明确，推迟第 4 位患者的入组时间
      if (!toxicity_2patients_no_occurred && !toxicity_occurred) {

        # 计算前三位患者的观察完成时间
        completed_times <- ifelse(
          previous_patients$Y_T == 1,
          previous_patients$t_entry + previous_patients$T_T,  # 出现毒性
          previous_patients$t_entry + U_T                     # 无毒性
        )

        # 分毒性情况判定推迟时间 (000), (001), (011), (111)
        if (sum(previous_patients$Y_T == 0) == 3) {
          # (000) 所有患者均无毒性
          sorted_times <- sort(completed_times, decreasing = TRUE)
          delay_time <- sorted_times[2]  # 第二大的值
        } else if (sum(previous_patients$Y_T == 0) == 2 && sum(previous_patients$Y_T == 1) == 1) {
          # (001) 一位患者有毒性，其他两位无毒性
          max_no_toxicity_time <- max(completed_times[previous_patients$Y_T == 0])
          single_toxicity_time <- min(completed_times[previous_patients$Y_T == 1])
          delay_time <- min(max_no_toxicity_time, single_toxicity_time)
        } else if (sum(previous_patients$Y_T == 1) == 2) {
          # (011) 两位患者有毒性
          delay_time <- min(completed_times[previous_patients$Y_T == 1])
        } else if (sum(previous_patients$Y_T == 1) == 3) {
          # (111) 三位患者有毒性
          delay_time <- min(completed_times[previous_patients$Y_T == 1])
        } else {
          stop("Unexpected case: Toxicity data does not fit expected patterns.")
        }

        # 推迟到“前三个患者部分结果出现的时间+0.05”
        # print(dat[[j]])
        enrolltimes_long_table[row,1] <- delay_time + 0.05
        # message(paste("第", j, "组前三个患者数据尚不明确，第 4 位患者的入组时间已推迟至", enrolltimes_long_table[row,1], "天"))

        # 推迟入组时间后重新排序表格并重新分配 PatientID
        enrolltimes_long_table <- reorder_enrolltimes_table(enrolltimes_long_table)

        # 跳过当前循环，回到while继续当前第row行的患者处理！不执行row+1的操作！
        # 推迟入组后，row+1行的患者会变为row行，所以继续对第row行患者处理！
        # 即使推迟入组后的患者正巧还在第row行，那还是要再对他进行判定！
        next

      } else { # 如果前三位患者数据满足两个条件之一，则根据毒性结果进行剂量决策

        # 一个极端巧合情况，同时2位患者无毒性观察好+1位患者观察到毒性，时间相同，故应该优先判断毒性
        if (toxicity_occurred) {
          # 如果前三位患者中出现毒性，保持最低剂量
          dose.ind[j] <- 1
          # message(paste("第", j, "组前三位患者出现毒性结果，第4位患者保持最低剂量", dose.ind[j]))
        } else if (toxicity_2patients_no_occurred) {
          # 如果前三位患者为两位无毒性+一位无毒性/NA，则提升剂量
          dose.ind[j] <- min(h.dose[j] + 1, nd[j])  # 确保不超过最大剂量
          h.dose[j] <- dose.ind[j]
          # message(paste("第", j, "组前三位患者两位观察完毕无毒性+一位在观察中/无毒性，第4位患者提升剂量至", dose.ind[j]))
        }
      }

    } else if (patient_id != 1 && patient_id %% cohortsize == 1) {

      # mcmc运算前，将data_no_yT_NA初始化为一个长度为ngroup的list
      data_no_yT_NA <- vector("list", ngroup)  # 简化成list更清晰，不用array(list())
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
        # 考虑某一些极度病态的情况：万一在某一组即将入组第7个患者时，存在另一组仍未有患者，此时同样给该组一个无毒性的数据，否则mcmc3rd在某组无患者的情况下无法进行
        # dat保留着真实的数据，我们后续只能对dat_simulation增加模拟数据和进行mcmc3rd
        if(is.null(dat[[kk]]) || nrow(dat[[kk]]) == 0){
          # message(paste("第", kk, "组还没有患者到来，添加模拟数据来对第",j,"组的第",patient_id,"位患者进行第三轮mcmc和剂量决策"))
          # 添加模拟数据（未观察到毒性，观察期结束）
          u <- umat[[kk]]
          out.temp <- outcome(dose.ind[kk], cohortsize = 1, pt.mat[[kk]], vpt[kk, ], u, type, day = test_day)
          if (is.null(out.temp) || nrow(out.temp) == 0) {
            stop("out.temp 为空或无有效数据")
          }
          out.temp <- as.data.frame(out.temp)
          colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")

          # 添加 group 列
          out.temp$group <- kk
          out.temp$Y_T = 0                      # 未观察到毒性
          out.temp$Y_E = 0                      # 未观察到疗效
          out.temp$dose = dose.ind[kk]          # 使用最低剂量
          out.temp$T_T = Inf                    # 毒性观察期结束
          out.temp$T_E = Inf                    # 疗效观察期结束
          out.temp$y_T = 0                      # 毒性结果为无毒性
          out.temp$y_E = 0                      # 疗效结果为无疗效
          out.temp$t_entry = test_day-U_T       # 入组时间推到观察期开始前
          out.temp$t_follow = U_T               # 累积观察时间等于观察期

          dat_simulation[[kk]] <- out.temp
        }
      }

      # 此时第j组至少是即将到来第7个患者，前三个患者数据肯定有已知的,data_no_yT_NA[[j]]不会是空的
      # 病态情况下有一组会来的格外的慢！第kk组有患者在治疗中但是yT结果未知；data_no_yT_NA需要临时补充一行无毒数据，才能进行后续的一阶段mcmc
      # 检查 data_no_yT_NA 是否所有元素都不是 NULL 且行数大于0，不是的话则需要为其他组补充一行无毒数据
      all_elements_not_null_nrow_more_than0 <- all(sapply(data_no_yT_NA, function(x) !is.null(x) && nrow(x) > 0))
      # print(all_elements_not_null_nrow_more_than0)

      for (kk in 1:ngroup){
        if(is.null(data_no_yT_NA[[kk]]) || nrow(data_no_yT_NA[[kk]]) == 0){
          # 添加模拟数据（未观察到毒性，观察期结束）
          # message(paste("第", kk, "组没有患者完成毒性观察，添加模拟数据来对第",j,"组的第",patient_id,"位患者进行第一轮mcmc和剂量决策"))
          u <- umat[[kk]]
          out.temp <- outcome(dose.ind[kk], cohortsize = 1, pt.mat[[kk]], vpt[kk, ], u, type, day = test_day)
          if (is.null(out.temp) || nrow(out.temp) == 0) {
            stop("out.temp 为空或无有效数据")
          }
          out.temp <- as.data.frame(out.temp)
          colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")

          # 添加 group 列
          out.temp$group <- kk
          out.temp$Y_T = 0                      # 未观察到毒性
          out.temp$Y_E = 0                      # 未观察到疗效
          out.temp$dose = dose.ind[kk]          # 使用最低剂量
          out.temp$T_T = Inf                    # 毒性观察期结束
          out.temp$T_E = Inf                    # 疗效观察期结束
          out.temp$y_T = 0                      # 毒性结果为无毒性
          out.temp$y_E = 0                      # 疗效结果为无疗效
          out.temp$t_entry = test_day-U_T       # 入组时间推到观察期开始前
          out.temp$t_follow = U_T               # 累积观察时间等于观察期

          data_no_yT_NA[[kk]] <- out.temp
        }
      }

      # print(dat_simulation)
      # print(data_no_yT_NA)
      all_dat <- do.call(rbind, dat_simulation)

      # 第一步：对data_no_yT_NA进行mcmc循环得到mcmc.temp_1st
      # 并提取出初步的后验参数beta_0_1st,beta_1_1st,用于下一步的似然函数
      ## 目标函数的调用
      mcmc.temp_1st <- mcmc.arms_Delayed_BPCC_1st(data_no_yT_NA, vbeta=vbeta, vgamma=vgamma, ngroup=ngroup, type)
      matrix.uti <- mcmc.temp_1st  # 提取 matrix.uti # 假设目标函数的输出是 matrix.uti
      beta_0_1st <- numeric(ngroup) # 初始化存储 beta_0 和 beta_1 均值的向量
      beta_1_1st <- numeric(ngroup)

      for (kk in 1:ngroup) { # 遍历每组数据，计算 beta_0 和 beta_1 的均值
        # 检查 matrix.uti[[kk]] 是否为有效矩阵
        if (!is.null(matrix.uti[[kk]]) && nrow(matrix.uti[[kk]]) > 0) {
          # 计算均值
          beta_0_1st[kk] <- mean(matrix.uti[[kk]][, 1])  # 第 1 列是 beta_0
          beta_1_1st[kk] <- mean(matrix.uti[[kk]][, 2])  # 第 2 列是 beta_1
        } else {
          warning(paste("mcmc.temp_1st第", kk, "组的 matrix.uti[[kk]] 为空或无有效数据"))
        }
      }

      # 第二步，书写并求解极大似然函数
      ## 全体已知数据放一起，得到第test_day天时的all_dat

      ################
      ## 全体已入组的毒性相关的数据的似然函数
      # 新的优化后的log_L函数
      log_L <- function(params, all_dat, beta_0_1st, beta_1_1st) {
        # 参数检查
        if (length(params) != 3) stop("params must have 3 elements: c, eta_0, eta_1")
        if (!all(c("y_T", "T_T", "t_follow", "dose", "group") %in% names(all_dat)))
          stop("all_dat must contain columns: y_T, T_T, t_follow, dose, group")

        # 提取参数
        c <- params[1]
        eta_0 <- params[2]
        eta_1 <- params[3]

        # 提取数据
        y_T <- all_dat$y_T
        tilde_T_T <- all_dat$T_T / U_T
        tilde_t_follow <- all_dat$t_follow / U_T # 只有y_T=NA时,会用到tilde_t_follow, y_T=1时直接用tilde_T_T
        d <- all_dat$dose
        group <- all_dat$group

        # 计算每个患者的 F 值(按组别提取 beta 参数). 只有F值是和组别有关的，wT和WT的值都与组别无关
        F_vals <- pnorm(beta_0_1st[group] + beta_1_1st[group] * d)

        # 计算 w_T 和 W_T（支持向量化）在wT，WT定义里有说，即使tilde_t输入值为Inf也不会报错！
        w_T_vals <- w_T(tilde_T_T, d, c, eta_0, eta_1)
        W_T_vals <- W_T(tilde_t_follow, d, c, eta_0, eta_1)

        # 检查向量长度一致性
        if (length(w_T_vals) != length(y_T) || length(W_T_vals) != length(y_T) || length(F_vals) != length(y_T)) {
          warning(paste("Vector lengths of w_T_vals, W_T_vals, F_vals, and y_T must match"))
        }

        # 情况 (1,1): 毒性反应 -> g = w^T * F
        log_g <- rep(0, length(y_T))
        valid_idx <- !is.na(y_T) & y_T == 1 # 筛选 y_T 不为 NA 且等于 1 的索引
        log_g[valid_idx] <- log(w_T_vals[valid_idx] * F_vals[valid_idx]) # 仅对有效位置赋值，其他位置保持为0

        # # 打印有多少毒性事件患者进入似然
        # cat("有毒性发生的患者数：", sum(valid_idx), "\n")

        # 情况 (0,0): 毒性未知 -> 1 - W^T * F
        log_survival <- rep(0, length(y_T))
        log_survival[is.na(y_T)] <- log(1 - W_T_vals[is.na(y_T)] * F_vals[is.na(y_T)]) # 仅计算 y_T = NA 的情况，其他位置保持为0

        # 汇总对数似然
        log_L_total <- sum(log_g+log_survival, na.rm = TRUE) #求和时忽略缺失值

        # 12.11加上
        # 假设先验均值与方差（需根据你的实际情况设定）
        # mu_c <- 1.5; sigma_c <- test_day/100
        # mu_eta0 <- 0; sigma_eta0 <- test_day/100
        # mu_eta1 <- 1.5; sigma_eta1 <- test_day/100

        # 惩罚项（高斯先验的负对数似然对参数的贡献）
        # penalty <- ((c - mu_c)^2/(2*sigma_c^2)) + ((eta_0 - mu_eta0)^2/(2*sigma_eta0^2)) + ((eta_1 - mu_eta1)^2/(2*sigma_eta1^2))
        penalty = 0

        return(log_L_total-penalty)
      }

      # 优化 log_L，求得 c, eta_0, eta_1
      optim_res <- optim(
        par = c(1.5, 0, 1.5),  # 初始值 (c, eta_0, eta_1)
        fn = function(params) -log_L(params, all_dat, beta_0_1st, beta_1_1st),  # 使用合并后的数据框
        method = "L-BFGS-B",
        lower = c(0.1, -1, 0.5),  # 三个参数分别的搜索下界：c > 0，暂定的这三个下界应该没问题，毕竟真值分别是1，0，2
        upper = c(4, 2, 4),
      )

      # 打印优化结果
      # cat("优化后的参数:\n")
      # print(optim_res$par)

      # 保存优化结果
      c_hat <- optim_res$par[1]
      eta_0_hat <- optim_res$par[2]
      eta_1_hat <- optim_res$par[3]


      ################

      # 第三步：再一次mcmc循环，使用mcmc_3rd函数得到mcmc.temp_3rd，用于第4步的剂量决策
      mcmc.temp_3rd <- mcmc.arms_Delayed_BPCC_3rd(dat_simulation, vbeta=vbeta, vgamma=vgamma, ngroup=ngroup, type,
                                                  c_hat, eta_0_hat, eta_1_hat)


      # 第四步：基于后一次的mcmc循环结果mcmc.temp_3rd，对即将入组新患者的这个第k组更新剂量决策
      u=umat[[j]]
      summ.temp<-summary.mcmc(matrix.uti=mcmc.temp_3rd[[j]],n.dose=nd[j], u=u, type)
      tox.summ <- summ.temp$tox

      if (nrow(dat[[j]]) == n.cohortsize * cohortsize){h.dose[j]=nd[j]} # 防止病态数据在假患者处才第一次尝试最高剂量
      if ((tox.summ[h.dose[j]]>.5) & (h.dose[j]!=nd[j])){ # escalation conditions
        dose.ind[j] <- h.dose[j]+1
        h.dose[j] <- dose.ind[j]
        # message(paste("当前已试验的最高剂量毒性低于0.3的概率高于0.5且未达试验设计的最大剂量,第",j,"组提升至下一个剂量",dose.ind[j]))
      }else{

        eff.summ <- summ.temp$eff
        capA.ind <- sort(which((tox.summ>C_T)&(eff.summ>C_E)))

        if(length(capA.ind)==0){ # empty
          dose.ind[j] = 0
          k.status[j] = 0
          finish_day[j] = test_day
          # message(paste('第',j,'组在第',patient_id,'位患者即将到来时提前结束，不选出OBD'))
        }else{# non empty

          uti <- summ.temp$uti[capA.ind]
          # 抓取所有的允许剂量指标的得分，即将用于随机化试验，分配不同剂量被选到的概率
          # if(patient_id==n.cohortsize*cohortsize+1){ # reached the maximum sample size
          #   print(paste("第", j,"组的模拟得分是", uti))
          #   dose.ind[j] = capA.ind[which.max(uti)]
          #   #这里的dose.ind有两个作用，一个是担任分配剂量时的临时中间值，随后将分配给outcome中入组患者的剂量
          #   #另一个作用是是在所有的ncohort都治疗结束后，担任展现最后的最优剂量
          #   k.status[j] = 0
          #   print(paste("第", j,"组已经入满了", n.cohortsize*cohortsize,"名患者，不再入组"))
          # }else{
          # randomize patients within admissible doses
          # 后续这里可以改为选取最高得分的剂量！！！～～～～～～～～～～
          uti.normalize <- uti/sum(uti)
          cumu.uti <- uti

          # 每个剂量的累积得分函数cumu.uti，用于随机化实验，看看随机数U在哪个区间
          for (i in 1:length(uti)) cumu.uti[i] <- sum(uti.normalize[1:i])
          r <- runif(1,0,1)
          dose.ind[j] <- min(h.dose[j]+1,capA.ind[min(which(r<cumu.uti))]) # 随机化分配时不能直接高过两个最高剂量分配
          # message(paste("第",j,"组在第",test_day,"天对第",patient_id,"位患者根据得分执行自适应剂量决策，剂量为",dose.ind[j]))


          # 更新最高剂量
          if(dose.ind[j]>h.dose[j]){
            dose.ind[j] <- h.dose[j]+1
            h.dose[j] <- dose.ind[j]
          }

          # 若patient_id=n.cohortsize*cohortsize+1,则赋值dose.ind[j]为第k组的最优剂量并跳过生成患者数据
          if (nrow(dat[[j]]) == n.cohortsize * cohortsize){
            dose.ind[j] = capA.ind[which.max(uti)]
            k.status[j] = 0
            finish_day[j] = test_day
            # message(paste("第", j,"组已经入满了", n.cohortsize*cohortsize,"名患者且观察到了所有数据，最优剂量为",dose.ind[j]))
            # message(paste("第", j,"组各剂量的得分分别是",summ.temp$uti))
            # message(paste("第", j,"组试验进行天数为",finish_day[j],"天"))
          }
        }
      }
    }  # 第j组的剂量决策步骤至此结束


    #下面开始对这第patient_id位患者入组，并更新dat[[j]]

    # 若当前组被鉴定为无效，或者即将到来的是最后一个假患者，则跳过生成患者数据，到下一行患者入组判定
    if ( k.status[j] == 0) { row = row + 1; next }

    u <- umat[[j]]
    out.temp <- outcome(dose.ind[j], cohortsize = 1, pt.mat[[j]], vpt[j, ], u, type, day = test_day)

    if (is.null(out.temp) || nrow(out.temp) == 0) {
      stop("out.temp 为空或无有效数据")
    }

    out.temp <- as.data.frame(out.temp)
    colnames(out.temp) <- c("Y_T", "Y_E", "dose", "T_T", "T_E", "y_T", "y_E", "t_entry", "t_follow")

    # 添加 group 列
    out.temp$group <- j

    # 初始化或追加 dat[[j]]
    if (is.null(dat[[j]]) || nrow(dat[[j]]) == 0) {
      dat[[j]] <- out.temp
    } else {
      dat[[j]] <- rbind(dat[[j]], out.temp)
    }

    # 如果当前入组的是第 n.cohortsize * cohortsize 个患者，则此时可以对假患者的入组时间赋值
    # 计算第 j 小组所有病人观察到疗效结果的最大值

    if (patient_id==n.cohortsize * cohortsize){
      j_group_efficacy_times <- ifelse(
        dat[[j]]$Y_E == 1,
        dat[[j]]$t_entry + dat[[j]]$T_E,  # 出现疗效
        dat[[j]]$t_entry + U_E            # 无疗效
      )

      # 找到第j小组假患者在 enrolltimes_long_table 中的位置

      fake_patient_row <- which(
        enrolltimes_long_table$Group == j &
          enrolltimes_long_table$PatientID == ( n.cohortsize * cohortsize + 1 )
      )

      if (length(fake_patient_row) == 1) {
        # message(paste("第", j, "组假患者的入组时间已赋值为", max(j_group_efficacy_times, na.rm = TRUE)))
        enrolltimes_long_table[fake_patient_row,1] <- max(j_group_efficacy_times, na.rm = TRUE)
        enrolltimes_long_table <- reorder_enrolltimes_table(enrolltimes_long_table)

      } else {
        # 打印调试信息
        # print((enrolltimes_long_table))  # 打印表格查看是否初始化问题
        # print(paste("Group:", j, "PatientID:", n.cohortsize * cohortsize + 1))  # 打印未找到的组和患者ID
        stop(paste("未找到第", j, "组的假患者记录，可能数据有误"))
        # message(paste("未找到第", j, "组的假患者记录，可能数据有误"))
      }
    }
    # 当前row的患者入组完毕，跳到下一行患者入组判定
    row = row + 1

    # 测试断点代码：print(paste('dat[[j]]:', dat[[j]]))

  } # enrolltimes_long_table里每一个患者的决策判定+入组步骤至此结束

  # 全体试验流程已结束，下面开始数据整理

  ### 循环结束，输出最后循环的参数均值和剂量分配序列
  # 把长表格还原为4*n.cohortsize的表格，最后一行命名为finish.day，在return里输出
  enrolltimes_final = convert_to_named_original_format(enrolltimes_long_table, ngroup, total_patients)

  # 参数均值
  # 这里的 mcmc.temp_3rd 是 mcmc.arms_Delayed_BPCC_3rd 函数的最后一次输出，所有患者数据全已知
  matrix_means <- lapply(mcmc.temp_3rd, function(mat) {
    if (!is.null(mat) && nrow(mat) > 0) {
      colMeans(mat)  # 对每个矩阵的列计算均值
    } else {
      NA  # 如果矩阵为空，返回 NA
    }
  })

  # 转换为数据框并添加列名
  matrix_means_df <- do.call(rbind, lapply(seq_along(matrix_means), function(k) {
    data.frame(
      Indicator = k,
      t(matrix_means[[k]])
    )
  }))

  # 添加列名到结果数据框
  colnames(matrix_means_df) <- c(
    "Indicator",  # 指标编号
    "beta0", "beta1", "gamma0", "gamma1", "alpha", "rho"  # 参数列名
  )

  # 打印结果
  # print(matrix_means_df)

  # 输出每个第k组的OBD剂量选择和所有轮次的剂量分配序列
  res <- vector("list", ngroup)
  for(k in 1:ngroup){
    dtemp <- NULL
    # length(dat[[k]][,3]):试验结束时，第k组参与治疗患者的数量（有些组可能中途结束试验）
    for(id in 1:length(dat[[k]][,3])){
      dtemp[id] = vd[[k]][which(umat[[k]]==dat[[k]][,3][id])]
      # dat[[k]][,3][id],获取第k组第id位患者的标准化剂量值。
      # which 返回满足条件的索引，即 umat[[k]] 中与当前患者的标准化剂量值相等的位置。
      # 查找与标准化剂量umat[[k]](index)对应的原始剂量vd[[k]](index)，并将其赋值给dtemp[id]
    }
    res[[k]] = list(d.select = dose.ind[k], d.alloc=dtemp)
  }
  return(list(dose_pick=res,parameter_mean=matrix_means_df,finish_day=finish_day,enroll=enrolltimes_final,dat=dat))
}
