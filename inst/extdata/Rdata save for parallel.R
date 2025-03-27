##### 真值赋值区
#################

# 新补充的真值赋值：cTE,eta_0TE,eta_1TE真值,T/E窗口期的赋值
cT_true=1; eta_0T_true=0; eta_1T_true=2;
cE_true=2; eta_0E_true=0; eta_1E_true=2;

# 第6种情况，更改WT，WE的真实值，再延迟一会毒性效应，这应该会延长几天的试验时间，主要是在第四天的时候
# cT_true=1; eta_0T_true=-0.5; eta_1T_true=2;
# cE_true=2; eta_0E_true=-0.5; eta_1E_true=2;

U_T=30; U_E=90;
C_T = 0.05; C_E = 0.05;

# 原先的真值代码
ngroup = 4 # number of indications
pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialization


# n.cohortsize = 14 # the total number of cohorts  14 for 场景一～场景三 12 for 场景四
cohortsize = 3 # the cohort size
phi.T = 0.3 # upper limit of the toxicity rate
phi.E = 0.3 # lower limit of the efficacy rate, binary Y_E
type = 'binary' # Y^E is binary

######## 需要改动的场景1234
# 真实输入 场景一OBD1234
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

#
# # 真实输入 场景二OBD：22(23)3
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

#
# # 真实输入 场景三OBD:0120
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


# # 真实输入 场景四OBD:0(12)23
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


# 初始化 prior.pt.mat 和 prior.pe.mat
prior.pt.mat <- vector("list", length(vd))
prior.pe.mat <- vector("list", length(vd))

# 毒性/疗效先验值
for (i in 1:ngroup){
  if(length(vd[[i]])==4){
    prior.pt.mat[[i]] = c(0.10, 0.20, 0.30, 0.40)
    prior.pe.mat[[i]] = c(0.15, 0.30, 0.45, 0.40)
  } else if (length(vd[[i]])==3) {
    prior.pt.mat[[i]] = c(0.10, 0.20, 0.30)
    prior.pe.mat[[i]] = c(0.15, 0.30, 0.45)
  }
}

# 保存真值数据到 RData 文件
save(
  cT_true, eta_0T_true, eta_1T_true, cE_true, eta_0E_true, eta_1E_true, 
  U_T, U_E, C_T, C_E,
  ngroup, cohortsize, phi.T, phi.E, type,
  vd, pt.mat, pe.mat, uti.mat, prior.pt.mat, prior.pe.mat, vpt,
  file = "~/Documents/BPDD/inst/extdata/true_values.RData"
)