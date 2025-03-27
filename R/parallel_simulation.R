#' Perform Parallel Simulations for Bayesian Platform Data-augmentation
#'
#' This function performs multiple parallel simulations of the Bayesian Platform
#' Data-augmentation process for delayed outcomes in Phase I/II clinical trials.
#' The function executes a series of simulations using the main simulation function
#' `main_BPDD` in parallel across multiple cores, generating various statistics and
#' saving results to CSV files.
#'
#' @param num_simulations An integer specifying the number of simulations to run.
#'        Default is 4.
#' @param num_cores An integer specifying the number of CPU cores to use for parallel processing.
#'        Default is 4.
#'
#' @return This function does not return a value directly. It saves results to CSV files
#'         in the working directory.
#' @export
parallel_simulation <- function(num_simulations = 4, num_cores = 4) {
  # 加载必要的包
  library(parallel)
  library(tictoc)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(1004)

  # 定义并行任务函数
  run_simulation <- function(seed, ngroup, num_simulations, num_cores) {
    set.seed(seed)

    # === 加载真值数据 ===
    true_values_path <- system.file("extdata", "true_values.RData", package = "BPDD")
    if (true_values_path == "") {
      stop("Failed to find true values file in the package.")
    }

    # 加载真值数据
    load(true_values_path)  # 这将加载所有的真值，包括 cT_true, eta_0T_true, eta_1T_true, vd, pt.mat, prior.pt.mat, prior.pe.mat 等

    # 运行主函数
    result <- main_BPDD(
      n.cohortsize = 14, cohortsize = 3, phi.T = 0.3, phi.E = 0.3,
      pt.mat = pt.mat, prior.pt.mat = prior.pt.mat, prior.pe.mat = prior.pe.mat,
      vpt = vpt, vd = vd, type = "binary", C_T = 0.05, C_E = 0.05
    )

    return(result)
  }

  # === 并行计算 ===
  tic.clearlog()  # 清空计时日志
  tic("Parallel Simulations")

  # 使用 mclapply 并行运行
  all_results <- mclapply(
    1:num_simulations,
    function(seed) {
      tryCatch(
        {
          # 生成并记录当前种子
          .Random.seed <- parallel::nextRNGStream(.Random.seed)  # 独立随机数流
          current_seed <- .Random.seed  # 保存当前的随机种子

          # 打印当前种子信息到控制台
          cat("Running simulation for seed:", seed, "\n")

          # 保存当前种子到文件，用于调试
          seed_log <- paste("Seed for simulation", seed, ":", paste(current_seed, collapse = " "), "\n")
          write(seed_log, file = "all_seeds.log", append = TRUE)

          # 调用并行任务函数并传递必要的变量
          run_simulation(seed, ngroup = 4, num_simulations = num_simulations, num_cores = num_cores)
        },
        error = function(e) {
          # 记录错误的种子和信息
          error_log <- paste("Error in simulation for seed:", seed, "\n", e$message, "\n")
          write(error_log, file = "error_seeds.log", append = TRUE)
          return(NULL)
        }
      )
    },
    mc.cores = num_cores,
    mc.preschedule = FALSE
  )

  execution_time <- toc()  # 获取运行时间

  # 保存运行时间为 CSV 文件
  execution_times <- data.frame(
    Task = execution_time$msg,
    Elapsed_Time_Sec = execution_time$toc - execution_time$tic
  )
  write.csv(
    execution_times,
    "execution_times_summary.csv",
    row.names = FALSE
  )
  cat("\nExecution times saved to 'execution_times_summary.csv'\n")


  # 保存 all_results 到 Rdata 文件
  save(all_results, file = "all_results.Rdata")
  cat("\nAll simulation results saved to 'all_results.Rdata'\n")

  # === 结果统计 ===
  cat("\n=== Result Statistics ===\n")

  # 初始化统计变量
  dose_counts <- matrix(0, nrow = ngroup, ncol = 4)  # 剂量选择统计
  finish_days <- matrix(NA, nrow = num_simulations, ncol = ngroup)  # 存储所有模拟的 finish_day
  all_enrollment_times <- list()  # 存储所有模拟的 enrollment_times
  all_dat <- list()  # 存储所有模拟的 dat

  # 遍历结果提取数据
  for (seed in 1:num_simulations) {
    if (!is.null(all_results[[seed]]) && is.list(all_results[[seed]]$dose_pick)) {
      for (group in 1:ngroup) {
        selected_dose <- all_results[[seed]]$dose_pick[[group]]$d.select
        if (!is.null(selected_dose) && selected_dose > 0) {
          dose_counts[group, selected_dose] <- dose_counts[group, selected_dose] + 1
        }
      }
      # 保存 finish_day 信息
      finish_days[seed, ] <- all_results[[seed]]$finish_day

      # 保存 enrollment_times 信息
      enrollment_time <- as.data.frame(all_results[[seed]]$enroll)
      enrollment_time$Simulation <- seed
      all_enrollment_times[[seed]] <- enrollment_time

      # 保存 dat 信息
      dat_info <- lapply(seq_along(all_results[[seed]]$dat), function(group_id) {
        dat_group <- all_results[[seed]]$dat[[group_id]]
        if (!is.null(dat_group) && nrow(dat_group) > 0) {
          dat_group$Simulation <- seed
          dat_group$Group <- group_id
          return(dat_group)
        } else {
          return(NULL)
        }
      })
      all_dat[[seed]] <- do.call(rbind, dat_info)
    }
  }

  # 整理所有 enrollment_times 并保存为 CSV 文件
  all_enrollment_times_df <- do.call(rbind, all_enrollment_times)
  write.csv(all_enrollment_times_df, "enrollment_times_summary.csv", row.names = FALSE)
  cat("\nEnrollment times saved to 'enrollment_times_summary.csv'\n")

  # 整理所有 dat 数据并保存为 CSV 文件
  all_dat_df <- do.call(rbind, all_dat)
  write.csv(all_dat_df, "dat_summary.csv", row.names = FALSE)
  cat("\nDat information saved to 'dat_summary.csv'\n")

  # 保存 finish_day 信息
  write.csv(data.frame(Simulation = 1:num_simulations, finish_days), "finish_days_summary.csv", row.names = FALSE)
  cat("\nFinish days saved to 'finish_days_summary.csv'\n")

  # 打印并保存剂量统计结果
  for (group in 1:ngroup) {
    cat("\nGroup", group, "dose selection counts:\n")
    print(dose_counts[group, ])
  }

  # 将剂量选择统计保存为 CSV 文件
  write.csv(
    as.data.frame(dose_counts),
    "dose_selection_counts_summary.csv",
    row.names = TRUE
  )
  cat("\nDose selection counts saved to 'dose_selection_counts_summary.csv'\n")

  # 提取所有参数均值
  all_means <- do.call(rbind, lapply(1:num_simulations, function(seed) {
    if (!is.null(all_results[[seed]]) && is.data.frame(all_results[[seed]]$parameter_mean)) {
      data <- all_results[[seed]]$parameter_mean
      data$Simulation <- seed  # 添加模拟编号
      return(data)
    } else {
      return(NULL)
    }
  }))

  # 保存参数均值结果
  write.csv(all_means, "parameter_means_summary.csv", row.names = FALSE)
  cat("\nParameter means saved to 'parameter_means_summary.csv'\n")

  # === 整理剂量分配数据并保存 ===
  all_data <- data.frame()
  for (k in 1:num_simulations) {
    if (!is.null(all_results[[k]]$dose_pick)) {
      for (j in 1:length(all_results[[k]]$dose_pick)) {
        group_data <- all_results[[k]]$dose_pick[[j]]
        group_df <- data.frame(
          Simulation = k,
          Group = j,
          d.select = group_data$d.select,
          d.alloc = group_data$d.alloc
        )
        all_data <- rbind(all_data, group_df)
      }
    }
  }

  # 保存剂量分配结果
  write.csv(all_data, "dose_allocation_summary.csv", row.names = FALSE)
  cat("\nDose allocation saved to 'dose_allocation_summary.csv'\n")

  cat("\nParallel results saved to CSV files.\n")
}
