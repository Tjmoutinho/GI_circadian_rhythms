### Functions ###

read_kronos_batch <- function(file_names, sheet_nums, starts, ends, steps){
  n = 1
  all_cond <- list()
  control <- read_kronos(file_names[1], sheet_nums[1], starts[1], ends[1], steps[1])
  all_cond[[n]] <- mutate(control, lum_norm = y / max(abs(y))) #normalizing data
  n=n+1
  for (i in 2:length(file_names)){
    temp_df <- read_kronos(file_names[i], sheet_nums[i], starts[i], ends[i], steps[i])
    all_cond[[n]] <- mutate(temp_df, lum_norm = y / max(abs(control$y))) #normalizing data
    n = n+1
  }
  return(all_cond)
}

read_kronos <- function(file_name, sheet_num, start, end, step){
  df_A <- read_excel(file_name, sheet = sheet_num) 
  names(df_A) <- c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8")
  df_A <- df_A[, c("col_3","col_4")] %>% 
    .[19:length(.$col_4),] %>% 
    .[complete.cases(.),]
  df_A_interp <- as.data.frame(spline(df_A, xout = seq(from = start, to = end, by = step)))
  # df_A_2 <- df_A[ , c("X__1","X__2")]
  # df_A_3 <- df_A_2[20:length(df_A_2$X__1), ]
  # df_A_4 <- df_A_3[complete.cases(df_A_3),]
  #df_F <- read_excel("20171003per2 metabolites.xlsx", sheet = 6) %>% .[, c("X__1","X__2")] %>% .[20:length(df_A_2$X__2),] %>% .[complete.cases(df_A_3),]
  
  # #normalizing data
  # df_A_norm <- mutate(df_A_interp, lum_norm = y / max(abs(y)))
  return(df_A_interp)
}

PRT <- function(df_norm){
  #pull out peaks and troughs in data
  peaks <- findPeaks(df_norm$y)-1
  # peaks_df <- df_A_norm[peaks, ]
  
  troughs <- findValleys(df_norm$y)-1
  
  #pull out roots
  n <- 1
  roots <- c()
  for (i in 1:(length(df_norm$y)-1)){
    if (df_norm$y[i] < 0 && df_norm$y[i+1] > 0 ){
      roots[n] <- i 
      n <- n+1
    }
    else if (df_norm$y[i] > 0 && df_norm$y[i+1] < 0){
      roots[n] <- i 
      n <- n+1
    }
  }
  
  # Add the logic Column
  df_norm$logic <- rep(NA, length(df_norm$x))
  df_norm$logic[peaks] <- 1
  df_norm$logic[troughs] <- -1
  df_norm$logic[roots] <- 0
  
  prt <- do.call("rbind", list(df_norm[peaks, ], df_norm[troughs, ], df_norm[roots, ])) %>% .[order(.$x), ]
  
  return(prt)
}

PRT_batch <- function(all_cond){
  all_prt <- list()
  for(j in 1:length(all_cond)){
    prt_single <- PRT(all_cond[[j]])
    all_prt[[j]] <- prt_single
  }
  return(all_prt)
}

AMP <- function(df_norm){
  #pull out peaks and troughs in data
  peaks_data <- df_norm$lum_norm[df_norm$logic == 1]
  peak_amplitude <- ((sum(peaks_data)))/length(peaks_data)
  
  troughs_data <- df_norm$lum_norm[df_norm$logic == -1]
  trough_amplitude <- abs(((sum(troughs_data)))/length(troughs_data))
  
  #calculate amplitude
  all_amplitude <- ((sum(abs(c(peaks_data ,troughs_data))))/length(c(peaks_data ,troughs_data)))
  
  #Write values into new dataframe
  amplitude_data <- data.frame(amp = c(all_amplitude, peak_amplitude, trough_amplitude))
  
  #add a column for the names
  amplitude_data$names <- c("Ave Amplitude", "Ave Peak Amp", "Ave Trough Amp")
  
  return(amplitude_data)
}

AMP_batch <- function(all_prt){
  all_AMP <- list()
  for(j in 1:length(all_prt)){
    AMP_single <- AMP(all_prt[[j]])
    all_AMP[[j]] <- AMP_single
  }
  return(all_AMP)
}

PSHIFT <- function(control.samp, exp.samp) {
  ## Detailed notes on this code can be found in Lab Notebooke #1 for Deborah H. Luzader on pages 132 and 133
  len <- c(length(exp.samp$x),length(control.samp$x)) %>% min() #this will later be used to make the lengths the same
  r <- ccf(exp.samp$logic, control.samp$logic, lag.max = NULL, 
           type = c("correlation"),
           plot = FALSE)    #Change to TRUE to see the logic graph
  ind.max <- r$acf %>% abs() %>% which.max()
  sign.max <- r$acf[ind.max] / abs(r$acf[ind.max])
  if(r$lag[ind.max] == 0 && sign.max == 1){
    lag.max <- 0    #this is the case where the peaks roots and troughs already align
  } else if((r$lag[ind.max] == -1 && sign.max == 1) || (r$lag[ind.max] == 1 && sign.max == -1)) {
    lag.max <- 1    #this is the case where the peaks roots and troughs are shifted to the left relative to the control. (EX: control data is root, peak, root, trough. sample data is peak, root, trough, root)
  } else if((r$lag[ind.max] == -1 && sign.max == -1) || (r$lag[ind.max] == 1 && sign.max == 1)) {
    lag.max <- -1   #this is the case where the peaks roots and troughs are shifted to the right relative to the control. (EX: control data is root, peak, root, trough. sample data is trough, root, peak, root)
  } else if(r$lag[ind.max] == 0 && sign.max == -1){lag.max <- 2}   #this is the case where the peaks and troughs are reversed relative to the control, so the data is antiphasic. (EX: control data is root, peak, root, trough. sample data is root, trough, root, peak)
  # exp.samp.peaks <- subset(exp.samp, logic == "1")
  # control.samp.peaks <- subset(control.samp, logic == "1")
  if(lag.max == 0){
    phase.shift <- exp.samp$x[1:len] - control.samp$x[1:len]
  } else if(lag.max == 1){
    phase.shift <- exp.samp$x[1:(len-lag.max)] - control.samp$x[(1+lag.max):len]
  } else if(lag.max == -1){
    phase.shift <- exp.samp$x[(1-lag.max):len] - control.samp$x[1:(len+lag.max)]
  } else if(abs(lag.max) == 2){
    
    phase.shift.pos <- exp.samp$x[(1+lag.max):(len)] - control.samp$x[1:(len-lag.max)]
    phase.shift.neg <- control.samp$x[(1+lag.max):(len)]- exp.samp$x[1:(len-lag.max)]
    # phase.shift.neg <- exp.samp$x[1:(len-lag.max)] - control.samp$x[(1+lag.max):(len)]
    if((sum(phase.shift.pos)) < (sum(phase.shift.neg))){
      phase.shift <- phase.shift.pos
    } else if ((sum(phase.shift.neg)) < (sum(phase.shift.pos))){
      phase.shift <- exp.samp$x[1:(len-lag.max)] - control.samp$x[(1+lag.max):(len)]
    }
    # print("ERROR")
  }
  # plot(phase.shift)
  # exp.samp.peaks$x[1]-control.samp.peaks$x[1]
  # exp.samp.peaks$x[2]-control.samp.peaks$x[2]
  return(phase.shift)
}

WRITE_DATA_FILES <- function(file_name, sheet_nums, start_value, end_value, output_file1, output_file2, output_file3, output_file4){
  
  number_of_sheets <- length(sheet_nums)
  file_names <- rep(file_name, number_of_sheets)
  #
  starts <- rep(start_value, number_of_sheets)
  #
  ends <- rep(end_value, number_of_sheets)
  steps <- rep(0.01, number_of_sheets)
  
  #### Read data from excel into data frames ####
  all_cond <- read_kronos_batch(file_names, sheet_nums, starts, ends, steps)
  
  lum_norm_all <- data.frame("X" = all_cond[[1]][,"x"])
  for (i in 1:number_of_sheets) {
    lum_norm_all[,i+1] <- all_cond[[i]][,"lum_norm"]
  }
  colnames(lum_norm_all) <- c("X", as.character(sheet_nums))
  write.xlsx2(lum_norm_all, output_file1, sheetName = "Norm_Lum_data", row.names=FALSE)
  
  raw_lum_all <- data.frame("X" = all_cond[[1]][,"x"])
  for (i in 1:number_of_sheets) {
    raw_lum_all[,i+1] <- all_cond[[i]][,"y"]
  }
  colnames(raw_lum_all) <- c("X", as.character(sheet_nums))
  write.xlsx2(raw_lum_all, output_file1, sheetName = "Raw_Lum_data", row.names=FALSE, append = TRUE)
  
  #calculate peaks, roots, troughs
  all_prt <- PRT_batch(all_cond)
  
  # print all_prt to output_file2
  write.xlsx2(all_prt[[1]], output_file2, sheetName = as.character(sheet_nums[1]), row.names=FALSE)
  for (i in 2:number_of_sheets) {
    df <- all_prt[[i]] 
    write.xlsx2(df, output_file2, sheetName = as.character(sheet_nums[i]), row.names=FALSE, append = TRUE)
  }
  
  # Test and fix prt dataframes
  prt_sheets_to_remove <- c()
  for (j in 1:length(all_prt)) {
    temp_logic <- all_prt[[j]]$logic
    for(k in 1:(length(temp_logic)-1)){
      if (temp_logic[k]*temp_logic[k+1] != 0){
        # print out prt[[j]]
        View(all_prt[[j]])
        n <- as.integer(readline(prompt="Do you want to keep this data? Enter 1 for yes or 2 for no: "))
        if (n == 1){
          # Remove rows specified by user
          rows_to_remove <- as.numeric(unlist(strsplit((readline(prompt="Which single row do you want to delete? (Will prompt again if more than one row needs to be deleted.) ")), split=" ")))
          all_prt[[j]] <- all_prt[[j]][-rows_to_remove,]
        }
        else if(n == 2){
          prt_sheets_to_remove <- c(prt_sheets_to_remove, j)
          break
        }
      }
    }
  }

  # delete prt[[prt_sheets_to_remove]]
  if(length(prt_sheets_to_remove) != 0){
    print(c("Deleted prt:",prt_sheets_to_remove))
    all_prt[[prt_sheets_to_remove]] <- NULL
    # edit sheet numbers
    sheet_nums <- sheet_nums[-prt_sheets_to_remove]
  }
  number_of_sheets <- length(sheet_nums)
  # Add period to each prt table
  for (k in 1:length(all_prt)) {
    for (j in 1:(length(all_prt[[k]]$x)-1)){
      all_prt[[k]]$period[j] = 4* (all_prt[[k]]$x[j+1]-all_prt[[k]]$x[j])
    }
    all_prt[[k]]$period[length(all_prt[[k]]$period)] <- NA
    
  }

  # print all_prt to output_file2
  write.xlsx2(all_prt[[1]], output_file4, sheetName = as.character(sheet_nums[1]), row.names=FALSE)
  for (i in 2:number_of_sheets) {
    df <- all_prt[[i]] 
    write.xlsx2(df, output_file4, sheetName = as.character(sheet_nums[i]), row.names=FALSE, append = TRUE)
  }
  
  #calculating phase shift
  phase_avg_for_prism <- data.frame("blank" = 0)
  for (i in 2:number_of_sheets) {
    phase_avg_for_prism[,i-1] <- PSHIFT(all_prt[[1]], all_prt[[i]]) %>%  mean()
  }
  colnames(phase_avg_for_prism) <- c(as.character(sheet_nums[-1]))
  write.xlsx2(phase_avg_for_prism, output_file1, sheetName = "Phase_avg_Data", row.names=FALSE, append = TRUE)

  if (file.exists("All_Phase_Data.xlsx") == 1){ file.remove("All_Phase_Data.xlsx") }

  for (i in 2:number_of_sheets) {
    temp_vector <- PSHIFT(all_prt[[1]], all_prt[[i]])
    write.xlsx2(temp_vector, output_file3, sheetName = as.character(sheet_nums[i]), row.names=FALSE, append = TRUE)
  }

  # calculate period
  period_avg_for_prism <- data.frame("blank" = 0)
  
  for (i in 1:number_of_sheets) {
    period_avg_for_prism[,i] <- all_prt[[i]][,"period"] %>% .[complete.cases(.)] %>% mean()
  }
  colnames(period_avg_for_prism) <- c(as.character(sheet_nums))
  write.xlsx2(period_avg_for_prism, output_file1, sheetName = "Period_avg_Data", row.names=FALSE, append = TRUE)

  #calculate average amplitude (total abs average, peak average, trough average)
  all_amp <- AMP_batch(all_prt)

  amp <- data.frame("X" = all_amp[[1]][,"names"])
  for (i in 1:number_of_sheets) {
    amp[,i+1] <- all_amp[[i]][,"amp"]
  }
  colnames(amp) <- c("Variables", as.character(sheet_nums))
  write.xlsx2(amp, output_file1, sheetName = "All_amplitude_data", row.names=FALSE, append = TRUE)
}

#### Unused Functions ####

# WRITE_DATA_FILES_old <- function(file_name, sheet_nums, start_value, end_value, output_file1, output_file2, output_file3){
#   
#   number_of_sheets <- length(sheet_nums)
#   file_names <- rep(file_name, number_of_sheets)
#   #
#   starts <- rep(start_value, number_of_sheets)
#   #
#   ends <- rep(end_value, number_of_sheets)
#   steps <- rep(0.01, number_of_sheets)
#   
#   #### Read data from excel into data frames ####
#   all_cond <- read_kronos_batch(file_names, sheet_nums, starts, ends, steps)
#   
#   lum_norm_all <- data.frame("X" = all_cond[[1]][,"x"])
#   for (i in 1:number_of_sheets) {
#     lum_norm_all[,i+1] <- all_cond[[i]][,"lum_norm"]
#   }
#   colnames(lum_norm_all) <- c("X", as.character(sheet_nums))
#   write.xlsx2(lum_norm_all, output_file1, sheetName = "Norm_Lum_data", row.names=FALSE)
#   
#   raw_lum_all <- data.frame("X" = all_cond[[1]][,"x"])
#   for (i in 1:number_of_sheets) {
#     raw_lum_all[,i+1] <- all_cond[[i]][,"y"]
#   }
#   colnames(raw_lum_all) <- c("X", as.character(sheet_nums))
#   write.xlsx2(raw_lum_all, output_file1, sheetName = "Raw_Lum_data", row.names=FALSE, append = TRUE)
#   
#   #calculate peaks, roots, troughs
#   all_prt <- PRT_batch(all_cond)
#   
#   period_avg_for_prism <- data.frame("blank" = 0)
#   for (i in 1:number_of_sheets) {
#     period_avg_for_prism[,i] <- all_prt[[i]][,"period"] %>% .[complete.cases(.)] %>% mean()
#   }
#   colnames(period_avg_for_prism) <- c(as.character(sheet_nums))
#   write.xlsx2(period_avg_for_prism, output_file1, sheetName = "Period_avg_Data", row.names=FALSE, append = TRUE)
#   
#   
#   write.xlsx2(all_prt[[1]][,c("x","lum_norm","period")], output_file2, 
#               sheetName = as.character(sheet_nums[1]), row.names=FALSE)
#   for (i in 2:number_of_sheets) {
#     df <- all_prt[[i]][,c("x","lum_norm","period")]
#     write.xlsx2(df, output_file2, sheetName = as.character(sheet_nums[i]), row.names=FALSE, append = TRUE)
#   }
#   
#   #calculate average amplitude (total abs average, peak average, trough average)
#   all_amp <- AMP_batch(all_cond)
#   
#   amp <- data.frame("X" = all_amp[[1]][,"names"])
#   for (i in 1:number_of_sheets) {
#     amp[,i+1] <- all_amp[[i]][,"amp"]
#   }
#   colnames(amp) <- c("Variables", as.character(sheet_nums))
#   write.xlsx2(amp, output_file1, sheetName = "All_amplitude_data", row.names=FALSE, append = TRUE)
#   
#   #calculating phase shift
#   return(all_prt)
#   # phase_avg_for_prism <- data.frame("blank" = 0)
#   # for (i in 2:number_of_sheets) {
#   #   phase_avg_for_prism[,i-1] <- PSHIFT(all_prt[[1]], all_prt[[i]]) %>%  mean()
#   # }
#   # colnames(phase_avg_for_prism) <- c(as.character(sheet_nums[-1]))
#   # write.xlsx2(phase_avg_for_prism, output_file1, sheetName = "Phase_avg_Data", row.names=FALSE, append = TRUE)
#   # 
#   # 
#   # if (file.exists("All_Phase_Data.xlsx") == 1){ file.remove("All_Phase_Data.xlsx") }
#   # 
#   # for (i in 2:number_of_sheets) {
#   #   temp_vector <- PSHIFT(all_prt[[1]], all_prt[[i]])
#   #   write.xlsx2(temp_vector, output_file3, sheetName = as.character(sheet_nums[i]), row.names=FALSE, append = TRUE)
#   # }
#   
# }
# 
# shift_axis <- function(p, y=0){
#   g <- ggplotGrob(p)
#   dummy <- data.frame(y=y)
#   ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
#   p + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))), 
#                         ymax=y, ymin=y) +
#     geom_hline(aes(yintercept=y), data = dummy) +
#     theme(axis.text.x = element_blank(), 
#           axis.ticks.x=element_blank())
#   
# }
# 
# lm_eqn <- function(df_1){
#   m <- lm(abs_lum_norm ~ x, df_1);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(coef(m)[1], digits = 2), 
#                         b = format(coef(m)[2], digits = 2), 
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));                 
# }
# 
# AMP_old <- function(df_norm){
#   #pull out peaks and troughs in data
#   peaks <- findPeaks(df_norm$lum_norm)-1
#   peaks_data <- do.call("rbind", list(df_norm[peaks, ]))
#   peak_amplitude <- ((sum(peaks_data$lum_norm)))/length(peaks_data$lum_norm)
#   
#   # peaks_df <- df_A_norm[peaks, ]
#   
#   troughs <- findValleys(df_norm$lum_norm)-1
#   troughs_data <- do.call("rbind", list(df_norm[troughs, ]))
#   trough_amplitude <- abs(((sum(troughs_data$lum_norm)))/length(troughs_data$lum_norm))
#   
#   #calculate amplitude
#   pt_data <- do.call("rbind", list(df_norm[peaks, ], df_norm[troughs, ]))
#   all_amplitude <- ((sum(abs(pt_data$lum_norm)))/length(pt_data$lum_norm))
#   
#   #Write values into new dataframe
#   amplitude_data <- data.frame(amp = c(all_amplitude, peak_amplitude, trough_amplitude))
#   
#   #add a column for the names
#   amplitude_data$names <- c("Ave Amplitude", "Ave Peak Amp", "Ave Trough Amp")
#   
#   return(amplitude_data)
# }
# 
# 
# AMP_batch_old <- function(all_cond){
#   all_AMP <- list()
#   for(j in 1:length(all_cond)){
#     AMP_single <- AMP(all_cond[[j]])
#     all_AMP[[j]] <- AMP_single
#   }
#   return(all_AMP)
# }