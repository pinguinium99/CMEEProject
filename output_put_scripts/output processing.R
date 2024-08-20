

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
options(dplyr.summarise.inform = FALSE)
options(dplyr.group_by.inform = FALSE)

  rm_sp <- function(df){
    df <- df[, sapply(df, function(col) sum(col > 0)) > 3] 
    return(df)
  }


#for (i in 1:(ncol(my_df)-1)) {
#  p  <- ggplot(data = my_df, aes(x = time, y = my_df[,i]))+
#    geom_line()
#  print(p)
#}

rm_0 <- function(sp){
  df_new <- sp[sp!=0]
  print(df_new)
}

criteria_a <- function(data){
  fulldf <- data.frame(rep(0,((nrow(data))-3)))
  for(i in 1:(ncol(data)-1)){
    abundance_change <- c()
    percentage_change <- c()
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x_plus_3 <- data[(j+3),i]
      abundance_change[j] <- x - x_plus_3
      if(x == x_plus_3){
        percentage_dif <- 0
      } else{
        percentage_dif <- (((x - x_plus_3)/x)*100)
      }
      percentage_change[j] <- percentage_dif
      if (x == 0){
        status [j]<- 0
      }else if (percentage_dif <= -90){
        status [j] <- 1
      } else if (percentage_dif <= -70){
        status [j] <- 2
      } else if (percentage_dif <= -50){
        status [j] <- 3
      } else{
        status[j]<- 4 
      }
    }
    sp <- c(status)
    fulldf[,i] <- sp
  }
  return(fulldf)
}


criteria_c <- function(data){
  fulldf <- data.frame(rep(0,((nrow(data))-3)))
  for(i in 1:(ncol(data)-1)){
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x <- x * 25
      one_gen <- (((x - (data[j+1,i]))*100))
      two_gen <- (((x - (data[j+2,i]))*100))
      three_gen <- (((x - (data[j+3,i]))*100))
      if (x == 0){
        status [j]<- 0 
      }else if (x <= 250 & one_gen < -25){
        # print("critical")
        status [j] <- 1
      } else if (x <= 2500 & two_gen < -20){
        status [j] <- 2 
      } else if (x <= 10000 & three_gen < -10){
        status [j] <- 3
      }else{
        status[j]<- 4
      }
    }
    sp <- c(status)
    fulldf[,i] <- sp
  }
  return(fulldf)
}

criteria_d <- function(data){
  fulldf <- data.frame(rep(0,((nrow(data))-3)))
  for(i in 1:(ncol(data)-1)){
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x <- x * 25 
      if (x == 0){
        status [j]<- 0 
      } else if (x <= 50){
        status [j] <- 1
      } else if (x <= 250){
        status [j] <- 2
      } else if (x <= 1000){
        status [j] <- 3
      } else{
        status[j] <- 4
      }
    }
    sp <- c(status)
    fulldf[,i] <- sp
  }
  return(fulldf)
}

worst_status <- function(data){
  dat_a <- criteria_a(data)
  dat_c <- criteria_c(data)
  dat_d <- criteria_d(data)
  init <- data.frame(rep(0,((nrow(data))-3)))
  for (i in 1:length(dat_a)){
    x <- cbind(dat_a[i],dat_c[i],dat_d[i])
    min_values <- apply(x, 1, min)
    new_df <- data.frame(minimum_value = min_values)
    init <- cbind(init,new_df)
  }
  return(init[,-1])
}

output <- c(14,15,16)

for(i in output){
  
  filename <- paste0("HPCoutput_file", i, ".csv")
  my_df <- read.csv(filename, header = TRUE)


my_df <- (rm_sp(my_df))

name <- vector()
for (i in 1:ncol(my_df)){
  name[i] <- paste0("species",i)
  colnames(my_df) <- name
}
  
  criteria_a_results <- criteria_a(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_A", i, ".csv"))
  
  criteria_c_results <- criteria_c(my_df)
  fwrite(criteria_c_results, file = paste0("criteria_C", i, ".csv"))
  
  criteria_d_results <- criteria_d(my_df)
  fwrite(criteria_d_results, file = paste0("criteria_D", i, ".csv"))
  
  rm(criteria_a_results,criteria_c_results,criteria_d_results)
  
  worst_status_results <- worst_status(my_df)
  fwrite(worst_status_results, file = paste0("worst_status", i, ".csv"))
  
  
  gc()
  

}





rate_of_change <- function(dat){
  hm_dat <- data.frame(previous_status = character(),
                       status = character(),
                       count = numeric(),
                       stringsAsFactors = FALSE) 
  for(i in 1:ncol(dat)){
    sp_stat <- data.frame(dat[,i])
    colnames(sp_stat) <- "status"
    status_transition_counts <- sp_stat %>%
      group_by(previous_status = lag(status, default = first(status)), status) %>%
      summarize(count = n()) %>%
      filter(!is.na(previous_status))
    print(status_transition_counts)
    hm_dat <- rbind(hm_dat,status_transition_counts)
  }
  return(hm_dat)
}

heat_map_pre <- rate_of_change(worst_status(rm_sp(my_df)))

status_list <- c("extinct","critical","endangered","vulnrable","least_concern")
status_map <- setNames(status_list, 0:4)


heat_map_pre <- heat_map_pre %>% group_by(previous_status,status) %>% summarise(count = sum(count)) %>% filter(previous_status != 0)

heat_map_pre <- heat_map_pre %>%
  mutate(
    previous_status = status_map[as.character(previous_status)],
    status = status_map[as.character(status)]
  )

heat_map_pre <- heat_map_pre %>%
  mutate(
    previous_status = factor(previous_status, levels = status_list),
    status = factor(status, levels = status_list)
  )
total_count <- sum(heat_map_pre$count)
heat_map_pre <- heat_map_pre %>%
  mutate(count = count / total_count)

ggplot(heat_map_pre, aes(previous_status, status, fill= count)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlOrRd", trans = 'log',na.value = "black") 

ggplot(heat_map_pre, aes(previous_status, status, fill = count)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Greys" ,trans = 'log') +
  labs(title = "Rate of Status Transitions", x = "Previous Status", y = "Current Status")+
  theme_classic()



#------------------ Full run ---------------------------------------
# use gc and rm and large files/objects 



library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
options(dplyr.summarise.inform = FALSE)
options(dplyr.group_by.inform = FALSE)


rm_sp <- function(df){
  df <- df[, sapply(df, function(col) sum(col > 0)) > 3] 
  return(df)
}

rm_0 <- function(sp){
  df_new <- sp[sp!=0]
  print(df_new)
}

#for (i in 1:(ncol(my_df)-1)) {
#  p  <- ggplot(data = my_df, aes(x = time, y = my_df[,i]))+
#    geom_line()
#  print(p)
#}




#t <- (1:nrow(my_df))
#my_df <- cbind(my_df,t)
#name[length(name)+1] <- "time"

criteria_a <- function(data){
  fulldf <- data.table(matrix(0, nrow = nrow(data) - 3, ncol = ncol(data)))
  for(i in 1:(ncol(data))){
    abundance_change <- c()
    percentage_change <- c()
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x_plus_3 <- data[(j+3),i]
      abundance_change[j] <- x - x_plus_3
      if(x == x_plus_3){
        percentage_dif <- 0
      } else{
        percentage_dif <- (((x - x_plus_3)/x)*100)
      }
      percentage_change[j] <- percentage_dif
      if (x == 0){
        status [j]<- 0
      }else if (percentage_dif <= -90){
        status [j] <- 1
      } else if (percentage_dif <= -70){
        status [j] <- 2
      } else if (percentage_dif <= -50){
        status [j] <- 3
      } else{
        status[j]<- 4 
      }
    }
    sp <- c(status)
    fulldf[, paste0("species", i) := status]
  }
  return(fulldf)
}

criteria_c <- function(data, multiplier){
  fulldf <- data.table(matrix(0, nrow = nrow(data) - 3, ncol = ncol(data)))
  for(i in 1:(ncol(data)-1)){
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i] * multiplier
      one_gen <- (((x - (data[j+1,i]))*100))
      two_gen <- (((x - (data[j+2,i]))*100))
      three_gen <- (((x - (data[j+3,i]))*100))
      if (x == 0){
        status [j]<- 0 
      }else if (x <= 250 & one_gen < -25){
        status [j] <- 1
      } else if (x <= 2500 & two_gen < -20){
        status [j] <- 2 
      } else if (x <= 10000 & three_gen < -10){
        status [j] <- 3
      }else{
        status[j]<- 4
      }
    }
    sp <- c(status)
    fulldf[, paste0("species", i) := status]
  }
  return(fulldf)
}

criteria_d <- function(data,multiplier){
  fulldf <-  data.table(matrix(0, nrow = nrow(data) - 3, ncol = ncol(data)))
  for(i in 1:(ncol(data)-1)){
    status <- numeric(nrow(data) - 3)
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i] * multiplier
      if (x == 0){
        status [j]<- 0 
      } else if (x <= 50){
        status [j] <- 1
      } else if (x <= 250){
        status [j] <- 2
      } else if (x <= 1000){
        status [j] <- 3
      } else{
        status[j] <- 4
      }
    }
    sp <- c(status)
    fulldf[,paste0("species", j) := status]
  }
  return(fulldf)
}

worst_status <- function(data,multiplier){
  dat_a <- criteria_a(data)
  dat_c <- criteria_c(data,multiplier)
  dat_d <- criteria_d(data,multiplier)
  init <- data.table(matrix(0, nrow = nrow(data) - 3, ncol = ncol(data)))
  for (i in 1:length(dat_a)){
    x <- cbind(dat_a[i],dat_c[i],dat_d[i])
    min_values <- apply(x, 1, min)
    init[, paste0("species", i) := min_values]
  }
  return(init[,-1])
}






#output <- c(1,2,3,4,5,6,8,9,10,11,12,13,14,15,16)
output <- c(1,2)

for (i in output){
  filename <- paste0("HPCoutput_file",i,".csv")
  my_df <- fread(filename)
  my_df <-rm_sp(my_df)
  
  
  colnames(my_df) <- paste0("species", 1:ncol(my_df))
  
  criteria_a_results <- criteria_a(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_A", i, ".csv"))
  
  criteria_c_results <- criteria_c(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_C", i, ".csv"))
  
  criteria_d_results <- criteria_d(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_D", i, ".csv"))
  

  worst_status_results <- worst_status(my_df)
  fwrite(worst_status_results, file = paste0("worst_status", i, ".csv"))
  

  gc()

}


library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

# Remove sparse species
rm_sp <- function(df) {
  df <- df[, which(colSums(df > 0) > 3), with = FALSE]
  return(df)
}

# Criteria functions
criteria_a <- function(data) {
  fulldf <- matrix(0, nrow = nrow(data) - 3, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(nrow(data) - 3)) {
      x <- data[j, i, with = FALSE][[1]]
      x_plus_3 <- data[j + 3, i, with = FALSE][[1]]
      percentage_change <- if (x == x_plus_3) 0 else ((x - x_plus_3) / x) * 100
      fulldf[j, i] <- if (x == 0) 0 else if (percentage_change <= -90) 1 else if (percentage_change <= -70) 2 else if (percentage_change <= -50) 3 else 4
    }
  }
  return(as.data.frame(fulldf))
}

criteria_c <- function(data) {
  fulldf <- matrix(0, nrow = nrow(data) - 3, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(nrow(data) - 3)) {
      x <- data[j, i, with = FALSE][[1]] * 25
      one_gen <- ((x - data[j + 1, i, with = FALSE][[1]]) * 100)
      two_gen <- ((x - data[j + 2, i, with = FALSE][[1]]) * 100)
      three_gen <- ((x - data[j + 3, i, with = FALSE][[1]]) * 100)
      fulldf[j, i] <- if (x == 0) 0 else if (x <= 250 & one_gen < -25) 1 else if (x <= 2500 & two_gen < -20) 2 else if (x <= 10000 & three_gen < -10) 3 else 4
    }
  }
  return(as.data.frame(fulldf))
}

criteria_d <- function(data) {
  fulldf <- matrix(0, nrow = nrow(data) - 3, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(nrow(data) - 3)) {
      x <- data[j, i, with = FALSE][[1]] * 25
      fulldf[j, i] <- if (x == 0) 0 else if (x <= 50) 1 else if (x <= 250) 2 else if (x <= 1000) 3 else 4
    }
  }
  return(as.data.frame(fulldf))
}

# Combining results from criteria
worst_status <- function(data) {
  dat_a <- criteria_a(data)
  dat_c <- criteria_c(data)
  dat_d <- criteria_d(data)
  init <- matrix(0, nrow = nrow(data) - 3, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    x <- cbind(dat_a[, i], dat_c[, i], dat_d[, i])
    min_values <- apply(x, 1, min)
    init[, i] <- min_values
  }
  return(as.data.frame(init))
}

# Processing files
output <- c(1, 2)
for (i in output) {
  filename <- paste0("HPCoutput_file", i, ".csv")
  my_df <- fread(filename)
  
  # Remove sparse species
  my_df <- rm_sp(my_df)
  
  # Assign column names
  setnames(my_df, paste0("species", seq_len(ncol(my_df))))
  
  # Save results for criteria A
  criteria_a_results <- criteria_a(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_A", i, ".csv"))
  
  criteria_c_results <- criteria_c(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_C", i, ".csv"))
  
  criteria_d_results <- criteria_d(my_df)
  fwrite(criteria_a_results, file = paste0("criteria_D", i, ".csv"))
  
  
  
  # Save worst status results
  worst_status_results <- worst_status(my_df)
  fwrite(worst_status_results, file = paste0("worst_status", i, ".csv"))
  
  # Garbage collection
  gc()
}


library(dplyr)
library(ggplot2)
par(mfrow=c(2,2))
speciation <- c(0.0001, 0.00019, 0.00057, 0.00086)
setspish <- c(1,2,3,4,5,6,8,9,10,11,12,13,14,15,16)
order <- setspish[setspish %% 4+1]
datavec <- list()
for(i in 1:4){
  indices <- which(order == i)
  datavec[[i]] <- setspish[indices]
}
datavec
i = 1
for (listval in datavec){
  x <- listval
  empty_df <- data.frame(nrow = 597)
  for(sp in x){
    file_name <- paste0("./Data/HPCoutput_file/worst_status",sp,".csv")
    data <- read.csv(file_name)
    data <- rm_sp(data)
    data <- data[,-1:-100]
    empty_df <- cbind(empty_df, data)
  }
  heat_map_pre <- rate_of_change(empty_df)
  status_list <- c("extinct","critical","endangered","vulnrable","least_concern")
  status_map <- setNames(status_list, 0:4)
  heat_map_pre <- heat_map_pre %>% group_by(previous_status,status) %>% summarise(count = sum(count)) %>% filter(previous_status != 0)
  heat_map_pre <- heat_map_pre %>%
    mutate(
      previous_status = status_map[as.character(previous_status)],
      status = status_map[as.character(status)]
    )
  heat_map_pre <- heat_map_pre %>%
    mutate(
      previous_status = factor(previous_status, levels = status_list),
      status = factor(status, levels = status_list)
    )
  total_count <- sum(heat_map_pre$count)
  heat_map_pre <- heat_map_pre %>%
    mutate(count = count / total_count)
  write.csv(heat_map_pre,file = paste0("speciation_rate",speciation[i]))
  i <- i+1
}
