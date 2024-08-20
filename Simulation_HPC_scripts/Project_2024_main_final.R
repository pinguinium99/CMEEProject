name <- "Prannoy Thazhakaden"
preferred_name <- "Prannoy"
email <- "PVT23@imperial.ac.uk"
username <- "PVT23"


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)

matrix_start <- function(nspecies,population,row) {
  my_matrix <- matrix(sample( 
    1 : nspecies, population, replace = TRUE), nrow = row) # should i use sqrt for rows?
  return(my_matrix)
}

speciation_fun <- function(pro_mat,matrix,speciation_rate,deadx,deady){
  
  z <- max(unique(matrix),unique(pro_mat))
  pro_mat[deadx,deady]= z+1
  return(pro_mat)
}

pro_speciation <- function(mat, speciation_rate){
  protracted_matrix <- matrix(0,ncol = sqrt(length(mat)),nrow = sqrt(length(mat)))
  return(protracted_matrix)
}

zerolist <- function(time,n_sp){
  x <- list()
  for (i in 1:n_sp){
    x[i] <- list(rep(0,time))
  } 
  return(x)
}replacfun <- function(mat,protracted_mat,sd,speciation_rate){    
  deadx <- round(runif(1,1,sqrt(length(mat))))
  deady <- round(runif(1,1,sqrt(length(mat))))
  spec <- runif(1,0,1)
  x <- round(rnorm(1, mean = 0, sd = sd))

  y <- round(rnorm(1, mean = 0, sd = sd))

  if(sqrt(length(mat)) <= (deadx + x) | (deadx + x) <= 1 | sqrt(length(mat)) <= (deady + y) | (deady + y) <= 1 ){
    
    while(sqrt(length(mat)) <= (deadx + x) | (deadx + x) <= 1 | sqrt(length(mat)) <= (deady + y) | (deady + y) <= 1 ){
      x <- round(rnorm(1, mean = 0, sd = sd))
      y <- round(rnorm(1, mean = 0, sd = sd))
    }
    mat[deadx,deady] <- mat[deadx + x,deady + y]
    protracted_mat[deadx,deady] <- protracted_mat[deadx + x,deady + y]
    
  }else{
    mat[deadx,deady] <- mat[deadx + x,deady + y]
    protracted_mat[deadx,deady] <- protracted_mat[deadx + x,deady + y]
  }
  if(spec < speciation_rate){
    
    protracted_mat <- speciation_fun(protracted_mat,mat,0.03,deadx,deady)
  }
  
  output <- list(mat,protracted_mat)
  return(output)
}

update_matrix <- function(matrix1, matrix2, value) {
  positions <- which(matrix1 == value, arr.ind = TRUE)
  for (pos in 1:nrow(positions)) {
    row <- positions[pos, 1]
    col <- positions[pos, 2]
    matrix2[row, col] <- value
    matrix1[row,col] <- 0
  }
  return(list(matrix2,matrix1))
}

multigensp <- function(time,mat,sd,speciation_rate,speciation_length,file_name){
  pro_mat <- pro_speciation(mat)
  new_count <- vector()
  x <- zerolist(time,length(unique(as.vector(mat))))
  first <- table(mat)
  sp1 <- sort(unique(as.vector(mat)))
  for(k in sp1){
    name <- names(first)
    index <- which(name == k)
    x[[k]][1]<- first[index]
  }
  for(i in 1:time){
    for(j in 1:(length(mat)/2)){
      result <- replacfun(mat,pro_mat,sd,speciation_rate)
      mat <- result[[1]]
      pro_mat <- result[[2]]
    }
    new_sp <- unique(pro_mat)
    new_count <- append(new_count,new_sp[new_sp != 0])
    species_check <- table(new_count)
    
    if(any(species_check >= speciation_length)){
      species_surv <- as.numeric(names(species_check[species_check >= speciation_length]))
      for (new in species_surv){
        newly_speciated <- update_matrix(pro_mat, mat, new)
        mat <- newly_speciated[[1]]
        pro_mat <- newly_speciated[[2]]
        new_count <- new_count[new_count != new]
        print(new_count)
      }
    }
    y <- table(mat)
    sp <- sort(unique(as.vector(mat)))
    if (max(sp) > max(sp1)){
      newsp <- (max(sp) - max(sp1))
      x <- append(x,(zerolist(time,newsp)))
    }
    for(k in sp){
      name <- names(y)
      index <- which(name == k)
      x[[k]][i]<- y[index]
    }
    sp1 <- sp
  }
  write.csv(x,file = file_name)
  
}

dummydat <- read.csv("HPCoutput_file1.csv", header = TRUE)
my_df <- dummydat
name <- vector()
for (i in 1:ncol(my_df)){
  name[i] <- paste0("species",i)
} 
t <- (1:nrow(my_df))
my_df <- cbind(my_df,t)
name[length(name)+1] <- "time"
colnames(my_df) <- name

rm_sp <- function(df){
  df <- df[, sapply(df, function(col) sum(col > 0)) > 3] 
  return(df)
}

my_df<-rm_sp(my_df)

for (i in 1:(ncol(my_df)-1)) {
  p  <- ggplot(data = my_df, aes(x = time, y = my_df[,i]))+
    geom_line()
  print(p)
}

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
heat_map_pre <- rate_of_change(worst_status(my_df))

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