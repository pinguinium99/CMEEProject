criteria_d <- function(data,multiplier){
  fulldf <- data.frame(rep(0,((nrow(data))-3)))
  for(i in 1:(ncol(data)-1)){
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x <- x * multiplier 
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



criteria_c <- function(data,multiplier){
  fulldf <- data.frame(rep(0,((nrow(data))-3)))
  for(i in 1:(ncol(data)-1)){
    status <- c()
    for(j in 1:(nrow(data)-3)){
      x <- data[j,i]
      x <- x * multiplier
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