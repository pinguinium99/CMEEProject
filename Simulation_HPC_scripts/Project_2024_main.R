name <- "Prannoy Thazhakaden"
preferred_name <- "Prannoy"
email <- "PVT23@imperial.ac.uk"
username <- "PVT23"

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
}

replacfun <- function(mat,protracted_mat,sd,speciation_rate){    
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

multigensp <- function(time,mat,sd,speciation_rate,speciation_length,file_name,filename_mat){
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
  save(mat,pro_mat,file = filename_mat)
}



