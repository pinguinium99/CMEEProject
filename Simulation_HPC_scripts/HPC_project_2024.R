rm(list=ls())

name <- "Prannoy Thazhakaden"
preferred_name <- "Prannoy"
email <- "PVT23@imperial.ac.uk"
username <- "PVT23"

source("Project_2024_main_initial.R")
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) 
set.seed(iter)

# splitting into different speciation rates use equation
# use a suitable speciation rate for then managing the multiplier this is fairly clear on what will happen
# run multiple times
# 
end_mat <- paste0("sim",iter,".rda")
speciation <- c(0.0001, 0.00019, 0.00057, 0.00086)
community <- 250000
filename <- paste0("HPCoutput_file",iter,".csv")
start_matrix <- matrix_start(community,250000,500)
multigensp(600,start_matrix,10,speciation[iter %% 4+1],10,filename, end_mat)
