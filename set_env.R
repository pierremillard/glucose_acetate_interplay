###################
# Set environment #
###################

# we assume to be already in the working directory, otherwise you can update the path of the supplementary data folder here

wd <- getwd()
#wd <- "C:/Users/millard/Documents/GIT/glucose_acetate_interplay/glucose_acetate_interplay/"

# set model, data and results directory

model_dir <- file.path(wd, "model")
results_dir <- file.path(wd, "results")
data_dir <- file.path(wd, "data")

# go to working directory

setwd(wd)

# load misc.R

source("misc.R")
