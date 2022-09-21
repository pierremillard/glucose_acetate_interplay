##########
# README #
##########

# This code performs the computational analyses detailed in:
#
#   From toxic waste to beneficial nutrient: acetate boosts Escherichia coli growth at low glycolytic flux.
#   P. Millard, S. Uttenweiler-Joseph, B. Enjalbert
#   bioRxiv preprint, doi: 10.1101/2022.09.20.506926
#
# All scripts are available at https://github.com/MetaSys-LISBP/glucose_acetate_interplay
#
# Author: Pierre Millard (pierre.millard@insa-toulouse.fr)
# Copyright: INRAE, 2022


##########################
# INITIALIZE ENVIRONMENT #
##########################

# go to code directory (uncomment the following line)

#setwd("~/GIT/glucose_acetate_interplay/glucose_acetate_interplay/")

# load libraries and initialize environment

source("set_env.R")


################################################################################################
# SIMULATE THE RESPONSE OF E. COLI TO CHANGES OF GLYCOLYTIC ACTIVITY AND ACETATE CONCENTRATION #
################################################################################################
# Outputs:                                                                                     #
#   - Data_Figure_1.RData: Steady-state fluxes (glucose uptake, growth rate, acetate flux)     #
#                          simulated for all perturbations                                     #
#   - Figure_1B.pdf:       Fig. 1, panel B                                                     #
#   - Figure_1C.pdf:       Fig. 1, panel C                                                     #
#   - Figure_1D.pdf:       Fig. 1, panel D                                                     #
################################################################################################

# load the kinetic model detailed at doi: 10.7554/eLife.63661 and
# available from Biomodels database (https://www.ebi.ac.uk/biomodels/) with the identifier MODEL2005050001. 

setwd(model_dir)
loadModel("Millard2020_Ecoli_glc_ace_kinetic_model.cps")

# remove events and define concentrations of acetate, glucose and biomass as fixed

deleteEvent(getEvents()$key)
setSpecies(key="Ace_out", type="fixed")
setSpecies(key="Glc", type="fixed")
setSpecies(key="X", type="fixed")

# definition of the perturbations to simulate

#   number of steps (50)
n_step <- 50
#   acetate concentration (between 0.1 and 100 mM)
ace_concentration <- 10**seq(-1, 2, length.out = n_step)
#   glycolytic activity (between 20 and 100 % of initial Vmax)
glc_vmax <- getParameters(key="(glc_upt).V")$value*seq(0.2*n_step, n_step)/n_step

# create empty matrix to store simulation results

fluxes <- c("Values[v_growth_rate]", "Values[v_glc_uptake]", "Values[v_ace_net]")
simulation_results <- array(NA, dim=c(n_step*0.8+1, n_step, length(fluxes)+1), dimnames=list(inhibition=NULL, ace_concentration=NULL, specie=c("ace_concentration", fluxes)))

# run simulations

# for each glycolytic activity
for (k in seq(n_step*0.8+1)){
  
  # create empty matrix to store simulation results for a given glycolytic activity (k)
  res_ace_range <- matrix(NA, nrow=n_step, ncol=length(fluxes)+1, dimnames=list(r=NULL, c=c("ace_concentration", fluxes)))
  
  # update Vmax_glycolysis
  setParameters(key="(glc_upt).V", value = glc_vmax[k])
  applyInitialState()
  
  # for each acetate concentration
  for (j in seq(n_step)){
    
    # update acetate concentration
    setSpecies(key="Ace_out", initial_concentration=ace_concentration[j])
    applyInitialState()
    
    # calculate steady-state fluxes
    res_ss <- runSteadyState()
    
    # save steady-state fluxes for glycolytic activity k & acetate concentration j
    res_ace_range[j,] <- c(ace_concentration[j], unlist(res_ss$global_quantities[res_ss$global_quantities$key %in% fluxes, "value"]))
  }
  
  # save steady-state fluxes for all acetate concentrations at glycolytic activity k
  simulation_results[k,,] <- res_ace_range
  
}

# go to results directory

setwd(results_dir)

# save simulation results, can be loaded in R with:
#   simulation_results <- load("Data_Figure_1.RData")

save(simulation_results, file="Data_Figure_1.RData")

# plot simulation results (Figure 1)

pdf(file="Figure_1B.pdf", width = 5, height = 4)

  filled.contour(x=log10(ace_concentration), y=(glc_vmax/max(glc_vmax))*100, z = t(simulation_results[,,2]), nlevels = 10,
                 plot.axes = { axis(1, at=c(-1, 0, 1, 2), label=c(0.1, 1, 10, 100));
                 axis(2)})
dev.off()
  
  
pdf(file="Figure_1C.pdf", width = 5, height = 4)
  filled.contour(x=log10(ace_concentration), y=(glc_vmax/max(glc_vmax))*100, z = t(simulation_results[,,3]), nlevels = 10,
                 plot.axes = { axis(1, at=c(-1, 0, 1, 2), label=c(0.1, 1, 10, 100));
                 axis(2)})
dev.off()

pdf(file="Figure_1D.pdf", width = 5, height = 4)
  cols = colorRampPalette(c('royalblue4','lightblue'))
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

  filled.contour(x=log10(ace_concentration), y=(glc_vmax/max(glc_vmax))*100, z = t(simulation_results[,,4]), nlevels = 10, 
                 col = c(cols(5), color.palette(4)),
                 plot.axes = { axis(1, at=c(-1, 0, 1, 2), label=c(0.1, 1, 10, 100));
                 axis(2)})
dev.off()


##############################################################################################
# GLYCOLYTIC CONTROL OF THE ACETATE FLUX, IN ABSENCE OR PRESENCE OF ACETATE (30 mM)          #
##############################################################################################
# Outputs:                                                                                   #
#   - Data_Figure_2.txt: All data, including the total carbon uptake flux (qGLC+qACE)        #
#                        and the contribution of acetate to carbon uptake (qACE/(qGLC+qACE)) #
#   - Figure_2.pdf:      Fig. 2                                                              #
##############################################################################################

# load experimental data

setwd(data_dir)
exp_data <- read.table("data_E_coli_K12_MG1655_Figure_2.txt", header=TRUE, sep="\t")

# go to results directory

setwd(results_dir)

# calculate the total carbon uptake flux (qGLC+qACE)

Ctotal <- -exp_data$qACE*2 + exp_data$qGLC*6

# calculate the contribution of acetate to carbon uptake (qACE/(qGLC+qACE))

contribution_acetate <- -exp_data$qACE*2 / (-exp_data$qACE*2 + exp_data$qGLC*6) * 100

# add calculated metrics to measurements dataframe

exp_data <- cbind(exp_data, Ctotal, contribution_acetate)

# save all data, can be loaded in R with exp_data <- read.table(file="Data_Figure_2.txt", sep="\t")

write.table(exp_data, file="Data_Figure_2.txt", quote=FALSE, sep="\t")

# plot all data as function of aMG concentration

pdf(file = "Figure_2.pdf", width = 6.3, height = 8.1)
  
  par(mfrow=c(4,3))
  
  filter <- (exp_data$ACE<2)
  
  # Glucose uptake flux as function of aMG concentration
  # fit polynomial model
  polyfit_no_ace <- polynomial_model(exp_data[filter,], "aMG", "qGLC", 2, lpred=seq(0,100,by=0.1))
  polyfit_ace <- polynomial_model(exp_data[!filter,], "aMG", "qGLC", 2, lpred=seq(0,100,by=0.1))
  # plot fit and confidence intervals
  plot(polyfit_no_ace$fited_curve_x, polyfit_no_ace$fited_curve_y, ylim=c(0,11), las=1, yaxs="i", ylab="qGLC (mmol/gDW/h)", xlab="[aMG] (mM)", col="#1F4E79", lwd=1.5, type="l")
  lines(polyfit_ace$fited_curve_x, polyfit_ace$fited_curve_y, col="#9DC3E6", lwd=1.5)
  # NOTE: uncomment the two following lines to plot 95 % confidence intervals on the polynomial fits
  #polygon(c(polyfit_no_ace$fited_curve_x, rev(polyfit_no_ace$fited_curve_x)), c(polyfit_no_ace$confidence_intervals[,1], rev(polyfit_no_ace$confidence_intervals[,2])), col=adjustcolor("#1F4E79", alpha.f=0.4), border = NA)
  #polygon(c(polyfit_ace$fited_curve_x, rev(polyfit_ace$fited_curve_x)), c(polyfit_ace$confidence_intervals[,1], rev(polyfit_ace$confidence_intervals[,2])), col=adjustcolor("#9DC3E6", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$aMG[filter], exp_data$qGLC[filter], lwd=1.0, pch=21, bg="#1F4E79", col="black", cex=1)
  points(exp_data$aMG[!filter], exp_data$qGLC[!filter], pch=21, bg="#9DC3E6", col="black", cex=1, lwd=1.0)
  
  # Acetate flux as function of aMG concentration
  # fit polynomial model
  polyfit_no_ace <- polynomial_model(exp_data[filter,], "aMG", "qACE", 2, lpred=seq(0,100,by=0.1))
  polyfit_ace <- polynomial_model(exp_data[!filter,], "aMG", "qACE", 2, lpred=seq(0,100,by=0.1))
  # plot fit and confidence intervals
  plot(polyfit_no_ace$fited_curve_x, polyfit_no_ace$fited_curve_y, ylim=c(-11,3), las=1, yaxs="i", ylab="qACE (mmol/gDW/h)", xlab="[aMG] (mM)", col="#548235", lwd=1.5, type="l")
  abline(h=0)
  lines(polyfit_ace$fited_curve_x, polyfit_ace$fited_curve_y, col="#A9D18E", lwd=1.5)
  #polygon(c(polyfit_no_ace$fited_curve_x, rev(polyfit_no_ace$fited_curve_x)), c(polyfit_no_ace$confidence_intervals[,1], rev(polyfit_no_ace$confidence_intervals[,2])), col=adjustcolor("#548235", alpha.f=0.4), border = NA)
  #polygon(c(polyfit_ace$fited_curve_x, rev(polyfit_ace$fited_curve_x)), c(polyfit_ace$confidence_intervals[,1], rev(polyfit_ace$confidence_intervals[,2])), col=adjustcolor("#A9D18E", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$aMG[filter], exp_data$qACE[filter], lwd=1.0, pch=21, bg="#548235", col="black", cex=1)
  points(exp_data$aMG[!filter], exp_data$qACE[!filter], pch=21, bg="#A9D18E", col="black", cex=1, lwd=1.0)
  
  # Growth rate as function of aMG concentration
  # fit polynomial model
  polyfit_no_ace <- polynomial_model(exp_data[filter,], "aMG", "mu", 2, lpred=seq(0,100,by=0.1))
  polyfit_ace <- polynomial_model(exp_data[!filter,], "aMG", "mu", 2, lpred=seq(0,100,by=0.1))
  # plot fit and confidence intervals
  plot(polyfit_no_ace$fited_curve_x, polyfit_no_ace$fited_curve_y, ylim=c(0,0.7), las=1, yaxs="i", ylab="growth rate (h-1)", xlab="[aMG] (mM)", col="#B63A2D", lwd=1.5, type="l")
  lines(polyfit_ace$fited_curve_x, polyfit_ace$fited_curve_y, col="#E49890", lwd=1.5)
  #polygon(c(polyfit_no_ace$fited_curve_x, rev(polyfit_no_ace$fited_curve_x)), c(polyfit_no_ace$confidence_intervals[,1], rev(polyfit_no_ace$confidence_intervals[,2])), col=adjustcolor("#B63A2D", alpha.f=0.4), border = NA)
  #polygon(c(polyfit_ace$fited_curve_x, rev(polyfit_ace$fited_curve_x)), c(polyfit_ace$confidence_intervals[,1], rev(polyfit_ace$confidence_intervals[,2])), col=adjustcolor("#E49890", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$aMG[filter], exp_data$mu[filter], lwd=1.0, pch=21, bg="#B63A2D", col="black", cex=1)
  points(exp_data$aMG[!filter], exp_data$mu[!filter], pch=21, bg="#E49890", col="black", cex=1, lwd=1.0)
  
  # Total carbon uptake as function of aMG concentration
  # fit polynomial model
  polyfit_no_ace <- polynomial_model(exp_data[filter,], "aMG", "Ctotal", 2, lpred=seq(0,100,by=0.1))
  polyfit_ace <- polynomial_model(exp_data[!filter,], "aMG", "Ctotal", 2, lpred=seq(0,100,by=0.1))
  # plot fit and confidence intervals
  plot(polyfit_no_ace$fited_curve_x, polyfit_no_ace$fited_curve_y, ylim=c(0, 60), las=1, yaxs="i", ylab="total carbon uptake (mmol/gDW/h)", xlab="[aMG] (mM)", col="#451E5D", lwd=1.5, type="l")
  lines(polyfit_ace$fited_curve_x, polyfit_ace$fited_curve_y, col="#C89CE4", lwd=1.5)
  #polygon(c(polyfit_no_ace$fited_curve_x, rev(polyfit_no_ace$fited_curve_x)), c(polyfit_no_ace$confidence_intervals[,1], rev(polyfit_no_ace$confidence_intervals[,2])), col=adjustcolor("#451E5D", alpha.f=0.4), border = NA)
  #polygon(c(polyfit_ace$fited_curve_x, rev(polyfit_ace$fited_curve_x)), c(polyfit_ace$confidence_intervals[,1], rev(polyfit_ace$confidence_intervals[,2])), col=adjustcolor("#C89CE4", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$aMG[filter], exp_data$Ctotal[filter], lwd=1.0, pch=21, bg="#451E5D", col="black", cex=1)
  points(exp_data$aMG[!filter], exp_data$Ctotal[!filter], pch=21, bg="#C89CE4", col="black", cex=1, lwd=1.0)
  
  # Acetate contribution to carbon uptake as function of aMG concentration
  # fit polynomial model
  polyfit <- polynomial_model(exp_data, "qGLC", "qACE", 2, lpred=seq(0,10,by=0.01))
  # note: to avoid a bug on windows platforms when the plotted polygon is bigger than the plot, see https://bugs.r-project.org/show_bug.cgi?id=18219
  #polyfit$confidence_intervals[polyfit$confidence_intervals > 3] <- 3
  #polyfit$confidence_intervals[polyfit$confidence_intervals < -11] <- -11
  # plot fit and confidence intervals
  plot(polyfit$fited_curve_x, polyfit$fited_curve_y, ylim=c(-11, 3), las=1, yaxs="i", ylab="qACE (mmol/gDW/h)", xlab="qGLC (mmol/gDW/h)", col="#BF9000", lwd=1.5, type="l")
  abline(h=0)
  #polygon(c(polyfit$fited_curve_x, rev(polyfit$fited_curve_x)), c(polyfit$confidence_intervals[,1], rev(polyfit$confidence_intervals[,2])), col=adjustcolor("#BF9000", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$qGLC[filter], exp_data$qACE[filter], lwd=1.0, pch=21, bg="#BF9000", col="black", cex=1)
  points(exp_data$qGLC[!filter], exp_data$qACE[!filter], lwd=1.0, pch=21, bg="#FFE699", col="black", cex=1)
  
  # Acetate contribution to carbon uptake as function of aMG concentration
  # fit polynomial model
  polyfit <- polynomial_model(exp_data, "qGLC", "contribution_acetate", 2, lpred=seq(0,10,by=0.01))
  # note: to avoid a bug on windows platforms when the plotted polygon is bigger than the plot, see https://bugs.r-project.org/show_bug.cgi?id=18219
  #polyfit$confidence_intervals[polyfit$confidence_intervals > 60] <- 60
  #polyfit$confidence_intervals[polyfit$confidence_intervals < -20] <- -20
  # plot fit and confidence intervals
  plot(polyfit$fited_curve_x, polyfit$fited_curve_y, ylim=c(-20, 60), las=1, yaxs="i", ylab="acetate contribution to C uptake (%)", xlab="qGLC (mmol/gDW/h)", col="#BF9000", lwd=1.5, type="l")
  abline(h=0)
  #polygon(c(polyfit$fited_curve_x, rev(polyfit$fited_curve_x)), c(polyfit$confidence_intervals[,1], rev(polyfit$confidence_intervals[,2])), col=adjustcolor("#BF9000", alpha.f=0.4), border = NA)
  # plot experimental data
  points(exp_data$qGLC[filter], exp_data$contribution_acetate[filter], lwd=1.0, pch=21, bg="#BF9000", col="black", cex=1)
  points(exp_data$qGLC[!filter], exp_data$contribution_acetate[!filter], lwd=1.0, pch=21, bg="#FFE699", col="black", cex=1)

dev.off()


###########################################################################################
# DETERMINATION OF THE ACETATE THRESHOLD CONCENTRATION AT WHICH THE ACETATE FLUX REVERSES #
###########################################################################################
# Output:                                                                                 #
#   - Figure_3.pdf:      Fig. 3                                                           #
###########################################################################################

# load experimental data

setwd(data_dir)
exp_data <- read.table("data_E_coli_K12_MG1655_Figure_3-4.txt", header=TRUE, sep="\t")

# go to results directory

setwd(results_dir)

# plot acetate flux as function of acetate concentration, at different glycolytic flux (aMG=0, 20, 40 or 100 mM)
pdf(file="Figure_3.pdf", width = 4, height = 4)

  cols <- rev(colorRampPalette(c("#2E75B6", "#D5E3F0"))(4))
  cols <- hcl.colors(6, "Oslo")

  # empty plot
  plot(NA, NA, type="p", log="x", xaxt='n', pch=20, las=1, xlim=c(0.1,100), ylim=c(-11,3), xlab="[acetate] (mM)", ylab="acetate flux (mmol/gDW/h)")
  abline(h=0)

  # aMG = 0 mM
  filter <- (exp_data[,"aMG"]==0)
  col <- cols[2]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"qACE"]

  points(x, y, bg=col, col="black", pch=21, lwd=1)
  # uncomment the following lines to run the fits, here we just plot the fitted curves
  # note: fitting is not stable, so it might have to be run several times to obtain a correct fit
  #res_fit <- fit_sigmoid(x, y, par_ini_sigmoid, lower_sigmoid, upper_sigmoid)
  #lines(seq(0.1,130,length.out=1000), sim_sigmoid(res_fit$par, seq(0.1,130,length.out=1000)), col=col, lwd=2)
  lines(seq(0.1,130,length.out=1000), sim_sigmoid(c(-0.2460494, -1.5989450, -0.9279984, 6.9209620), seq(0.1,130,length.out=1000)), col=col, lwd=2)
  
  # aMG = 20 mM
  filter <- (exp_data[,"aMG"]==20)
  col <- cols[3]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"qACE"]
  
  points(x, y, bg=col, col="black", pch=21, lwd=1)
  #res_fit <- fit_sigmoid(x, y, par_ini_sigmoid, lower_sigmoid, upper_sigmoid)
  #lines(seq(0.1,130,length.out=1000), sim_sigmoid(res_fit$par, seq(0.1,130,length.out=1000)), col=col, lwd=2)
  lines(seq(0.1,130,length.out=1000), sim_sigmoid(c(-2.454911e-01, -2.840827e+00, 5.006536e+01, 1.000000e+06), seq(0.1,130,length.out=1000)), col=col, lwd=2)

  # aMG = 40 mM
  filter <- (exp_data[,"aMG"]==40)
  col <- cols[4]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"qACE"]
  
  points(x, y, bg=col, col="black", pch=21, lwd=1)
  #res_fit <- fit_sigmoid(x, y, par_ini_sigmoid, lower_sigmoid, upper_sigmoid)
  #lines(seq(0.1,130,length.out=1000), sim_sigmoid(res_fit$par, seq(0.1,130,length.out=1000)), col=col, lwd=2)
  lines(seq(0.1,130,length.out=1000), sim_sigmoid(c(-0.1506843, -7.0591207, 3.5268378, 18.6152155), seq(0.1,130,length.out=1000)), col=col, lwd=2)
  
  # aMG = 100 mM
  filter <- (exp_data[,"aMG"]==100)
  col <- cols[5]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"qACE"]

  points(x, y, bg=col, col="black", pch=21, lwd=1)
  #res_fit <- fit_sigmoid(x, y, par_ini_sigmoid, lower_sigmoid, upper_sigmoid)
  #lines(seq(0.1,130,length.out=1000), sim_sigmoid(res_fit$par, seq(0.1,130,length.out=1000)), col=col, lwd=2)
  lines(seq(0.1,130,length.out=1000), sim_sigmoid(c(-0.2395667, -9.8451042, 36.7387431, 64514.4364655), seq(0.1,130,length.out=1000)), col=col, lwd=2)
  
  # add axes
  xlab_main <- c(0.1, 1, 10, 100)
  xlab_sec <- c(seq(0.2, 0.9, by=0.1), seq(2, 9, by=1), seq(20, 90, by=10))
  axis(side = 1, at = xlab_main, labels = TRUE)
  axis(side = 1, at = xlab_sec, labels = FALSE, tcl=-0.3)

dev.off()


###############################
# IMPACT OF ACETATE ON GROWTH #
###############################
# Output:                     #
#   - Figure_4.pdf: Fig. 4    #
###############################

# load experimental data

setwd(data_dir)
exp_data <- read.table("data_E_coli_K12_MG1655_Figure_3-4.txt", header=TRUE, sep="\t")

# go to results directory

setwd(results_dir)

# plot growth rate as function of acetate concentration, at high (aMG=0mM) and low (aMG=100mM) glycolytic flux

pdf(file="Figure_4.pdf", width = 4, height = 4)

  cols <- rev(colorRampPalette(c("#2E75B6", "#D5E3F0"))(4))
  cols <- hcl.colors(6, "Oslo")

  # empty plot
  plot(NA, NA, type="p", log="x", col=col, xaxt='n', pch=20, las=1, xlim=c(0.1,100), ylim=c(0, 0.7), xlab="[acetate] (mM)", ylab="growth rate (h-1)")
  
  # aMG = 0 mM
  filter <- (exp_data[,"aMG"]==0)
  col <- cols[2]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"mu"]
  
  points(x, y, bg=cols[2], col="black", pch=21, lwd=1)
  
  df <- data.frame(x=log10(x), y=y)
  df <- df[order(df$x),]
  
  model <- lm(y~poly(x,3), data = df)
  summary(model)
  range_ace <- seq(min(df$x), max(df$x), len=100)
  predy <- predict(model, data.frame(x=range_ace))
  lines(10**range_ace, predy, type="l", col=col, lwd=2)
  
  # aMG = 100 mM
  filter <- (exp_data[,"aMG"]==100)
  col <- cols[5]
  x <- exp_data[filter,"ACE"]
  y <- exp_data[filter,"mu"]

  points(x, y, bg=col, col="black", pch=21, lwd=1)
  df=data.frame(x=log10(x), y=y)
  df <- df[order(df$x),]
  
  model = lm(y~poly(x,3), data = df)
  summary(model)
  #plot(x=df$x, y=df$y, pch=20, col="grey", xlim=c(-1,2))
  #lines(10**df$x, predict(lm(y~poly(x,3), data=df)), type="l", col=cols[3], lwd=2)
  i = seq(min(df$x), max(df$x), len=100)       #  x-values for line
  predy = predict(model, data.frame(x=i))
  lines(10**i, predy, type="l", col=col, lwd=2)
  
  # add axes
  xlab_main <- c(0.1, 1, 10, 100)
  xlab_sec <- c(seq(0.2, 0.9, by=0.1), seq(2, 9, by=1), seq(20, 90, by=10))
  axis(side = 1, at = xlab_main, labels = TRUE)
  axis(side = 1, at = xlab_sec, labels = FALSE, tcl=-0.3)

dev.off()


########################################################################################################
# IMPACT OF ACETATE ON GROWTH ON DIFFERENT CARBON SOURCES (GLC, GLY, GAL) AND STRAINS (wt, Δpta, Δacs) #
########################################################################################################
# Outputs:                                                                                             #
#   - Statistical_results_glc_gly_gal.tsv: statistical results including mean ± sd of the growth rate  #
#                                          for each condition, relative change of growth rate, and     #
#                                          associated p-values                                         #
#   - Figure_6.pdf:   Fig. 6                                                                           #
#   - Figure_7A.pdf:  Fig. 7, panel A                                                                  #
#   - Figure_7B.pdf:  Fig. 7, panel B                                                                  #
########################################################################################################

# load experimental data

setwd(data_dir)
exp_data <- read.table("data_E_coli_K12_BW25113_Figure_6-7.txt", header=TRUE, sep="\t")

# go to results directory

setwd(results_dir)

# perform statistical analysis

statistical_results <- array(NA, dim=c(20, 11), dimnames=list(NULL, c("carbon source", "strain", "aMG", "ace", "mu_mean", "mu_sd", "n_replicates", "p-val (vs wt, same ace)", "p-val (vs 0 mM ace, same strain)", "relative growth change", "rgc_sd")))
k <- 0
# for each carbon source
for (c_source in unique(exp_data$carbon_source)){
  # for each strain
  for (strain in unique(exp_data$genotype)){
    # for each aMG concentration
    for (aMG in rev(unique(exp_data[((exp_data$genotype == strain) & (exp_data$carbon_source == c_source)), "aMG"]))){
      # for each acetate concentration
      for (ace in unique(exp_data$ACE)){
        # perform all statistical tests (Student’s t-tests with two-tailed distributions)
        # to compare the different strains/conditions, as detailed below
        k <- k+1
        mu_strain <- exp_data[((exp_data$genotype == strain) & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG)), "mu"]
        mu_mean <- mean(mu_strain)
        mu_sd <- sd(mu_strain)
        cat("**********************************************\n")
        cat("Statistical analysis for the following condition:\n")
        cat("\nstrain: ", strain, ", carbon source: ", c_source , ", [ace]=", ace, "mM, [aMG]=", aMG, "mM\n", sep="")
        cat("**********************************************\n")
        t_test_factor_strain <- list(p.value=NA)
        if (strain != "wt"){
          n_rep_com <- Reduce(intersect, list(exp_data[((exp_data$genotype == "wt") & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG)), "experiment_number"],
                                              exp_data[((exp_data$genotype == strain) & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG)), "experiment_number"]))
          mu_wt <- exp_data[((exp_data$genotype == "wt") & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG) & (exp_data$experiment_number %in% n_rep_com)), "mu"]
          mu_strain <- exp_data[((exp_data$genotype == strain) & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG) & (exp_data$experiment_number %in% n_rep_com)), "mu"]
          t_test_factor_strain <- t.test(mu_wt, mu_strain)
          cat("\n\nStatistical test:", paste("\n  wt vs ", strain, ", carbon source: ", c_source , ", [ace]=", ace, "mM\n", sep=""), sep="\n")
          cat("\nwt: ", paste(mu_wt, collapse=", "), "\n", sep="")
          cat(strain, ": ", paste(mu_strain, collapse=", "), "\n", sep="")
          print(t_test_factor_strain)
        }
        t_test_factor_ace <- list(p.value=NA)
        relative_growth_change <- NA
        relative_growth_change_sd <- NA
        if (ace == 10){
          mu_no_ace <- exp_data[((exp_data$genotype == strain) & (exp_data$ACE == 0) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG)), "mu"]
          mu_ace <- exp_data[((exp_data$genotype == strain) & (exp_data$ACE == ace) & (exp_data$carbon_source == c_source) & (exp_data$aMG == aMG)), "mu"]
          relative_growth_change <- mean(mu_ace/mu_no_ace)-1
          relative_growth_change_sd <- sd(mu_ace/mu_no_ace)
          # here we compared directly the relative change of growth rate induced by 10 mM acetate
          t_test_factor_ace <- t.test(mu_ace/mu_no_ace, mu=1)
          # to compare directly the growth rates in presence and absence of acetate (with paired tests, because for
          # each strain/carbon source cultivations were carried out in parallel with and without acetate), just
          # uncomment the following line (results are similar)
          # t_test_factor_ace <- t.test(mu_ace, mu_no_ace, paired=TRUE)
          cat("\n\nStatistical test:", paste("\n  strain: ", strain, ", carbon source: ", c_source , ", [ace]=0 mM vs [ace]=10 mM\n", sep=""), sep="\n")
          cat("relative change of growth rates: ", paste(mu_ace/mu_no_ace-1, collapse=", "), "\n", sep="")
          print(t_test_factor_ace)
        }
        statistical_results[k, ] <- c(c_source, strain, aMG, ace, mu_mean, mu_sd, length(mu_strain), t_test_factor_strain$p.value, t_test_factor_ace$p.value, relative_growth_change, relative_growth_change_sd)
      }
    }
  }
}

# save all results

write.table(statistical_results,"Statistical_results_glc_gly_gal.tsv", quote=FALSE, sep="\t", row.names=FALSE) 

# plot relative change in growth rate for each condition on glucose

pdf(file = "Figure_6.pdf", width = 4, height = 3.5)

  data_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "glucose") & (statistical_results[,"ace"] == "10"), "relative growth change"])
  names(data_to_plot) <- statistical_results[(statistical_results[,"carbon source"] == "glucose") & (statistical_results[,"ace"] == "10"), "strain"]
  data_sd_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "glucose") & (statistical_results[,"ace"] == "10"), "rgc_sd"])
  ze_barplot <- barplot(data_to_plot*100, ylim=c(-10,15), beside=TRUE, legend.text=FALSE, las=1, ylab="relative growth change (%)", xaxt='n')
  abline(h=0)
  error.bar(ze_barplot, data_to_plot*100, data_sd_to_plot*100)
  axis(side=1,line=0,at=ze_barplot,labels=c("wt", "wt","pta","acs"),mgp=c(0,0.5,0),tck=-0.02)
  box()

dev.off()

# plot relative change in growth rate for each condition on glycerol

pdf(file = "Figure_7A.pdf", width = 4, height = 3.5)

  data_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "glycerol") & (statistical_results[,"ace"] == "10"), "relative growth change"])
  names(data_to_plot) <- statistical_results[(statistical_results[,"carbon source"] == "glycerol") & (statistical_results[,"ace"] == "10"), "strain"]
  data_sd_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "glycerol") & (statistical_results[,"ace"] == "10"), "rgc_sd"])
  ze_barplot <- barplot(data_to_plot*100, ylim=c(-10,10), beside=TRUE, legend.text=FALSE, las=1, ylab="relative growth change (%)", xaxt='n')
  abline(h=0)
  error.bar(ze_barplot, data_to_plot*100, data_sd_to_plot*100)
  axis(side=1,line=0,at=ze_barplot,labels=c("wt", "pta","acs"),mgp=c(0,0.5,0),tck=-0.02)
  box()

dev.off()

# plot relative change in growth rate for each condition on galactose

pdf(file = "Figure_7B.pdf", width = 4, height = 3.5)

  data_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "galactose") & (statistical_results[,"ace"] == "10"), "relative growth change"])
  names(data_to_plot) <- statistical_results[(statistical_results[,"carbon source"] == "galactose") & (statistical_results[,"ace"] == "10"), "strain"]
  data_sd_to_plot <- as.numeric(statistical_results[(statistical_results[,"carbon source"] == "galactose") & (statistical_results[,"ace"] == "10"), "rgc_sd"])
  ze_barplot <- barplot(data_to_plot*100, ylim=c(-20,60), beside=TRUE, legend.text=FALSE, las=1, ylab="relative growth change (%)", xaxt='n')
  abline(h=0)
  error.bar(ze_barplot, data_to_plot*100, data_sd_to_plot*100)
  axis(side=1,line=0,at=ze_barplot,labels=c("wt", "pta","acs"),mgp=c(0,0.5,0),tck=-0.02)
  box()

dev.off()

