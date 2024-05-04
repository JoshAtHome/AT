rm(list = ls())

wd <- '/Users/josh/Documents/Laptop/SS24/AdvancedTopics/codes2'
setwd(wd)

# install.packages("HonestDiD")
# install.packages("remotes")
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
# remotes::install_github("asheshrambachan/HonestDiD")
# install.packages("fixest")

library(fixest)
library(HonestDiD)
library(did)
library(dplyr)
library(ggplot2)
set.seed(4)

################################################################################
#############################SIMULATE DATA######################################
################################################################################

sim_data = function(heterogeneity, parallel){
  periods = 7
  params <- reset.sim(time.periods = periods, n = 1000, ipw = TRUE, reg = FALSE)
  params$te <- 10
  h <- heterogeneity
  p <- parallel
  if (parallel==1){
    params$theu <- params$thet + rnorm(periods, mean = 0, sd = 2)  # violate pt
    # params$betu <- params$bett + rnorm(periods, mean = 0, sd = 1)  # violate cpt
  }
  if (heterogeneity==1){
    het_eff <- rnorm(periods, mean = 10, sd = 5)
    params$te.bet.X <- het_eff
  }
  data <- build_sim_dataset(sp_list=params, panel = TRUE)
  data2 <- data[data$G %in% c(0, 7), ]
  data2 <- mutate(data2,
                  group = ifelse(treat == 1, 'treatment', 'control'),
                  control = ifelse(treat != 1, 1, 0),
                  X1 = Y^2 / 2,
                  effect_act = ifelse(heterogeneity == 1, mean(het_eff), 10))
  return(data2)
}

sim_data(0,0)

show <- function(dat, label = "", show.means = TRUE) {
  gdat <- dat %>% group_by(group, period) %>% summarize(y = mean(Y))
  gg <- ggplot(gdat, aes(y = y, x = period, color = group)) + geom_line() + geom_vline(xintercept = 6.5)+theme_bw()
  print(gg)}


################################################################################
#############################ATT-CALCULATION####################################
################################################################################

att_standard = function(heterogeneity, parallel){
  dat <- sim_data(heterogeneity, parallel) # simulates the data accordingly
  att <- att_gt(yname = "Y",
                tname = "period",
                idname = "id",
                gname = "G",
                xformla = ~1,
                data = dat) # creates the att-values
 # p <- as.numeric(att$Wpval) 
  att_est <- att$att[6]
  att_act <- dat$effect_act[1]
  beta <- att$att
  sigma <- as.matrix(att$V_analytical)
  return(list(att_est=att_est, att_act=att_act, beta=beta,sigma=sigma))
}

att_conditional = function(heterogeneity, parallel){
  dat <- sim_data(heterogeneity, parallel) # simulates the data accordingly
  att <- att_gt(yname = "Y",
                tname = "period",
                idname = "id",
                gname = "G",
                xformla = ~X,
                data = dat) # creates the att-values
#  p <- as.numeric(att$Wpval) 
  att_est <- att$att[6]
  att_act <- dat$effect_act[1]
  beta <- att$att
  sigma <- as.matrix(att$V_analytical)
  return(list(att_est=att_est, att_act=att_act, beta=beta,sigma=sigma))
}

a<- att_standard(1,0)
a
################################################################################
#############################SIMULATION#########################################
################################################################################

sim <- function(amount, conditional, heterogeneity, parallel) {
  # amount <- number of simulations
  # conditional <- condition on covariate
  # heterogeneity <- passed through (homogoneous/heterogenous treatment effect)
  # parallel <- passed through (violation of parallel trends)
  effect_est <- numeric(amount) # all estimated effects
  effect_act <- numeric(amount)  # all actual effects
  p_values <- numeric(amount)
  if (conditional==0){
    for (i in 1:amount) {
      res <- att_standard(heterogeneity, parallel)
      effect_est[i] <- res$att_est
      effect_act[i] <- res$att_act
 #     p_values[i] <- res$p
      
    }}
  else if (conditional==1){
    for (i in 1:amount) {
      res <- att_conditional(heterogeneity, parallel)
      effect_est[i] <- res$att_est
      effect_act[i] <- res$att_act
 #     p_values[i] <- res$p
    }}
  return(list(effect_est=effect_est, effect_act=effect_act, p_values=p_values))
}

################################################################################
#############################Results############################################
################################################################################

# sim(amount, conditional, heterogeneity, parallel)
n <- 100

standard_hom <- sim(n,0,0,0)
standard_het <- sim(n,0,1,0)
cond_hom <- sim(n,1,0,0)
cond_het <- sim(n,1,1,0)
standard_hom_vio <- sim(n,0,0,1)
standard_het_vio <- sim(n,0,1,1)
cond_hom_vio <- sim(n,1,0,1)
cond_het_vio <- sim(n,1,1,1)


library(knitr)
names <- c("standard_hom", "standard_het", "cond_hom", "cond_het",
                   "standard_hom_vio", "standard_het_vio", "cond_hom_vio", "cond_het_vio")
biases <- numeric(length(names))
for (i in seq_along(names)) {
  biases[i] <- abs(mean(unlist(get(names[i])$effect_est), na.rm = TRUE) -
                     mean(unlist(get(names[i])$effect_act), na.rm = TRUE))
}
bias_table <- data.frame(
  Parallel_Trends = c('no violation','no violation','no violation','no violation',
                      'violation','violation','violation','violation'),
  DiD_Model = c('Standard','Standard','Conditional','Conditional',
                'Standard','Standard','Conditional','Conditional'),
  Treatment_effect = c('homogeneous','heterogeneous','homogeneous','heterogeneous',
              'homogeneous','heterogeneous','homogeneous','heterogeneous'),
  Bias = biases
)

bias_table %>%
  kable(caption = "Summary of biases",
        format = "latex",
        col.names = c("Parallel_Trends","DiD_Model", "Treatment_effect", "Bias"),
        align = "c")

kable(bias_table, caption = "Biases for violations")

################################################################################
# p_values_df <- data.frame(standard_hom = standard_hom$p,
#                           standard_het = standard_het$p,
#                           cond_hom = cond_hom$p,
#                           cond_het = cond_het$p)
# library(ggplot2)
# library(gridExtra)
# make_hist <- function(p_vals){
#   plots <- lapply(names(p_vals), function(col_name) {
#     ggplot(p_vals, aes(x = .data[[col_name]])) +
#       geom_histogram(binwidth = 0.0025, fill = "skyblue", color = "black") +
#       geom_vline(xintercept = 0.05, color = "red", size = 0.5) +
#       labs(title = col_name,
#            x = "p-values", y = "freq") 
#   })
#   grid.arrange(grobs = plots, ncol = 2)
# }
# make_hist(p_values_df)
################################################################################






################################################################################
############################sensitivity#########################################
################################################################################
# devtools::install_github("jonathandroth/pretrends")
# install_github("pedrohcgs/CS_RR", dependencies = TRUE)

get_betasigma <- function(conditional, heterogeneity, parallel){
  if (conditional==0){
    obj <- att_standard(heterogeneity, parallel)
  }
  else if (conditional==1){
    obj <- att_conditional(heterogeneity, parallel)
  }
  beta <- obj$beta
  sigma <- obj$sigma
  return(list(beta=beta,sigma=sigma))
}

get_rm <- function(conditional, heterogeneity, parallel){
  results <- replicate(100, get_betasigma(conditional, heterogeneity, parallel), simplify = FALSE)
  beta_sum <- Reduce(`+`, lapply(results, function(x) x$beta))
  beta_mean <- beta_sum / length(results)
  sigma_sum <- Reduce(`+`, lapply(results, function(x) x$sigma))
  sigma_mean <- sigma_sum / length(results)
  beta <- beta_mean
  sigma <- sigma_mean
  original  <- constructOriginalCS(betahat=beta, sigma=sigma, numPrePeriods=5, numPostPeriods=1)
  robust <- createSensitivityResults_relativeMagnitudes(betahat=beta,sigma=sigma,
           numPrePeriods=5, numPostPeriods=1, gridPoints=100, Mbarvec=seq(0.5, 2, by = 0.25))
  p <- createSensitivityPlot_relativeMagnitudes(robust, original)
  return(p)
#  return(list(original=original, robust=robust))
}


rm1 <- get_rm(0,0,0)
rm2 <- get_rm(0,1,0)
rm3 <- get_rm(1,0,0)
rm4 <- get_rm(1,1,0)
grid.arrange(rm1, rm2, rm3, rm4, ncol = 2)

# Starting off with no violation in parallel trends. Only under homogenous treatment effects and using the 
# standard model are the results from above are robust.  
# Here, just under an equal size violation of parallel trends is allowed compared to the maximum
# in pre-treatment periods.

set.seed(123)
rm5 <- get_rm(0,0,1)
rm6 <- get_rm(0,1,1)
rm7 <- get_rm(1,0,1)
rm8 <- get_rm(1,1,1)
grid.arrange(rm5, rm6, rm7, rm8, ncol = 2)



