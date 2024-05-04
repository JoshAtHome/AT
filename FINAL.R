rm(list = ls())

set.seed(181)
library(dplyr)
library(ggplot2)

T = 7
N = 20
shock <- T-1
################################################################################
# param hom decides if homogenous (0) or heterogenous (1) treatment effect 
# param par decides if a violation in parallel trends occurs

sim_data = function(heterogeneity, violation){
  dat = expand.grid(t = 1:T,i = 1:N) 
  normal_random <- rnorm(T, mean = 0.5, sd = 0.5)
  trend <- cumsum(normal_random)
  
  dat$effect <- rnorm(N*T, mean = 3, sd = 3)
  if (heterogeneity!=1){
    dat$effect <- rnorm(N*T, mean = heterogeneity, sd = 3)
  }
  if (heterogeneity==0){
    dat$effect <- 3
  }
  
  dat$violation <- rnorm(N*T, mean = violation, sd = 3)
  if (violation==1){
    dat$violation <- rnorm(N*T, mean = 0, sd = 3)
  }
  if (violation==0){
    dat$violation <- 0
  }
  
  dat <- mutate(dat,
                group = ifelse(i > N/2,"treat","control"),
                G = shock*(group == "treat"),
                treat = 1L*(group == "treat"),
                control = 1L*(group == "control"), 
                current = 1L*(G == t), 
                pre = 1L*(t < T-1),
                post = 1L*(t >= T-1),
                pre_tre = pre*treat, pre_con = pre*control,
                post_tre = post*treat, post_con = post*control,
                tre = trend[t],
                X1 = rep(rnorm(N, mean = 3, sd = 1), each = T),
                y = tre +              # the basic trend
                  treat * 2 * X1 -   # the basic differences between treat / control
                  control * X1  +    # the basic differences between treat / control
                  violation * post_tre +          # the violation of parallel trends
                  current * effect) # the treatment effect
  dat
}

show_graph = function(dat, label = "", show.means = TRUE) {
  gdat = dat %>% group_by(group, t, post, treat) %>% summarize(y = mean(y))
  gg = ggplot(gdat, aes(y = y, x = t, color = group)) + geom_line() +  geom_vline(xintercept = T-1) +
    theme_bw()
  gg
}
d <- sim_data(0,1)
plot = show_graph(d)
plot


################################################################################
################################################################################

att_estimate <- function(heterogeneity, violation, x_formula){
  dat <- sim_data(heterogeneity, violation) # simulates the data accordingly
  att <- att_gt(yname = "y",
                tname = "t",
                idname = "i",
                gname = "G",
                xformla = x_formula,
                data = dat) # creates the att-values
  att_est <- att$att[5]
  att_act <- mean(dat$effect[dat$treat == 1])
  beta <- att$att
  sigma <- as.matrix(att$V_analytical)
  return(list(att_est=att_est, att_act=att_act, beta=beta,sigma=sigma))
}
att_estimate(1,1,~X1)

sim <- function(amount, conditional, heterogeneity, violation) {
  effect_est <- numeric(amount) # all estimated effects
  effect_act <- numeric(amount)  # all actual effects
  p_values <- numeric(amount)
  if (conditional==0){
    for (i in 1:amount) {
      res <- att_estimate(heterogeneity, violation, ~1)
      effect_est[i] <- res$att_est
      effect_act[i] <- res$att_act
    }}
  else if (conditional==1){
    for (i in 1:amount) {
      res <- att_estimate(heterogeneity, violation, ~X1)
      effect_est[i] <- res$att_est
      effect_act[i] <- res$att_act
    }}
  return(list(effect_est=effect_est, effect_act=effect_act, p_values=p_values))
}
################################################################################
################################################################################
################################################################################
################################################################################
library(knitr)

n<-200

# sim(n, conditional, heterogeneity, violation) --- REAL ATT = 3
# @param conditional 0,1
# @param heterogeneity 0: effect = 3, x: effect = random, mean x
# @param violation 0: no violation, 1: mean violation of 0, x: mean violation of x

#################### no violation ##############################################
namev0 <- c('sim000','sim010','sim050','sim100','sim110','sim150')
sim000 <- sim(n,0,0,0) # standard / homogeneity
sim010 <- sim(n,0,1,0) # standard / mean heterogeneity = true effect
sim050 <- sim(n,0,5,0) # standard / mean heterogeneity = 5
sim100 <- sim(n,1,0,0) # conditional / homogeneity
sim110 <- sim(n,1,1,0) # conditional / mean heterogeneity = true effect
sim150 <- sim(n,1,5,0) # conditional / mean heterogeneity = 5

#################### violation mean 0 ##########################################
namev1 <- c('sim001','sim011','sim051','sim101','sim111','sim151')
sim001 <- sim(n,0,0,1) # standard / homogeneity
sim011 <- sim(n,0,1,1) # standard / mean heterogeneity = true effect
sim051 <- sim(n,0,5,1) # standard / mean heterogeneity = 5
sim101 <- sim(n,1,0,1) # conditional / homogeneity
sim111 <- sim(n,1,1,1) # conditional / mean heterogeneity = true effect
sim151 <- sim(n,1,5,1) # conditional / mean heterogeneity = 5

#################### violation mean 3 ##########################################
namev3 <- c('sim003','sim013','sim053','sim103','sim113','sim153')
sim003 <- sim(n,0,0,3) # standard / homogeneity
sim013 <- sim(n,0,1,3) # standard / mean heterogeneity = true effect
sim053 <- sim(n,0,5,3) # standard / mean heterogeneity = 5
sim103 <- sim(n,1,0,3) # conditional / homogeneity
sim113 <- sim(n,1,1,3) # conditional / mean heterogeneity = true effect
sim153 <- sim(n,1,5,3) # conditional / mean heterogeneity = 5


create_bias_table <- function(names, vio) {
  biases <- numeric(length(names))
  for (i in seq_along(names)) {
    biases[i] <- abs(mean(unlist(get(names[i])$effect_est), na.rm = TRUE) -
                       mean(unlist(get(names[i])$effect_act), na.rm = TRUE))}
  bias_table <- data.frame(
    DiD_Model = c('Standard', 'Standard', 'Standard', 'Conditional', 'Conditional', 'Conditional'),
    Treatment_effect = c('homogeneous','heterogeneous (mean=3)','heterogeneous (mean=5)'),
    Bias = biases)
  ans <- ifelse(vio == 0, 'No violation', ifelse(vio == 1, 'Violation (mean=0)', 
                                                 ifelse(vio == 3, 'Violation (mean=3)', 0)))
  #kable(bias_table, caption = paste("Biases for", ans, "of the Parallel Trends Assumption"))
  bias_table
}

tab0 <- create_bias_table(namev0, 0)
tab1 <- create_bias_table(namev1, 1)
tab3 <- create_bias_table(namev3, 3)

tab0
tab1
tab3


tab0 %>%
  kable(caption = "Biases for No Violation of the Parallel Trends Assumption",
        format = "latex",
        col.names = c("DiD_Model", "Treatment_effect", "Bias"),
        align = "c")

################################################################################
################################################################################







