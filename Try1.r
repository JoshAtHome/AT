rm(list = ls())
# install.packages("did")
library('did')

set.seed(1814)
library(dplyr)
library(ggplot2)

T = 10
N = 20
shock <- T/2+2
################################################################################
# param hom decides if homogenous (0) or heterogenous (1) treatment effect 
# param par decides if a violation in parallel trends occurs

sim_data = function(hom, par){
  dat = expand.grid(t = 1:T,i = 1:N) 
  normal_random <- rnorm(T, mean = 1, sd = 0.5)
  trend <- cumsum(normal_random)
  effect <- rnorm(N*T, mean = 3, sd = 1)
  if (hom==0){
    effect <- 3
  }
  dat <- mutate(dat,
                group = ifelse(i > N/2,"treat","control"),
                G = shock*(group == "treat"),
                treat = 1L*(group == "treat"),
                control = 1L*(group == "control"), 
                pre = 1L*(t < T/2+2),
                post = 1L*(t >= T/2+2),
                pre_tre = pre*treat, pre_con = pre*control,
                post_tre = post*treat, post_con = post*control,
                tre = trend[t],
                X1 = rnorm(N*T , mean = 3, sd = 3),
                effect = effect,
                y = tre +              # the basic trend
                    treat * 5 * X1  +  # the basic differences between treat / control
                    control * X1  +    # the basic differences between treat / control
                    par * 2 * post_tre +          # the violation of parallel trends
                    post_tre * effect) # the post treatment actual effect
  dat
}

show_graph = function(dat, label = "", show.means = TRUE) {
  gdat = dat %>% group_by(group, t, post, treat) %>% summarize(y = mean(y))
  gg = ggplot(gdat, aes(y = y, x = t, color = group)) + geom_line() +  geom_vline(xintercept = T/2 + 1) +
    theme_bw()
  gg
}
d <- sim_data(0)
plot = show_graph(d)
plot

################################################################################
################################################################################
eff = function(hom){
  dat <- sim_data(hom) # simulates the data accordingly
  data <- dat %>% mutate(group = as.factor(ifelse(group == 'treat', 1, 0)))
  att <- att_gt(yname = "y",
                tname = "t",
                idname = "i",
                gname = "G",
                xformla = ~1,
                data = data) # creates the att-values
  
  p <- as.numeric(att$Wpval)  # gets the p-value from att object
  att_est <- att$att[8] # this works for T=10. Only important that POST treatment!
  att_act <- mean(data$effect) # gets the mean effect
  return(list(att_est=att_est, att_act=att_act))
}

eff_cond = function(hom){
  dat <- sim_data(hom) # simulates the data accordingly
  data <- dat %>% mutate(group = as.factor(ifelse(group == 'treat', 1, 0)))
  att <- att_gt(yname = "y",
                tname = "t",
                idname = "i",
                gname = "G",
                xformla = ~X1, 
                data = data) # creates the att-values
  
  p <- as.numeric(att$Wpval)  # gets the p-value from att object
  att_est <- att$att[8] # this works for T=10. Only important that POST treatment!
  att_act <- mean(data$effect) # gets the mean effect
  return(list(att_est=att_est, att_act=att_act))
}
eff_cond(1)
################################################################################
################################################################################
sim_eff <- function(amount, hom) {
  # amount <- number of sims
  # hom <- passed through (homogoneous/heterogenous treatment effect)
  effect_est <- numeric(amount) # all estimated effects
  effect_act <- numeric(amount)  # all actual effects
  for (i in 1:amount) {
      res <- eff(hom)
      effect_est[i] <- res$att_est
      effect_act[i] <- res$att_act
      
  }
  return(list(effect_est=effect_est, effect_act=effect_act))
}
sim_eff_cond <- function(amount, hom) { 
  # amount <- number of sims
  # hom <- passed through (homogoneous/heterogenous treatment effect)
  effect_est <- numeric(amount) 
  effect_act <- numeric(amount) 
  for (i in 1:amount) {
    res <- eff_cond(hom)
    effect_est[i] <- res$att_est
    effect_act[i] <- res$att_act
    
  }
  return(list(effect_est=effect_est, effect_act=effect_act))
}
################################################################################
################################################################################
his = function(p_vals){
  df <- data.frame(p_vals = p_vals)
  
  ggplot(df, aes(x = p_vals)) +
    geom_histogram(binwidth = 0.25, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 3, color = "red", size = 0.5) +
    geom_vline(xintercept = mean(p_vals), color = "yellow", size = 0.5) +
    labs(title = paste(length(p_vals), " simulations of ATT"),
         x = "ATT", y = "Frequency") 
}

n <- 1000


# when 



# when hom = 1, the homogeneity assumption is violated
basic_hom <- sim_eff(n,0)
basic_het <- sim_eff(n,1)
cond_hom <- sim_eff_cond(n,0)
cond_het <- sim_eff_cond(n,1)

bias1 <- abs(mean(unlist(basic_hom$effect_est), na.rm = TRUE) - mean(unlist(basic_hom$effect_act), na.rm = TRUE))
bias2 <- abs(mean(unlist(basic_het$effect_est), na.rm = TRUE) - mean(unlist(basic_het$effect_act), na.rm = TRUE))
bias3 <- abs(mean(unlist(cond_hom$effect_est), na.rm = TRUE) - mean(unlist(cond_hom$effect_act), na.rm = TRUE))
bias4 <- abs(mean(unlist(cond_het$effect_est), na.rm = TRUE) - mean(unlist(cond_het$effect_act), na.rm = TRUE))

bias1 # 0.03377619
bias2 # 0.007066046
bias3 # 0.01644108
bias4 # 0.02884108


his(cond_hom$effect_est)
################################################################################
################################################################################

# The dataset is set up using T=10, and N=20 units. Half of them are part of the treatment group, 
# and the other half are the control group. The treatment occurs in period 7 and only once. I created
# dummy variables that effect the outcome depending on pre/post and treatment/control.
# The covariates are random yet correlated with the outcome, so that the effect can be shown. The
# outcome is dependent on, and therefore explanable by the covariates. In theory, this suggests that 
# given a violation in pre-treatment parallel trends, conditioning on these covariates should decrease the 
# violation in parallel trends. By doing so, I should be able to use the simulated data to demonstrate the plausability 
# of the relaxed parallel trends assumption, the conditional parallel trends.

################################################################################
















# After creating the data, the estimation takes place. Using the att_gt function () I am able estimate the
# average group-time treatment effect on the treated. The effect itself is for my purposes less relevant,
# however I include the results in figure X for transparency and insight. The function
# get_p_value and get_p_value_cond return a set of p-values and ATTs, for the basic model,
# and the model which includes covariates, respectively.

################################################################################
################################################################################
# returns p-value and ATT in post-treatment period 
# param h carried over for sim_data depending on type of treatment effect
get_p_value = function(h){
  dat <- sim_data(h) # simulates the data accordingly
  data <- dat %>% mutate(group = as.factor(ifelse(group == 'treat', 1, 0)))
  att <- att_gt(yname = "y",
                tname = "t",
                idname = "i",
                gname = "G",
                xformla = ~1,
                data = data) # creates the att-values
                
  p <- as.numeric(att$Wpval)  # gets the p-value from att object
  effect <- att$att[8] # this works for T=10. Only important that POST treatment!
  return(list(p = p, effect = effect))
}

# returns CONDITIONAL p-value and ATT in post-treatment period 
get_p_value_cond = function(h){ 
  dat <- sim_data(h) # simulates the data accordingly
  data <- dat %>% mutate(group = as.factor(ifelse(group == 'treat', 1, 0)))
  att <- att_gt(yname = "y",
                tname = "t",
                idname = "i",
                gname = "G",
                xformla = ~X1,
                data = data) # creates the att-values
  
  p <- as.numeric(att$Wpval)  # gets the p-value from att object
  effect <- att$att[8] # this works for T=10. Only important that POST treatment!
  return(list(p = p, effect = effect))
}
################################################################################
################################################################################

# This process is repeated 10000 times. It yields a set of p-values indicating rejections / non-rejections. 
# The distribution of these p-values shows us the amount of times our test would reject 
# a violation of the parallel trends assumption, as seen in figure 3. 
# Figure 3a displays the test results under no conditioning on covariates. 
# 60 percent of the the tests are rejected indicating that in 60 percent of the cases our simlation was
#  able to reject a violation in parallel trends. Comparing this to the conditional model reveals a higher, 70 percent
# rejection. 


# finish
# A higher rejection rate for the conditional model indicates that, theory indeed holds true, and that 
# conditioning on covariates, can lead to a more plausable (condtitional) parallel trends assumption.

################################################################################
################################################################################

# simulates the process amount times.
# returns a vector of (conditional) p-values
# returns a vector of the (conditional) effects
################################################################################
sim_p_values <- function(amount, with, h) { 
  # amount <- number of sims
  # with <- 0 normal, 1 conditional on X1
  # h <- passed through (homogoneous/heterogenous treatment effect)
  p_vector <- numeric(amount) 
  effect_vector <- numeric(amount) 
  if (with == 0) {
    for (i in 1:amount) {
      res <- get_p_value(h)
      p_vector[i] <- res$p
      effect_vector[i] <- res$effect
    }
  } else {
    for (i in 1:amount) {
      res <- get_p_value_cond(h)
      p_vector[i] <- res$p
      effect_vector[i] <- res$effect
    }
  }
  return(list(p_vector = p_vector, effect_vector = effect_vector))
}

################################################################################
# creates a histogram of all the simulated p-values
make_hist = function(p_vals){
 hist(p_vals, breaks = 200, col = "skyblue", border = "black",
      main = paste(length(p_vals), ' p-values'),
      xlab = "p-values", ylab = "Frequency")
 abline(v = 0.05, col = "red", lwd = 2)
 above <- sum(p_vals > 0.05)
 mtext(paste("Paralle trends:", 100*(above/length(p_vals)), ' %'), side = 1, line = 2, at = 0.5, col = "blue", cex = 0.8)
}
library(ggplot2)

library(grid)
make_hist2 = function(p_vals){
  df <- data.frame(p_vals = p_vals)
  
  ggplot(df, aes(x = p_vals)) +
    geom_histogram(binwidth = 0.0025, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.5) +
    labs(title = paste(length(p_vals), " p-values"),
         x = "p-values", y = "Frequency") +
 # annotate("text", x = 0.05, y = 100, label = paste("Paralle trends:", round(100 * sum(p_vals > 0.05) / length(p_vals), 2), "%"),
 #            vjust = 0, hjust = -0.2, color = "blue", size = 4)

    theme(aspect.ratio = 1) 
 # grid.text("Hello", x = unit(0.5, "npc"), y = unit(0, "npc"), just = "bottom", gp = gpar(fontsize = 12, fontface = "bold")) + 
  grid.rect(gp = gpar(col = "black", fill = NA)) + 
  
  # Arrange the plot and rectangle box
  grid.arrange(hist_plot, heights = c(4, 1))
  
}
n <- 20
sim_0 <- sim_p_values(n, 0, 1)  
make_hist2(sim_0$p_vector)


# Example usage:
# p_values <- runif(1000)
# make_hist(p_values)

################################################################################

# plot the comparison between the basic ATTs and the conditional ATTs
compare_att <- function(n, h){
  sim_0 <- sim_p_values(n, 0, h) # simulate basic
  sim_1 <- sim_p_values(n, 1, h) # simulate conditional
  eff_0 <- sim_0$effect_vector # retrieve basic ATT vector
  eff_1 <- sim_1$effect_vector # retrieve conditional ATT vector
  ci_0 <- t.test(eff_0)$conf.int # calculate basic confidence interval
  ci_1 <- t.test(eff_1)$conf.int # calculate conditional confidence interval
  
  # put into dataframe
  data <- data.frame(group = c("basic", "conditional"),
                     mean = c(mean(eff_0), mean(eff_1)),
                     lower = c(ci_0[1], ci_1[1]),
                     upper = c(ci_0[2], ci_1[2]))
  # plot nicely
  ggplot(data, aes(x = group, y = mean)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, color = "black") +
    labs(title = "Simulated ATTs ",
         x = "Group",
         y = "ATT") +
    theme_minimal()
}

show_graph = function(dat, label = "", show.means = TRUE) {
  gdat = dat %>% group_by(group, t, post, treat) %>% summarize(y = mean(y))
  gg = ggplot(gdat, aes(y = y, x = t, color = group)) + geom_line() +  geom_vline(xintercept = T/2 + 1) +
    theme_bw()
  gg
}
d <- sim_data(1)
plot = show_graph(d)
plot

################################################################################



library(gridExtra)

n <- 20
sim_0 <- sim_p_values(n, 0, 1)  
sim_1 <- sim_p_values(n, 1, 1) 
make_hist(sim_1$p_vector)

a <- compare_att(10, 1)
b <- compare_att(12, 0)
grid.arrange(a,b, nrow = 2)


c <- make_hist2(sim_0$p_vector)
d <- make_hist2(sim_1$p_vector)
c


grid.arrange(c, d, nrow=2)

df <- data.frame(
  p=sim_0$p_vector
)
head(df)





p<-ggplot(df, aes(x=p)) + 
   geom_histogram(binwidth=0.005, color="black", fill="skyblue", alpha=0.9) + 
   geom_vline(aes(xintercept=0.05), color="red", linetype="solid", size=0.5)
   ggtitle("Bin size = 3") 
p






