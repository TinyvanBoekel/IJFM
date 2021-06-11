# This script file contains the R code for the paper "To pool or not to pool: that is the question in microbial kinetics" by M.A.J.S. van Boekel,
# published in the International Journal of Food Microbiology, 2021. https://doi.org/10.1016/j.ijfoodmicro.2021.109283

# packages used:
library(papaja)
library(ggplot2)
library(dplyr)
library(brms)
library(tidyverse)
library(rstan)
library(broom)
library(patchwork)
library(here)
library(GGally)
library(papaja)
library(tidybayes)
library(ggridges)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

knitr::opts_chunk$set(echo = FALSE, 
                      fig.align = "center",
                      fig.width = 6, 
                      fig.height = 4, 
                      dev = "png",
                      warning = FALSE,
                      message = FALSE,
                      comment = NA)

theme_set(theme_bw())

# Loading the data for the first case study followed by a plot (Figure 2 in main article):
mattick65 <- read.csv(here("data","Mattick_65_080.csv"), header=TRUE, sep=";")

mattick65$trial <- as.factor(mattick65$trial)

(dataplot1 <- mattick65 %>% ggplot(aes(x=time, y=logN, colour=trial))+
    geom_point(aes(y=logN))+
    geom_line()+
    labs(x = "time (h)", y = expression(paste("log"[10], "N")))+
    theme(legend.position = "none")
)

# Loading the data for the second case study:
salmo <- read.csv(file=here("data","Combase_data_set.csv"), header=TRUE, sep=",")
salmo$trial <- as.factor(salmo$trial)
salmo$exp_no <- as.factor(salmo$exp_no)
#preparing the data at each temp:
salmo80 <- subset(salmo, temp ==80, select=trial:logN)
salmo78 <- subset(salmo, temp ==78, select=trial:logN)
salmo76 <- subset(salmo, temp ==76, select=trial:logN)
salmo74 <- subset(salmo, temp ==74, select=trial:logN)
salmo72 <- subset(salmo, temp ==72, select=trial:logN)
salmo70 <- subset(salmo, temp ==70, select=trial:logN)
salmo65 <- subset(salmo, temp ==65, select=trial:logN)
salmo60 <- subset(salmo, temp ==60, select=trial:logN)
salmo55 <- subset(salmo, temp ==55, select=trial:logN)

# prior predictive check (Supplemental material Figure S1)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha", lb=0),
           prior(normal(0.01,0.01), nlpar="beta",lb=0),
           prior(cauchy(0,10), class="sigma")
)

mattick_ppc <- brm(formula=nlform, data=mattick65, family = gaussian(), sample_prior="only", prior = nlprior, chains= 4, warmup=1000, iter = 2000, control = list(adapt_delta = 0.9), file=here("fits", "mattick_ppc"))

fe_only <- tibble(time=seq(min(mattick65$time),max(mattick65$time), length.out = 100)) %>% add_fitted_draws(mattick_ppc, re_formula = NA, n=200)

ggplot(fe_only, aes(x=time, y=.value, group=.draw))+
  geom_line(colour="blue")+
  ylim(0,8)+
  labs(x="time (h)", y=expression(paste(log[10],"N")))


# Bayesian Regression with brms on pooled data, case study 1:
nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(5,5), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick_pooled2 <- brm(formula=nlform, data=mattick65, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick_pooled2"))

# collecting the posterior samples in a dataframe:
mattick_pooled2_post <- posterior_samples(mattick_pooled2)
print(mattick_pooled2, digits=4)
pairs(mattick_pooled2)

# trace plots (Supplemental material Figure S2)
bayesplot::mcmc_trace(mattick_pooled2_post)+labs(x="number of post-warmup samples")

# pairs and correlation plot (Figure 3 in main article):
mattick_cor1 <- dplyr::select(mattick_pooled2_post,b_logN0_Intercept:sigma)
mattick_cor1 <- setNames(mattick_cor1, c(expression(paste("log"[10],"N")), 
                                         expression(alpha),expression(beta), expression(sigma[e])))


(mattick_corplot1 <-mattick_cor1  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")) +
    theme(axis.text.x = element_text(angle=70, hjust=1))
)

# Preparing for Table 1, case study 1 using papaja:
mattick_summary1 <- summary(mattick_pooled2)
mattick_summary2 <- rbind(data.frame(mattick_summary1$fixed)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(mattick_summary1$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

rownames(mattick_summary2) <- c("$\\log N_0$", "$\\alpha$", "$\\beta (h)$", "$\\sigma_e$")
colnames(mattick_summary2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  mattick_summary2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:pooled)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(3,3,3,3)
  ),
  escape = FALSE
)

# Case study 1: preparing the data for averaging, followed by Bayesian regression using brms
t0 <- mattick65 %>% dplyr::filter(time==0)
t1 <- mattick65 %>% dplyr::filter(time==0.03)
t2 <- mattick65 %>% dplyr::filter(time==0.07)
t3 <- mattick65 %>% dplyr::filter(time==0.1)
t4 <- mattick65 %>% dplyr::filter(time==0.13)
t5 <- mattick65 %>% dplyr::filter(time==0.17)
t6 <- mattick65 %>% dplyr::filter(time==0.2)
t7 <- mattick65 %>% dplyr::filter(time==0.23)
t8 <- mattick65 %>% dplyr::filter(time==0.27)
t9 <- mattick65 %>% dplyr::filter(time==0.3)
time_avg <- c(0,0.03,0.07,0.1,0.13,0.17,0.2,0.23,0.27,0.3)
mean_avg <- c(mean(t0$logN), mean(t1$logN),mean(t2$logN),mean(t3$logN),mean(t4$logN),mean(t5$logN),mean(t6$logN),mean(t7$logN), mean(t8$logN), mean(t9$logN))
sd_avg <- c(sd(t0$logN), sd(t1$logN),sd(t2$logN),sd(t3$logN),sd(t4$logN),sd(t5$logN),sd(t6$logN),sd(t7$logN), sd(t8$logN), sd(t9$logN))
avg_df <- data.frame(time_avg, mean_avg, sd_avg)

nlform<-bf(mean_avg ~ logN0-(1/2.303)*(time_avg/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(5,5), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick_avg2 <- brm(formula=nlform, data=avg_df, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick_avg2"))
# collecting the posterior samples in a dataframe:
mattick_avg2_post <- posterior_samples(mattick_avg2)
# numerical regression results (summary of the posterior)
print(mattick_avg2)
# a pairs plot of the marginal posteriors
pairs(mattick_avg2)

# Plotting the fits to the data:

# first calculating the time range to plot:

time.seq <- data.frame(time_avg = seq(from = 0, to = 0.3, by = 0.01))

#fitted is about the mean mu:
muSummary <-
  fitted(mattick_avg2, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

#predict is about future individual values:
pred.mattick_avg <-
  predict(mattick_avg2,
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

#plot of fitted and predicted averaged values
avg_df_plot <- avg_df %>%
  ggplot(aes(x = time_avg, y = mean_avg)) +
  geom_ribbon(data = pred.mattick_avg, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_errorbar(aes(ymin=mean_avg-sd_avg, ymax=mean_avg+sd_avg))+
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.5) +
  geom_line(data = muSummary, aes(y = Estimate)) +
  geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3) +
  labs(x="time (h)", y=expression(paste("log N"[10],"N")), subtitle = "B")

# Pairs and correlation plots for the pooled and averaged data (Supplement Figure S3):
mattick_cor1 <- dplyr::select(mattick_avg2_post,b_logN0_Intercept:sigma)
mattick_cor1 <- setNames(mattick_cor1, c(expression(paste("log"[10],"N")), 
                                         expression(alpha),expression(beta), expression(sigma[e])))


(mattick_corplot1 <-mattick_cor1  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")) +
    theme(axis.text.x = element_text(angle=70, hjust=1))
)

#Preparing to plot the pooled and averaged data (Figure 4 in main manuscript)

# first calculating the time range to plot:

time.seq <- data.frame(time = seq(from = 0, to = 0.3, by = 0.005))

#fitted is about mu:
muSummary <-
  fitted(mattick_pooled2, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

#predict is about future individual values:
pred.mattick1 <-
  predict(mattick_pooled2,
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# regression line + data points
mattick_plot1 <- mattick65 %>%
  ggplot(aes(x = time, y = logN)) +
  geom_line(data = muSummary, aes(y = Estimate)) +
  geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3) +
  labs(x="time (h)")

# regression line + data points + 95% credible interval
mattick_plot2 <- mattick_plot1 + 
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.5)

# regression line + data points + 95% prediction interval
mattick_plot3 <- mattick_plot1+
  geom_ribbon(data = pred.mattick1, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue", alpha=0.4)

mattick_plot4 <- mattick_plot3 +
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.5) +
  labs(subtitle="A")

mattick_plot4 + avg_df_plot

# Code for Table S1 in Supplement:
mattick_summary1 <- summary(mattick_avg2)
mattick_summary2 <- rbind(data.frame(mattick_summary1$fixed)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(mattick_summary1$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

rownames(mattick_summary2) <- c("$\\log_{10} N_0$", "$\\alpha$", "$\\beta (h)$", "$\\sigma_e$")
colnames(mattick_summary2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  mattick_summary2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:pooled-avg)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,2,2,2)
  ),
  escape = FALSE
)

# Regression codes for the 18 no-pooling cases in case-study 1.

# trial 1:
trial1 <- filter(mattick65, trial==1)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.01,0.01), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial1 <- brm(formula=nlform, data=trial1, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial1"))

post_trial1 <- posterior_samples(mattick2_trial1)

print(mattick2_trial1, digits=4)

# composing a dataframe containing fitted values for the individual regressions (no pooling):

fit_trial1 <- cbind(time.seq, fitted(mattick2_trial1, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial1$trial <- as.factor(rep(1, length(nrow(time.seq))))

# calculating the prediction intervals and regression line (Figure S4 in Supplement)
predline1 <-
  predict(mattick2_trial1, 
          newdata = time.seq) %>%
 as_tibble() %>%
  bind_cols(time.seq)


regrline1 <- ggplot(data = trial1, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline1, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial1, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 1")


# trial 2
trial2 <- filter(mattick65, trial==2)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.03,0.03), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick2_trial2 <- brm(formula=nlform, data=trial2, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial2"))

post_trial2 <- posterior_samples(mattick2_trial2)

print(mattick2_trial2)

fit_trial2 <- cbind(time.seq, fitted(mattick2_trial2, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial2$trial <- as.factor(rep(2, length(nrow(time.seq))))

# calculating the prediction intervals and regression line for each no-pooled case (piling up in Figure S4 in Supplement)
predline2 <-
  predict(mattick2_trial2, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline2 <- ggplot(data = trial2, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline2, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial2, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 2")


# Trial 3:
trial3 <- filter(mattick65, trial==3)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick2_trial3 <- brm(formula=nlform, data=trial3, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial3"))

post_trial3 <- posterior_samples(mattick2_trial3)

print(mattick2_trial3)

fit_trial3 <- cbind(time.seq, fitted(mattick2_trial3, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial3$trial <- as.factor(rep(3, length(nrow(time.seq))))

predline3 <-
  predict(mattick2_trial3, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline3 <- ggplot(data = trial3, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial3, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 3")


# trial 4:
trial4 <- filter(mattick65, trial==4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(5,5), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial4 <- brm(formula=nlform, data=trial4, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial4"))

post_trial4 <- posterior_samples(mattick2_trial4)

print(mattick2_trial4, digits=4)

fit_trial4 <- cbind(time.seq, fitted(mattick2_trial4, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial4$trial <- as.factor(rep(4, length(nrow(time.seq))))

predline4 <-
  predict(mattick2_trial4, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline4 <- ggplot(data = trial4, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline4, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial4, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 4")

# trial 5:

trial5 <- filter(mattick65, trial==5)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial5 <- brm(formula=nlform, data=trial5, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial5"))

post_trial5 <- posterior_samples(mattick2_trial5)

print(mattick2_trial5)

fit_trial5 <- cbind(time.seq, fitted(mattick2_trial5, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial5$trial <- as.factor(rep(5, length(nrow(time.seq))))

predline5 <-
  predict(mattick2_trial5, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline5 <- ggplot(data = trial5, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline5, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial5, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 5")


# trial 6:
trial6 <- filter(mattick65, trial==6)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial6 <- brm(formula=nlform, data=trial6, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial6"))

post_trial6 <- posterior_samples(mattick2_trial6)

print(mattick2_trial6)

fit_trial6 <- cbind(time.seq, fitted(mattick2_trial6, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial6$trial <- as.factor(rep(6, length(nrow(time.seq))))

predline6 <-
  predict(mattick2_trial6, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline6 <- ggplot(data = trial6, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial6, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 6")

#trial 7:
trial7 <- filter(mattick65, trial==7)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.055,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial7 <- brm(formula=nlform, data=trial7, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial7"))

post_trial7 <- posterior_samples(mattick2_trial7)

print(mattick2_trial7)

fit_trial7 <- cbind(time.seq, fitted(mattick2_trial7, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial7$trial <- as.factor(rep(7, length(nrow(time.seq))))

predline7 <-
  predict(mattick2_trial7, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline7 <- ggplot(data = trial7, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline7, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial7, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 7")


#trial 8:

trial8 <- filter(mattick65, trial==8)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.055,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial8 <- brm(formula=nlform, data=trial8, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial8"))

post_trial8 <- posterior_samples(mattick2_trial8)

print(mattick2_trial8)

fit_trial8 <- cbind(time.seq, fitted(mattick2_trial8, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial8$trial <- as.factor(rep(8, length(nrow(time.seq))))

predline8 <-
  predict(mattick2_trial8, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline8 <- ggplot(data = trial8, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline8, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial8, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 8")


# trial 9:

trial9 <- filter(mattick65, trial==9)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial9 <- brm(formula=nlform, data=trial9, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial9"))

post_trial9 <- posterior_samples(mattick2_trial9)

print(mattick2_trial9)

fit_trial9 <- cbind(time.seq, fitted(mattick2_trial9, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial9$trial <- as.factor(rep(9, length(nrow(time.seq))))

predline9 <-
  predict(mattick2_trial9, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline9 <- ggplot(data = trial9, 
                    aes(x = time, y = logN)) +
  geom_ribbon(data = predline9, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial9, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 9")


# trial 10:

trial10 <- filter(mattick65, trial==10)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial10 <- brm(formula=nlform, data=trial10, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial10"))

post_trial10 <- posterior_samples(mattick2_trial10)

print(mattick2_trial10)

fit_trial10 <- cbind(time.seq, fitted(mattick2_trial10, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial10$trial <- as.factor(rep(10, length(nrow(time.seq))))

predline10 <-
  predict(mattick2_trial10, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline10 <- ggplot(data = trial10, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline10, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial10, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 10")


# trial 11:
trial11 <- filter(mattick65, trial==11)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial11 <- brm(formula=nlform, data=trial11, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial11"))

post_trial11 <- posterior_samples(mattick2_trial11)
#
print(mattick2_trial11)

fit_trial11 <- cbind(time.seq, fitted(mattick2_trial11, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial11$trial <- as.factor(rep(11, length(nrow(time.seq))))

predline11 <-
 predict(mattick2_trial11, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline11 <- ggplot(data = trial11, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline11, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial11, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 11")


#trial 12:
trial12 <- filter(mattick65, trial==12)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(5,5), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial12 <- brm(formula=nlform, data=trial12, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial12"))

post_trial12 <- posterior_samples(mattick2_trial12)

print(mattick2_trial12)

fit_trial12 <- cbind(time.seq, fitted(mattick2_trial12, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial12$trial <- as.factor(rep(12, length(nrow(time.seq))))

#predline12 <-
#  predict(mattick2_trial12, 
#          newdata = time.seq) %>%
#  as_tibble() %>%
#  bind_cols(time.seq)

#saveRDS(predline12, file=here("fits","predline12"))
predline12 <- readRDS(file=here("fits", "predline12"))


regrline12 <- ggplot(data = trial12, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline12, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial12, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 12")

#trial 13:
trial13 <- filter(mattick65, trial==13)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial13 <- brm(formula=nlform, data=trial13, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial13"))

post_trial13 <- posterior_samples(mattick2_trial13)

print(mattick2_trial13)

fit_trial13 <- cbind(time.seq, fitted(mattick2_trial13, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial13$trial <- as.factor(rep(13, length(nrow(time.seq))))

predline13 <-
  predict(mattick2_trial13, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline13 <- ggplot(data = trial13, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline13, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial13, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 13")

# trial 14:
trial14 <- filter(mattick65, trial==14)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick2_trial14 <- brm(formula=nlform, data=trial14, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial14"))

post_trial14 <- posterior_samples(mattick2_trial14)
print(mattick2_trial14)
fit_trial14 <- cbind(time.seq, fitted(mattick2_trial14, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial14$trial <- as.factor(rep(14, length(nrow(time.seq))))

predline14 <-
  predict(mattick2_trial14, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline14 <- ggplot(data = trial14, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline14, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial14, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 14")


# trial 15:
trial15 <- filter(mattick65, trial==15)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial15 <- brm(formula=nlform, data=trial15, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial15"))

post_trial15 <- posterior_samples(mattick2_trial15)

print(mattick2_trial15)

fit_trial15 <- cbind(time.seq, fitted(mattick2_trial15, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial15$trial <- as.factor(rep(15, length(nrow(time.seq))))

predline15 <-
  predict(mattick2_trial15, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline15 <- ggplot(data = trial15, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline15, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial15, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 15")


# trial 16:
trial16 <- filter(mattick65, trial==16)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)


mattick2_trial16 <- brm(formula=nlform, data=trial16, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial16"))

post_trial16 <- posterior_samples(mattick2_trial16)
print(mattick2_trial16)

fit_trial16 <- cbind(time.seq, fitted(mattick2_trial16, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial16$trial <- as.factor(rep(16, length(nrow(time.seq))))

predline16 <-
  predict(mattick2_trial16, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline16 <- ggplot(data = trial16, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline16, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial16, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 16")

# trial 17:
trial17 <- filter(mattick65, trial==17)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick2_trial17 <- brm(formula=nlform, data=trial17, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial17"))

post_trial17 <- posterior_samples(mattick2_trial17)

print(mattick2_trial17)

fit_trial17 <- cbind(time.seq, fitted(mattick2_trial17, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial17$trial <- as.factor(rep(17, length(nrow(time.seq))))

predline17 <-
  predict(mattick2_trial17, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline17 <- ggplot(data = trial17, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline17, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial17, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 17")

# trial 18:
trial18 <- filter(mattick65, trial==18)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta"),
           prior(cauchy(0,25), class="sigma")
)

mattick2_trial18 <- brm(formula=nlform, data=trial18, family = gaussian(), prior = nlprior, chains= 4, warmup=2000, iter = 4000, control = list(adapt_delta = 0.999), file=here("fits", "mattick2_trial18"))

post_trial18 <- posterior_samples(mattick2_trial18)

print(mattick2_trial18)

fit_trial18 <- cbind(time.seq, fitted(mattick2_trial18, newdata=time.seq, re_formula = NA))[,-(3:5)]

fit_trial18$trial <- as.factor(rep(18, length(nrow(time.seq))))

predline18 <-
  predict(mattick2_trial18, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline18 <- ggplot(data = trial18, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = predline18, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_trial18, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "trial 18")


# collection of unpooled individual fits:
fit_unpooled_all <- rbind(fit_trial1, fit_trial2, fit_trial3, fit_trial4, fit_trial5, fit_trial6, fit_trial7, fit_trial8, fit_trial9, fit_trial10, fit_trial11, fit_trial12, fit_trial13, fit_trial14, fit_trial15, fit_trial16, fit_trial17, fit_trial18)

# printing the individual fits (Figure S4 in Supplement):
regrline1 + regrline2 + regrline3 + regrline4 +  regrline5 + regrline6 +regrline7 + regrline8+regrline9 + regrline10 +regrline11 +regrline12 +regrline13 + regrline14 +regrline15 +regrline16 +regrline17 +regrline18 +plot_layout(ncol=6)

#Forest plots for the no-pooled case (Figure 5 in main article):
p_alpha <-
  bind_rows(
    post_trial1,
    post_trial2,
    post_trial3,
    post_trial4,
    post_trial5,
    post_trial6,
    post_trial7,
    post_trial8,
    post_trial9,
    post_trial10,
    post_trial11,
    post_trial12,
    post_trial13,
    post_trial14,
    post_trial15,
    post_trial16,
    post_trial17,
    post_trial18
  )
iter <- 8000

p_alpha <- 
  p_alpha %>% 
  mutate(trial = rep(c("1","2","3","4", "5", "6", "7", "8", "9", "10",  "11", "12",  "13", "14",  "15", "16", "17", "18"), each = iter)) %>%   mutate(serial=fct_relevel(trial, "1","2","3","4", "5", "6", "7", "8", "9", "10",  "11", "12",  "13", "14",  "15", "16", "17", "18"))

alpha_plot_unpooled <- p_alpha %>% 
  ggplot(aes(x = b_alpha_Intercept, y = serial)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(alpha), y="trial")+
  coord_cartesian(xlim =c(0, 2.3))

logN0_plot_unpooled <- p_alpha %>% 
  ggplot(aes(x = b_logN0_Intercept, y = serial)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("log"[10],"N"[0])), y="trial")+
  coord_cartesian(xlim =c(5.5, 8.5))

beta_plot_unpooled <- p_alpha %>% 
  ggplot(aes(x = b_beta_Intercept, y = serial)) +
  geom_halfeyeh(fill = "green4", alpha=0.5,
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(beta), y="trial")+
  coord_cartesian(xlim =c(0, 0.13))

logN0_plot_unpooled + alpha_plot_unpooled + beta_plot_unpooled

# brms code for multilevel regression (partial pooling), case study 1
nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1+(1|ID|trial), 
           alpha~1+(1|ID|trial), 
           beta~1+(1|ID|trial),
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,0.5), nlpar="alpha"),
           prior(normal(0.03,0.01), nlpar="beta"),
           prior(cauchy(0,10), class="sd", nlpar="logN0"),
           prior(cauchy(0,10), class="sd", nlpar="alpha"),
           prior(cauchy(0,10), class="sd", nlpar="beta"),
           prior(cauchy(0,10), class="sigma"),
           prior(lkj(1), class="cor")
)

mattick2_multi <- brm(formula=nlform, data=mattick65, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999, max_treedepth=15), file=here("fits", "mattick2_multi"))

mattick2_post <- posterior_samples(mattick2_multi)

print(mattick2_multi, digits=4)

# pair and correlation plots (Figure S5 in Supplement)
mattick_cor2 <- dplyr::select(mattick2_post,b_logN0_Intercept:sigma)

mattick_cor2 <- setNames(mattick_cor2, c(expression(paste("log"[10],"N")), expression(alpha), expression(beta), expression(sigma[logN[0]]),expression(sigma[alpha]),expression(sigma[beta]),expression(rho[logN[0]-alpha]),expression(rho[logN[0]-beta]),expression(rho[alpha-beta]), expression(sigma[e])))


(mattick_corplot2 <-mattick_cor2  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 2, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))  +
    theme(axis.text.x = element_text(angle=70, hjust=1, size = 6))+
    theme(axis.text.y = element_text(angle=0, hjust=1, size = 6))
) 

#Code for Table 2 in the manuscript:

mattick2_multi_cor <- summary(mattick2_multi)

mattick2_multi_cor2 <- rbind(data.frame(mattick2_multi_cor$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(do.call(rbind,mattick2_multi_cor$random))%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(mattick2_multi_cor$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

rownames(mattick2_multi_cor2) <- c("$\\log N_0$", "$\\alpha$", "$\\beta  (h)$", "$\\sigma_u$","$\\sigma_w$", "$\\sigma_v$","$\\rho_{uw}$", "$\\rho_{uv}$","$\\rho_{vw}$","$\\sigma_{e}$")

colnames(mattick2_multi_cor2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  mattick2_multi_cor2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:num-summary-multi)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,3,2,2)
  ),
  escape = FALSE
)

# Code for the individual deviations from the population variables with mcmc_plot from bayesplot (Figure 6 in main article):

re_plot1 <- mcmc_plot(mattick2_multi, pars=c("^r_trial__logN0"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("trial 1","trial 2","trial 3","trial 4","trial 5","trial 6","trial 7","trial 8","trial 9","trial 10","trial 11","trial 12","trial 13","trial 14","trial 15","trial 16","trial 17","trial 18")) +labs(x=expression(paste("u (log"[10],"N"[0],")")))
re_plot2 <- mcmc_plot(mattick2_multi, pars=c("^r_trial__alpha"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("trial 1","trial 2","trial 3","trial 4","trial 5","trial 6","trial 7","trial 8","trial 9","trial 10","trial 11","trial 12","trial 13","trial 14","trial 15","trial 16","trial 17","trial 18")) + labs(x=expression(paste("w (",alpha,")")))
re_plot3 <- mcmc_plot(mattick2_multi, pars=c("^r_trial__beta"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c("trial 1","trial 2","trial 3","trial 4","trial 5","trial 6","trial 7","trial 8","trial 9","trial 10","trial 11","trial 12","trial 13","trial 14","trial 15","trial 16","trial 17","trial 18")) + labs(x=expression(paste("v (",beta,")")))
re_plot1+re_plot2+re_plot3

# Code for plotting the fits from the multilevel results, case study 1 (Figure 7 in main article):

newavg <- data.frame(time = seq(from = 0, to = 0.3, by = 0.005))

# calculating the credible interval

fitavg2 <- cbind(newavg, fitted(mattick2_multi, newdata = newavg, re_formula = NA)[,-2])
# store and read back the file if desired: 

#saveRDS(fitavg2, file=here("fits","fitavg2_multi"))
#fitavg2 <- readRDS(file=here("fits", "fitavg2_multi"))

names(fitavg2) <- c("time", "logN", "lower", "upper")

# same for the prediction interval
predavg2 <- cbind(newavg, predict(mattick2_multi, newdata = newavg, re_formula = NA)[,-2])

#saveRDS(predavg2, file=here("fits","predavg2_multi"))
#predavg2 <- readRDS(file=here("fits", "predavg2_multi"))
names(predavg2) <- c("time", "logN", "lower", "upper")

# Building up the plot with the multilevel fits:
mattick_plot6 <- mattick_plot2 + 
  geom_ribbon(data = fitavg2, aes(ymin = lower, ymax = upper), fill = "red", alpha =0.3)+
  labs(x="time (h)", y=expression(paste("log"[10],"N")), subtitle = "A: 95% CI")

mattick_plot7 <- mattick_plot3 +
  geom_ribbon(data = predavg2, aes(ymin = lower, ymax = upper), fill = "red", alpha =0.3)+
  labs(x="time (h)", y=expression(paste("log"[10],"N")), subtitle = "B: 95% PI")

mattick_plot6 + mattick_plot7

# The leave-one-out-cross-validation procedyre :
mattick2_multi <- add_criterion(mattick2_multi, c("loo"), file=here("fits","mattick2_multi"))
mattick_pooled2 <- add_criterion(mattick_pooled2, c("loo"), file=here("fits", "mattick_pooled2"))

loo_multi_mattick <- loo(mattick2_multi)
loo_pooled_mattick <- loo(mattick_pooled2)

loo_result_mattick <- loo_compare(loo_multi_mattick,loo_pooled_mattick)

elpd_diff_mattick <- c(loo_result_mattick[1,1], loo_result_mattick[2,1])
se_diff_mattick <- c(loo_result_mattick[1,2], loo_result_mattick[2,2])
loo_df_mattick <- data.frame(elpd_diff_mattick, se_diff_mattick)
colnames(loo_df_mattick) <- c("loo-cv-value","SE" )
rownames(loo_df_mattick) <- c("multilevel partially pooled","single level completely pooled")
# Putting the results in an apa-table:
apa_table(
  loo_df_mattick,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:loo-mattick)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(1,1)
  ),
  escape = FALSE)

# Partial pooling fits (Figure S6 in Supplement)
time.seq <- data.frame(time = seq(from = 0, to = 0.3, by = 0.005))

muSummary <-
  fitted(mattick_pooled2, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# regression line + data points
mattick_plot1 <- mattick65 %>%
  ggplot(aes(x = time, y = logN)) +
  geom_line(data = muSummary, aes(y = Estimate)) +
  geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3) +
  labs(x="time (h)", y=expression(paste(log[10],"N")))

mattick_plot8 <- mattick_plot1 + facet_wrap(~trial, ncol=5)

# collection of partially pooled fits:
fit_partpooled_all <- cbind(mattick65, fitted(mattick2_multi)[,-2])

# individual plots from pooled (black), partpooled (red) and not-pooled (blue):
(mattick_plot9 <- mattick_plot8 +  
    geom_line(data=fit_partpooled_all, aes(y=Estimate, colour="red"))+
    geom_line(data=fit_unpooled_all, aes(y=Estimate, colour="blue"))+
    theme(legend.position = "none")
)

# Partial pooling fits for trials 1 and 4 (Figure S7 in Supplement)
tmp <- filter(mattick65, trial %in% c(1, 4))
tmp1 <- filter(fit_partpooled_all, trial %in% c(1,4))
tmp2 <- filter(fit_unpooled_all, trial %in% c(1,4))
(mattick_plot8 %+% tmp +
    geom_line(data=tmp1, lty=2, size=1,aes(y=Estimate, colour="blue"))+
    geom_line(data=tmp2, lty=3,size=1, aes(y=Estimate, colour="blue"))+
    theme(legend.position = "none")
)

#Posterior predictive checks (Figure S8 in Supplement):
plot_ppc_A <- pp_check(mattick2_multi, nsamples=100)+labs(x="log N", subtitle = "A")
plot_ppc_B <- pp_check(mattick_pooled2, nsamples=100)+labs(x="log N", subtitle ="B")
plot_ppc_A + plot_ppc_B



# CASE STUDY 2

# Overview of the data (Figure 8 in the manuscript)
salmo %>%  mutate(temp=str_c("T=",temp)) %>% ggplot(aes(x=time, y=logN, color=trial))+
  geom_point(aes(y=logN))+
  geom_line()+
  facet_wrap(~temp, ncol=3, scales = "free")+
  labs(x = "time (h)", y=expression(paste(log[10],"N")))+
  theme(legend.position = "none")

# single level modeling at each temperature using the Weibull model with logN0 as parameter

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.01,0.01), nlpar="beta")
)


Salmo2_model55 <- brm(formula=nlform, data=salmo55, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits","Salmo2_model55"))

print(Salmo2_model55)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

Salmo2_model60 <- brm(formula=nlform, data=salmo60, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model60"))

print(Salmo2_model60)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

Salmo2_model65 <- brm(formula=nlform, data=salmo65, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits","Salmo2_model65"))

print(Salmo2_model65, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.05,0.05), nlpar="beta")
)

Salmo2_model70 <- brm(formula=nlform, data=salmo70, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model70"))

print(Salmo2_model70, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)


Salmo2_model72 <- brm(formula=nlform, data=salmo72, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model72"))

print(Salmo2_model65, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.1,0.1), nlpar="beta")
)

Salmo2_model74 <- brm(formula=nlform, data=salmo74, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model74"))

print(Salmo2_model74, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)


Salmo2_model76 <- brm(formula=nlform, data=salmo76, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model76"))

print(Salmo2_model76, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(1,1), nlpar="alpha"),
           prior(normal(0.005,0.001), nlpar="beta")
)

Salmo2_model78 <- brm(formula=nlform, data=salmo78, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model78"))

print(Salmo2_model78, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^alpha, 
           logN0~1, 
           alpha~1, 
           beta~1,
           nl=TRUE)

Salmo2_model80 <- brm(formula=nlform, data=salmo80, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model80"))

print(Salmo2_model80, digits=4)

# Fits of individual regression lines (Figure S9 in manuscript)
# composing a dataframe containing fitted values for the individual regressions:
# combines the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =1.5, by=0.1),temp=c(55))
newvary2 <- expand.grid(time=seq(from = 0, to =1.05, by=0.1),temp=c(60))
newvary3 <- expand.grid(time=seq(from = 0, to =0.3, by=0.01),temp=c(65))
newvary4 <- expand.grid(time=seq(from = 0, to =0.15, by=0.005),temp=c(70))
newvary5 <- expand.grid(time=seq(from = 0, to =0.08, by=0.004),temp=c(72))
newvary6 <- expand.grid(time=seq(from = 0, to =0.05, by=0.001),temp=c(74))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(76))
newvary8 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(78))
newvary9 <- expand.grid(time=seq(from = 0, to =0.02, by=0.0001),temp=c(80))


fit_55 <- cbind(newvary1, fitted(Salmo2_model55, newdata=newvary1, re_formula = NA))[,-4]

regrline55 <- ggplot(data = salmo55, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_55, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_55, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "55 C")

fit_60 <- cbind(newvary2, fitted(Salmo2_model60, newdata=newvary2, re_formula = NA))[,-4]

regrline60 <- ggplot(data = salmo60, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_60, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_60, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "60 C")

fit_65 <- cbind(newvary3, fitted(Salmo2_model65, newdata=newvary3, re_formula = NA))[,-4]

regrline65 <- ggplot(data = salmo65, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_65, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_65, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "65 C")

fit_70 <- cbind(newvary4, fitted(Salmo2_model70, newdata=newvary4, re_formula = NA))[,-4]

regrline70 <- ggplot(data = salmo70, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_70, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_70, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "70 C")

fit_72 <- cbind(newvary5, fitted(Salmo2_model72, newdata=newvary5, re_formula = NA))[,-4]

regrline72 <- ggplot(data = salmo72, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_72, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_72, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "72 C")

fit_74 <- cbind(newvary6, fitted(Salmo2_model74, newdata=newvary6, re_formula = NA))[,-4]

regrline74 <- ggplot(data = salmo74, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_74, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_74, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "74 C")

fit_76 <- cbind(newvary7, fitted(Salmo2_model76, newdata=newvary7, re_formula = NA))[,-4]

regrline76 <- ggplot(data = salmo76, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_76, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_76, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "76 C")

fit_78 <- cbind(newvary8, fitted(Salmo2_model78, newdata=newvary8, re_formula = NA))[,-4]

regrline78 <- ggplot(data = salmo78, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_78, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_78, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "78 C")

fit_80 <- cbind(newvary9, fitted(Salmo2_model80, newdata=newvary9, re_formula = NA))[,-4]

regrline80 <- ggplot(data = salmo80, 
                     aes(x = time, y = logN)) +
  geom_ribbon(data = fit_80, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_line(data = fit_80, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y=expression(paste(log[10],"N")), subtitle = "80 C")

regrline55+regrline60+regrline65+regrline70+regrline72+regrline74+regrline76+regrline78+regrline80+plot_layout(ncol=3)

# Example of a pairs and correlation plot for T=55 C (Figure S10 in Supplement)
Salmo2_model55 <- brm(formula=nlform, data=salmo55, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model55"))
post4 <- posterior_samples(Salmo2_model55)
mattick_cor55 <- dplyr::select(post4,b_logN0_Intercept:sigma)
mattick_cor55 <- setNames(mattick_cor55, c(expression(paste("log"[10],"N"[0])), 
                                           expression(alpha),expression(beta), expression(sigma[e])))

(corplot55 <-mattick_cor55  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")        
    )) 

# Fits ar each temperature resulting from completely pooled data (Figure S11 in Supplement)
nlform<-bf(logN ~ logN0-(1/2.303)*(time/(10^bref*10^(-Zinv*(temp-67.25))))^alpha,
           logN0~1, 
           Zinv~1, 
           bref~1,
           alpha~1, 
           nl=TRUE)
nlprior<-c(prior(normal(6.5,0.01), nlpar = "logN0"),
           prior(normal(0.7,0.1), nlpar="alpha"),
           prior(normal(0.1,0.1), nlpar="c1"),
           prior(normal(-2.5,1), nlpar="bref")
)

Salmo2_model_all_pooled <- brm(formula=nlform, data=salmo, family = gaussian(), prior = nlprior, warmup=4000, iter=8000, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_model_all_pooled"))
Salmo2_model_all_pooled_post <- posterior_samples(Salmo2_model_all_pooled)

# combines the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =1.5, by=0.1),temp=c(55))
newvary2 <- expand.grid(time=seq(from = 0, to =1.05, by=0.1),temp=c(60))
newvary3 <- expand.grid(time=seq(from = 0, to =0.3, by=0.01),temp=c(65))
newvary4 <- expand.grid(time=seq(from = 0, to =0.15, by=0.005),temp=c(70))
newvary5 <- expand.grid(time=seq(from = 0, to =0.08, by=0.004),temp=c(72))
newvary6 <- expand.grid(time=seq(from = 0, to =0.05, by=0.001),temp=c(74))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(76))
newvary8 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(78))
newvary9 <- expand.grid(time=seq(from = 0, to =0.02, by=0.0001),temp=c(80))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6, newvary7,newvary8, newvary9)

#Salmo2_pooled <- cbind(newvary, predict(Salmo2_model_all_pooled, newdata=newvary)[,-2])
#saveRDS(Salmo2_pooled, file=here("fits", "Salmo2_pooled.rds"))

Salmo2_pooled=readRDS(file=here("fits", "Salmo2_pooled.rds"))

names(Salmo2_pooled) <- c("time", "temp", "logN", "lower", "upper")

Salmo2_pooled$temp=as.factor(Salmo2_pooled$temp)

(pooled_fit <-  ggplot(salmo, aes(x=time, y=logN))+
    geom_point(shape=21, fill="blue", color="black")+
    facet_wrap(~temp, ncol=5, scales="free_x")+
    geom_line(data = Salmo2_pooled, aes(y = logN), size = 1, colour="blue") +
    geom_line(data = Salmo2_pooled, aes(y = lower), lty = 2) +
    geom_line(data = Salmo2_pooled, aes(y = upper), lty = 2)) +
  labs(x="time (h)", y=expression(paste(log[10],"N")))

# Pair plots and correlation for the completely pooled data (Figure S12 in Supplement)
pooled_cor1 <- dplyr::select(Salmo2_model_all_pooled_post,b_logN0_Intercept:sigma) %>% mutate(beta_ref=10^b_bref_Intercept) %>% mutate(Z=1/b_c1_Intercept) %>% mutate(sigma_e=sigma)

pooled_cor1a <- pooled_cor1 %>% select(-b_bref_Intercept,-b_c1_Intercept,-sigma)

pooled_cor1a <- setNames(pooled_cor1a, c(expression(paste("log"[10],"N"[0])), 
                                         expression(alpha),expression(beta[ref]),"Z", expression(sigma[e])))


(pooled_corplot1a <-pooled_cor1a  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")  +
            theme(axis.text.x = element_text(angle=70, hjust=1, size=6))
    ))


# collection of posterior results and plotting:
post1 <- posterior_samples(Salmo2_model55)
post2 <- posterior_samples(Salmo2_model60)
post3 <- posterior_samples(Salmo2_model65)
post4 <- posterior_samples(Salmo2_model70)
post5 <- posterior_samples(Salmo2_model72)
post6 <- posterior_samples(Salmo2_model74)
post7 <- posterior_samples(Salmo2_model76)
post8 <- posterior_samples(Salmo2_model78)
post9 <- posterior_samples(Salmo2_model80)


p <-
  bind_rows(
    post1,
    post2,
    post3,
    post4,
    post5,
    post6,
    post7,
    post8,
    post9
  )
iter <- 4000

p <- 
  p %>% 
  mutate(temperature = rep(c("55 C","60 C","65 C","70 C","72 C","74 C","76 C","78 C", "80 C"),each = iter))
nt_plot <- p %>% 
  ggplot(aes(x = b_alpha_Intercept, y = temperature)) +
  geom_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95)+
  labs(x=expression(alpha), subtitle = "B")

logN0_plot <- p %>% 
  ggplot(aes(x = b_logN0_Intercept, y = temperature)) +
  geom_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("log"[10],"N"[0])), subtitle = "A")

logkr_plot <- p %>% 
  ggplot(aes(x = log10(b_beta_Intercept), y = temperature)) +
  geom_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95) +
  labs(x=expression(paste("log"[10],beta,)), subtitle = "C")

logN0_plot + nt_plot + logkr_plot

# Calculation of reference temperature Tref:

salmo <- salmo %>% mutate(nominator=(logN*log(logN))^2) %>% mutate(denominator=(logN*log(logN))^2/temp)
Tref <- sum(salmo$nominator)/sum(salmo$denominator)

# Bayesian regression with use of the Bigelow model on the pooled data, analysis with Tref=67.25 C: 

nlform<-bf(logN ~ logN0-(1/2.303)*(time/(10^bref*10^(-c1*(temp-67.25))))^alpha,
           logN0~1, 
           c1~1, 
           bref~1,
           alpha~1, 
           nl=TRUE)
nlprior<-c(prior(normal(6.5,0.01), nlpar = "logN0"),
           prior(normal(0.6,0.06), nlpar="alpha"),
           prior(normal(0.09,0.009), nlpar="c1"),
           prior(normal(-2.5,1), nlpar="bref")
)

Salmo3_model_all_pooled <- brm(formula=nlform, data=salmo, family = gaussian(), prior = nlprior, warmup=4000, iter=8000, control = list(adapt_delta = 0.999), file=here("fits", "Salmo3_model_all_pooled"))
Salmo3_model_all_pooled_post <- posterior_samples(Salmo3_model_all_pooled)
print(Salmo3_model_all_pooled, digits=4)

# Putting the Bigelow pooled results in a Table (Table 4 in manuscript):
pooled_summary1 <- summary(Salmo3_model_all_pooled)
pooled_summary2 <- rbind(data.frame(pooled_summary1$fixed)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(pooled_summary1$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

pooled_cor1 <- dplyr::select(Salmo3_model_all_pooled_post,b_logN0_Intercept:sigma) %>% mutate(beta_ref=10^b_bref_Intercept) %>% mutate(Z=1/b_c1_Intercept) %>% mutate(sigma_e=sigma)

pooled_cor1a <- pooled_cor1 %>% select(-b_bref_Intercept,-b_c1_Intercept,-sigma)

pooled_summary2[2,1]=mean(pooled_cor1a$Z)
pooled_summary2[2,2]=sd(pooled_cor1a$Z)
pooled_summary2[2,3]=quantile(pooled_cor1a$Z, probs=0.025)
pooled_summary2[2,4]=quantile(pooled_cor1a$Z, probs=0.975)

pooled_summary2[3,1]=mean(pooled_cor1a$beta_ref)
pooled_summary2[3,2]=sd(pooled_cor1a$beta_ref)
pooled_summary2[3,3]=quantile(pooled_cor1a$beta_ref, probs=0.025)
pooled_summary2[3,4]=quantile(pooled_cor1a$beta_ref, probs=0.975)

rownames(pooled_summary2) <- c("$\\log N_0$", "$Z (^oC)$", "$\\beta_{ref} (h)$", "$\\alpha (-)$", "$\\sigma_e$")
colnames(pooled_summary2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  pooled_summary2,
  placement = "H",
  align = c("c", "c", "c", "c", "c"),
  caption = "(ref:summary-pooled)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,3,2,2)
  ),
  escape = FALSE
)


# same as above but now with multilevel modeling where nt is allowed to vary per temperature AND per trial:
#Tref=67.25C:

nlform<-bf(logN ~ logN0-(1/2.303)*(time/((10^bref)*10^(Zinv*(67.25-temp))))^alpha,
           logN0~1, 
           Zinv~1, 
           bref~1,
           alpha~1+(1|temp/trial), 
           nl=TRUE)

nlprior<-c(prior(normal(6.5,0.001), nlpar = "logN0"),
           prior(normal(0.67,0.001), nlpar="alpha"),
           prior(normal(0.09,0.01), nlpar="Zinv"),
           prior(normal(-2.0,0.1), nlpar="bref"),
           prior(exponential(1), class="sd", nlpar="alpha"),
           prior(exponential(1), class="sigma")
)

inits <- list(logN0=6.5, alpha=0.63,Zinv=0.09,bref=-1.8, sigma=0.5)
inits_list <- list(inits,inits,inits,inits)

Salmo2_all_multi3 <- brm(formula=nlform, data=salmo, family = gaussian(), inits = 0, prior = nlprior,warmup=4000, iter=8000, control = list(adapt_delta = 0.99, max_treedepth=12), file=here("fits", "Salmo2_all_multi3"))

Salmo2_all_multi3_post <- posterior_samples(Salmo2_all_multi3)
print(Salmo2_all_multi3)

# Pair plots and correlations resulting from partial pooling using the Bigelow model (Figure S13 in Supplememt); multilevel regression where nt is allowed to vary per temperature AND per trial:

nlform<-bf(logN ~ logN0-(1/2.303)*(time/((10^bref)*10^(c1*(70-temp))))^alpha,
           logN0~1, 
           c1~1, 
           bref~1,
           alpha~1+(1|temp/trial), 
           nl=TRUE)

inits <- list(logN0=6.5, alpha=0.63,c1=0.09,bref=-2.5, sigma=0.5)
inits_list <- list(inits,inits,inits,inits)

nlprior<-c(prior(normal(6.5,0.01), nlpar = "logN0"),
           prior(normal(0.63,0.05), nlpar="alpha"),
           prior(normal(0.09,0.01), nlpar="c1"),
           prior(normal(-2.5,0.1), nlpar="bref"),
           prior(exponential(1), class="sd", nlpar="alpha"),
           prior(exponential(1), class="sigma")
)

Salmo2_all_multi2 <- brm(formula=nlform, data=salmo, family = gaussian(), prior = nlprior, inits=inits_list,warmup=2000, iter=4000, control = list(adapt_delta = 0.9, max_treedepth=12), file=here("fits", "Salmo2_all_multi2"))

Salmo2_all_multi2_post <- posterior_samples(Salmo2_all_multi2)
print(Salmo2_all_multi2)

multi_cor2 <- dplyr::select(Salmo2_all_multi2_post,b_logN0_Intercept:sigma)%>% mutate(beta_ref=10^b_bref_Intercept) %>% mutate(Z=1/b_c1_Intercept) %>% mutate(sigma_e=sigma)

multi_cor2a <- multi_cor2 %>% select(-b_bref_Intercept,-b_c1_Intercept,-sigma)

multi_cor2a <- setNames(multi_cor2a, c(expression(paste("log"[10],"N"[0])), 
                                       expression(alpha),expression(sigma[alpha][-T]),expression(sigma[alpha][-trial]),expression(beta[ref]),"Z", expression(sigma[e])))

(multi_corplot2 <-multi_cor2a %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 3, stars=FALSE, title="corr. coef")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red"))  +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=6)
    )) 

# Plot of random deviations from the population level

re_plot1 <- mcmc_plot(Salmo2_all_multi3, pars=c("^r_temp__"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c(expression(paste(alpha, " , T=55")),expression(paste(alpha, " , T=60")),expression(paste(alpha, " , T=65")),expression(paste(alpha, " , T=70")),expression(paste(alpha, " , T=72")),expression(paste(alpha, " , T=74")),expression(paste(alpha, " , T=76")) ,expression(paste(alpha, " , T=78")), expression(paste(alpha," , T= 80"))))
re_plot2 <- mcmc_plot(Salmo2_all_multi3, pars=c("^r_temp:trial__"), point_est="mean", prob_outer=0.95, prob=0.5) +
  ggplot2::scale_y_discrete(labels=c(expression(paste(alpha, " , Trial 55.1")),expression(paste(alpha, " , Trial 55.2")),expression(paste(alpha, " , Trial 55.3")),expression(paste(alpha, " , Trial 60.1")),expression(paste(alpha, " , Trial 60.2")),expression(paste(alpha, " , Trial 60.3")),expression(paste(alpha, " , Trial 60.4")) ,expression(paste(alpha, " , Trial 65.1")), expression(paste(alpha," , Trial 65.2")),expression(paste(alpha, " , Trial 65.3")),expression(paste(alpha, " , Trial 65.4")),expression(paste(alpha, " , Trial 70.1")),expression(paste(alpha, " , Trial 70.2")),expression(paste(alpha, " , Trial 70.3")),expression(paste(alpha, " , Trial 72.1")),expression(paste(alpha, " , Trial 72.2")),expression(paste(alpha, " , Trial 72.3")),expression(paste(alpha, " , Trial 72.4")),expression(paste(alpha, " , Trial 74.1")),expression(paste(alpha, " , Trial 74.2")),expression(paste(alpha, " , Trial 74.3")),expression(paste(alpha, " , Trial 76.1")),expression(paste(alpha, " , Trial 76.2")),expression(paste(alpha, " , Trial 76.3")),expression(paste(alpha, " , Trial 78.1")),expression(paste(alpha, " , Trial 78.2")),expression(paste(alpha, " , Trial 78.3")),expression(paste(alpha, " , Trial 78.4")),expression(paste(alpha, " , Trial 80.1")),expression(paste(alpha, " , Trial 80.2")),expression(paste(alpha, " , Trial 80.3")),expression(paste(alpha, " , Trial 80.4"))))
re_plot1+re_plot2

# Plot of the fits resulting from multilevel modeling (Figure 11 in manuscript):

# combines the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =1.5, by=0.1),temp=c(55), trial=c(55.1,55.2,55.3))
newvary2 <- expand.grid(time=seq(from = 0, to =1.05, by=0.1),temp=c(60), trial=c(60.1,60.2,60.3,60.4))
newvary3 <- expand.grid(time=seq(from = 0, to =0.3, by=0.01),temp=c(65), trial=c(65.1,65.2,65.3,65.4))
newvary4 <- expand.grid(time=seq(from = 0, to =0.15, by=0.005),temp=c(70), trial=c(70.1,70.2,70.3))
newvary5 <- expand.grid(time=seq(from = 0, to =0.08, by=0.004),temp=c(72), trial=c(72.1,72.2,72.3,72.4))
newvary6 <- expand.grid(time=seq(from = 0, to =0.05, by=0.001),temp=c(74), trial=c(74.1,74.2,74.3))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(76), trial=c(76.1,76.2,76.3))
newvary8 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(78), trial=c(78.1,78.2,78.3))
newvary9 <- expand.grid(time=seq(from = 0, to =0.02, by=0.0001),temp=c(80), trial=c(80.1,80.2,80.3,80.4))
newvary_multi2 <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6, newvary7,newvary8, newvary9)

Salmo2_multi2_NA <- cbind(newvary_multi2, predict(Salmo2_all_multi2, newdata=newvary_multi2, re_formula = NA)[,-2])

names(Salmo2_multi2_NA) <- c("time", "temp", "trial", "logN", "lower", "upper")

Salmo2_multi2_NA$temp=as.factor(Salmo2_multi2_NA$temp)
Salmo2_multi2_NA$trial=as.factor(Salmo2_multi2_NA$trial)

(ind_fit_multi2_NA <-  ggplot(salmo, aes(x=time, y=logN))+
    geom_point(shape=1)+
    facet_wrap(~temp, ncol=5, scales="free_x")+
    geom_line(data = Salmo2_multi2_NA, aes(y = logN), size = 1, colour="blue") +
    geom_line(data = Salmo2_multi2_NA, aes(y = lower), lty = 2) +
    geom_line(data = Salmo2_multi2_NA, aes(y = upper), lty = 2) +
    labs(x="time (h)", y=expression(paste(log[10],"N")))
)

# Nunerical summaries of the posterior resulting from multilevel modeling (Table 5 in manuscript):
multi_summary1 <- summary(Salmo2_all_multi3)
multi_summary2 <- rbind(data.frame(multi_summary1$fixed)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(multi_summary1$random$temp)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS),data.frame(multi_summary1$random$`temp:trial`)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS),data.frame(multi_summary1$spec_pars) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

multi_cor2 <- dplyr::select(Salmo2_all_multi3_post,b_logN0_Intercept:sigma)%>% mutate(beta_ref=10^b_bref_Intercept) %>% mutate(Z=1/b_Zinv_Intercept) %>% mutate(sigma_e=sigma)

multi_cor2a <- multi_cor2 %>% select(-b_bref_Intercept,-b_Zinv_Intercept,-sigma)

multi_summary2[2,1]=mean(multi_cor2a$Z)
multi_summary2[2,2]=sd(multi_cor2a$Z)
multi_summary2[2,3]=quantile(multi_cor2a$Z, probs=0.025)
multi_summary2[2,4]=quantile(multi_cor2a$Z, probs=0.975)

multi_summary2[3,1]=mean(multi_cor2a$beta_ref)
multi_summary2[3,2]=sd(multi_cor2a$beta_ref)
multi_summary2[3,3]=quantile(multi_cor2a$beta_ref, probs=0.025)
multi_summary2[3,4]=quantile(multi_cor2a$beta_ref, probs=0.975)

rownames(multi_summary2) <- c("$\\log N_0$", "$Z (^oC)$", "$\\beta_{ref} (h)$", "$\\alpha (-)$","$\\sigma_{\\alpha - T}$", "$\\sigma_{\\alpha - trial}$","$\\sigma_e$")
colnames(multi_summary2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  multi_summary2,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:multi-parms)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,3,2,2)
  ),
  escape = FALSE
)

# Fits at the group level resulting from partial pooling (Figure S14 in Supplement), combines the data with predictions using the random effects with re_formula = NULL (Figure S14 in Supplement)
newvary1 <- expand.grid(time=seq(from = 0, to =1.5, by=0.1),temp=c(55), trial=c(55.1,55.2,55.3))
newvary2 <- expand.grid(time=seq(from = 0, to =1.05, by=0.1),temp=c(60), trial=c(60.1,60.2,60.3,60.4))
newvary3 <- expand.grid(time=seq(from = 0, to =0.3, by=0.01),temp=c(65), trial=c(65.1,65.2,65.3,65.4))
newvary4 <- expand.grid(time=seq(from = 0, to =0.15, by=0.005),temp=c(70), trial=c(70.1,70.2,70.3))
newvary5 <- expand.grid(time=seq(from = 0, to =0.08, by=0.004),temp=c(72), trial=c(72.1,72.2,72.3,72.4))
newvary6 <- expand.grid(time=seq(from = 0, to =0.05, by=0.001),temp=c(74), trial=c(74.1,74.2,74.3))
newvary7 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(76), trial=c(76.1,76.2,76.3))
newvary8 <- expand.grid(time=seq(from = 0, to =0.04, by=0.001),temp=c(78), trial=c(78.1,78.2,78.3))
newvary9 <- expand.grid(time=seq(from = 0, to =0.02, by=0.0001),temp=c(80), trial=c(80.1,80.2,80.3,80.4))
newvary_multi2 <- rbind(newvary1,newvary2,newvary3, newvary4, newvary5, newvary6, newvary7,newvary8, newvary9)

Salmo2_multi2_NULL=readRDS(file=here("fits", "Salmo2_multi2_NULL.rds"))

names(Salmo2_multi2_NULL) <- c("time", "temp", "trial", "logN", "lower", "upper")

Salmo2_multi2_NULL$temp=as.factor(Salmo2_multi2_NULL$temp)
Salmo2_multi2_NULL$trial=as.factor(Salmo2_multi2_NULL$trial)

(ind_fit_multi2_NULL <-  ggplot(salmo, aes(x=time, y=logN, group=trial))+
    geom_point(shape=1)+
    facet_wrap(~temp,ncol=5, scales="free_x")+
    geom_line(data = Salmo2_multi2_NULL, aes(y = logN), size = 1, colour="blue") +
    labs(x="time (h)", y=expression(paste(log[10],"N")))+
    coord_cartesian(ylim=c(0,7))
)

# Model comparison with loo (Table 6 in manuscript):
Salmo2_model_all_pooled <- add_criterion(Salmo2_model_all_pooled, c("loo"), file=here("fits","Salmo2_model_all_pooled"))
Salmo2_all_multi1 <- add_criterion(Salmo2_all_multi1, c("loo"), file=here("fits", "Salmo2_all_multi1"))
Salmo2_all_multi2 <- add_criterion(Salmo2_all_multi2, c("loo"), file=here("fits", "Salmo2_all_multi2"))

loo_pooled <- loo(Salmo3_model_all_pooled)
loo_multi3 <- loo(Salmo2_all_multi3)

loo_result <- loo_compare(loo_pooled,loo_multi3)

elpd_diff <- c(loo_result[1,1], loo_result[2,1])
se_diff <- c(loo_result[1,2], loo_result[2,2])
loo_df <- data.frame(elpd_diff, se_diff)
colnames(loo_df) <- c("loo-cv-value","SE" )
rownames(loo_df) <- c("multilevel partially pooled","single level completely pooled")

apa_table(
  loo_df,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:loo-values)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Calculation of prediction times, Weibull model with logN0 as parameter, but with alpha fixed at 0.67

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(0.01,0.01), nlpar="beta")
)

Salmo2_fixed55 <- brm(formula=nlform, data=salmo55, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits","Salmo2_fixed55"))

print(Salmo2_fixed55)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

Salmo2_fixed60 <- brm(formula=nlform, data=salmo60, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed60"))


print(Salmo2_model60)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

Salmo2_fixed65 <- brm(formula=nlform, data=salmo65, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits","Salmo2_fixed65"))

print(Salmo2_fixed65)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(0.05,0.05), nlpar="beta")
)

Salmo2_fixed70 <- brm(formula=nlform, data=salmo70, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed70"))

print(Salmo2_fixed70, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)


Salmo2_fixed72 <- brm(formula=nlform, data=salmo72, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed72"))

print(Salmo2_fixed65, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(0.1,0.1), nlpar="beta")
)

Salmo2_fixed74 <- brm(formula=nlform, data=salmo74, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed74"))


print(Salmo2_fixed74, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)


Salmo2_fixed76 <- brm(formula=nlform, data=salmo76, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed76"))

print(Salmo2_fixed76, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)

nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(0.005,0.001), nlpar="beta")
)

Salmo2_fixed78 <- brm(formula=nlform, data=salmo78, family = gaussian(), prior = nlprior, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed78"))

print(Salmo2_fixed78, digits=4)

nlform<-bf(logN ~ logN0-(1/2.303)*(time/beta)^0.67, 
           logN0~1, 
           beta~1,
           nl=TRUE)
nlprior<-c(prior(normal(7,1), nlpar = "logN0"),
           prior(normal(0.001,0.0001), nlpar="beta")
)

Salmo2_fixed80 <- brm(formula=nlform, data=salmo80, family = gaussian(), prior = nlprior, warmup=4000, iter=8000, control = list(adapt_delta = 0.999), file=here("fits", "Salmo2_fixed80"))

print(Salmo2_fixed80, digits=4)

# calculation of TDT with Tref=67.25

temp_T <- c(55-67.25,60-67.25,65-67.25,67.25-67.25,72-67.25,74-67.25,76-67.25,78-67.25,80-67.25)
log_beta <- c(log10(summary(Salmo2_fixed55)[["fixed"]][2,1]),log10(summary(Salmo2_fixed60)[["fixed"]][2,1]),log10(summary(Salmo2_fixed65)[["fixed"]][2,1]),log10(summary(Salmo2_fixed70)[["fixed"]][2,1]),log10(summary(Salmo2_fixed72)[["fixed"]][2,1]),log10(summary(Salmo2_fixed74)[["fixed"]][2,1]),log10(summary(Salmo2_fixed76)[["fixed"]][2,1]),log10(summary(Salmo2_fixed78)[["fixed"]][2,1]),log10(summary(Salmo2_fixed80)[["fixed"]][2,1]))

TDT_df <- data.frame(temp_T, log_beta)
TDT_fit <- 
  brm(data = TDT_df, family = gaussian,
      formula = log_beta ~ 1 + (temp_T),
      chains = 4, iter = 4000, warmup = 2000, 
      control = list(adapt_delta =0.95), 
      seed=15, file=here("fits", "TDT_fit"))
TDT_fit_post <- posterior_samples(TDT_fit)
print(TDT_fit, digits=4)

# Plot of the TDT (Figure S15 in Supplement):

(TDT_plot <- ggplot(TDT_df, aes(x=temp_T, y=log_beta))+
    geom_point()+
    geom_smooth(method = lm, se=TRUE)+
    labs(x=expression(T - T[ref]), y=expression(paste("log", beta)))
)

# Numerical summary of TDT regression (Table S2 in Supplement):

TDT_fit_post <- TDT_fit_post %>% mutate(beta_ref=10^b_Intercept) %>% mutate(Z=-1/b_temp_T)
TDT_beta <- c(mean(TDT_fit_post$beta_ref), sd(TDT_fit_post$beta_ref))
TDT_Z <- c(mean(TDT_fit_post$Z), sd(TDT_fit_post$Z))
TDT_summary <- data.frame(TDT_beta,TDT_Z)
TDT_summary <- t(TDT_summary)
colnames(TDT_summary) <- c("mean","SE" )
rownames(TDT_summary) <- c("$\\beta_{ref} (h)$","$Z (^o C)$")

apa_table(
  TDT_summary,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:TDT-summary)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(3,4)
  ),
  escape = FALSE
)

# Plotting of densities for predicted times (Figure 12 in manuscript)
Salmo2_all_multi3_post <- Salmo2_all_multi3_post %>% mutate(t6D75=6*(10^b_bref_Intercept)*10^(b_Zinv_Intercept*(67.25-75))*(log(10))^(1/b_alpha_Intercept))

TDT_fit_post <- TDT_fit_post %>% mutate(t6D75=6*(10^b_Intercept)*10^(-b_temp_T*(67.25-75))*(log(10))^(1/0.67)) 

(t75_density <- ggplot(data=Salmo2_all_multi3_post)+
    geom_density(fill="red", aes(x=t6D75))+
    geom_vline(xintercept=mean(Salmo2_all_multi3_post$t6D75), lty=2)+
    geom_density(data=TDT_fit_post, fill="turquoise", alpha=0.3, aes(x=t6D75))+
    geom_vline(xintercept=mean(TDT_fit_post$t6D75), lty=2)+
    geom_density(data=TDT_fit_post,fill="blue",alpha=0.1, aes(x=t6D75))+
    labs(x="time needed at 75 °C (h) for 6 decimal reductions")
)


