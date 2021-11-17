# clean up environment
rm(list = ls())
gc(reset = TRUE)

# set working directory
setwd("~/Desktop/Aalto/BDA/bda-project/R")

# install libraries
library(survival)
library(tidyverse)
library(loo)
library(ggplot2)
library(GGally)
# set number of cores
options(mc.cores = parallel::detectCores())

# read lung cancer data from `survival` library
data("cancer", package = "survival")
# build dataset
data <- cancer %>%
  drop_na()
# identify covariate labels
xtags <- data %>%
  dplyr::select(-status, -time) %>%
  colnames()
# investigate the underlying variate's distribution
data %>% 
  ggplot(aes(x=time)) +
  geom_density()
hist(data$time, breaks = 30)

# make transformed variate variable to measure linear correlation
data$logtime <- log(data$time)
cor.mat <- cor(data)
cor.mat["logtime",]

# split data into "died" and "censored", and measure correlations
# for the two independently
cens.data <- data[data$status == 1,]
cens.cor.mat <- cor(cens.data)
cens.cor.mat["logtime",]
died.data <- data[data$status == 2,]
died.data.mat <- cor(died.data)
died.data.mat["logtime",]
# we thus find that the correlation between some variables changes drastically
# depending on the value of the status column
# we thus must decide what we are looking to predict: time, status or time|status?

# plot some survival curves
plot(survfit(Surv(time, status) ~ 1, data = data))
plot(survfit(Surv(time, status) ~ age, data = data))
plot(survfit(Surv(time, status) ~ sex, data = data))
plot(survfit(Surv(time, status) ~ ph.ecog, data = data))

# pairs plot
GGally::ggpairs(data)
