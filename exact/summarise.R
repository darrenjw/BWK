#!/usr/bin/env Rscript
## summarise.R
## Summarise MCMC output

if (!require("pacman")) install.packages("pacman")
pacman::p_load("smfsb")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1)
    stop("Provide exactly one filename")

df = read.table(args[1], header=TRUE)
df = df[,-1]
mcmcSummary(df, bins=50)

df = df[,1:3]
mcmcSummary(df, bins=50)
image(cor(df)[ncol(df):1,])
pairs(df[sample(1:10000,1000),],pch=19,cex=0.2)

## eof

