## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(digits = 3)

## ----load_pkg, message=FALSE--------------------------------------------------
library(rags2ridges)

## ----createS------------------------------------------------------------------
p <- 6
n <- 20
X <- createS(n = n, p = p, dataset = TRUE)
head(X, n = 4) # Show 4 first of the n rows

## ----covML--------------------------------------------------------------------
S <- covML(X)
print(S)

## ----solveS-------------------------------------------------------------------
P <- solve(S)
print(P)

## ----solveS2, error=TRUE------------------------------------------------------
p <- 25
S2 <- createS(n = n, p = p)  # Direct to S
P2 <- solve(S2)

## ----ridgeP-------------------------------------------------------------------
P2 <- ridgeP(S2, lambda = 1.17)
print(P2[1:7, 1:7]) # Showing only the 7 first cols and rows

## ----optPenalty.LOOCV---------------------------------------------------------
Y <- createS(n, p, dataset = TRUE)
opt <- optPenalty.kCVauto(Y, lambdaMin = 0.001, lambdaMax = 100)
str(opt)

