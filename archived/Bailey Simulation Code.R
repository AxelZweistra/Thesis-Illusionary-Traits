############ Illusory traits simulation code ##############

# Prep
## random seed and packages 
set.seed(427)

library(mvtnorm)
library(psych)
library(ggplot2)
library(reshape)
library(gridExtra)
library(lavaan)
library(reshape2)

## Function definitions
row <- function( M ){
    if( ! ( is.matrix(M) && dim(M)[1] == dim(M)[2] && dim(M)[1]>0 ) ){
        rowM <- NULL
    } else {
        F <- dim(M)[1]
        rowM <- matrix( NA, nrow=F^2, ncol=1)
        for ( r in (1:F) ) {
            for( c in (1:F) ) {
                rowM[ c+F*(r-1), 1 ] <- M[r,c]
            }
        }
    }
    return( rowM )
}
irow <- function( rowM ){
    if( ! ( !is.null(rowM) && ( is.matrix(rowM) || (is.array(rowM) && length( which( dim(rowM) > 1 ) ) == 1) ) && dim(rowM)[1]>0 && dim(rowM)[2]==1 ) ){
        M <- NULL
    } else {
        F <- sqrt( length( rowM ) )
        M <- matrix( NA, nrow=F, ncol=F )
        for ( r in (1:F) ) {
            for ( c in (1:F) ) {
                M[r,c] <- rowM[ c+F*(r-1) ]
            }
        }
    }
    return( M )
}


# Set parameters for the data generation model
## Number of participants
N <- 20000 

## Number of generated time series
T <- 5

## Number of X-related confounders
c1 <- 5

## Number of Y-related confounders
c2 <- 5

## Autoregressive effects (for all variables)
alpha <- 0.3

## Cross-lagged effects of the focal effect (X and Y)
beta_focal <- 0.1

## Main confounding effects: Cross-lagged effects 
## (1) between X-related confounders and X
## (2) between Y-related confounders and Y
## (3) within X-related confounders
## (4) within Y-related confounders
## They are all set to be equal in this simulation but can be easily changed

beta_covfocal <- 0.1

## Trivial confounding effects: Cross-lagged effects 
## (1) between X-related confounders and Y
## (2) between Y-related confoudners and X
## (3) between X-related confounder and Y-related confounders
##  They are all set to be equal in this simulation but can be easily changed

beta_covnofocal <- 0.01

## Time-specific residual variance
sigma <- 0.5


# Create objects for the simulation

## total number of variables
v <- c1 + c2 + 2


## A: Matrix of path coefficients
##   A[1,2] - Effects of V2 at t-1 on V1 at t
##   A[2,1] - Effects of V1 at t-1 on V2 at t
##


A <- matrix(NA, v, v)
diag(A) <- alpha

A[lower.tri(A, diag = FALSE)] <- beta_covnofocal 
A[upper.tri(A, diag = FALSE)] <- beta_covnofocal

c1t <- 3:(2+c1)
c2t <- (3+c1):v

A[c1t, 1] <- beta_covfocal
A[1, c1t] <- beta_covfocal

A[c2t, 2] <- beta_covfocal
A[2, c2t] <- beta_covfocal


A[c1t, c1t][lower.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal
A[c1t, c1t][upper.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal

A[c2t, c2t][lower.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal
A[c2t, c2t][upper.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal

A[1,2] <- beta_focal
A[2,1] <- beta_focal




## Check the max eigenvalue (to check convergence --- see the manuscript)
# max(eigen(A)$values)

## Sigma: Var-cov matrix of time-specific residuals
Sigma <- diag(rep(sigma, v))

## Sigma1: Stationary Var-cov matrix for the first time point
Sigma1 <- irow( solve( diag(dim(A)[1]^2) - t(A) %x% t(A) ) %*% row(Sigma) )

## Mean1: Stationary means for first time point (as no process intercepts are defined/used, apparently 0)
Mean1 <- matrix( 0, nrow=dim(A)[1], ncol=1 )



# Generate data

## Wide-format "bivariate" (X and Y) data df: N x 2*T
## df stores the data from all time points
df <- matrix(NA, nrow = N, ncol = 2*T)
colnames(df) <- c(paste("x", 1:T, sep = ""), paste("y", 1:T, sep = ""))

## data at each time point D: N x v (total number of variables) 
## D stores the data for each time point only, i.e., D will be updated in a loop
## Generate the initial values for X and Y and store them in D
D <- rmvnorm(N, mean = Mean1, Sigma1)

## Store the initial data from D to df. 
## Store the initial X (the first column of D) and Y (the second column of D) in df
df[, 1] <- D[,1]
df[, 1+T] <- D[, 2]



## Update D and df till time T

for (i in 2:T){

D <-  D %*% t(A) + rmvnorm(N, sigma = Sigma)

df[, i] <- D[,1]
df[, i+T] <- D[, 2]

}



# Lavaan code
## CLPM

CLPM <- '
y1~~x1

x2~~ur*y2
x3~~ur*y3
x4~~ur*y4
x5~~ur*y5


x2~ax*x1
x3~ax*x2
x4~ax*x3
x5~ax*x4

y2~ay*y1
y3~ay*y2
y4~ay*y3
y5~ay*y4

x2~by*y1
x3~by*y2
x4~by*y3
x5~by*y4

y2~bx*x1
y3~bx*x2
y4~bx*x3
y5~bx*x4

'


## RI-CLPM

RICLPM <- '
RIx =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
RIy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5


wx1 =~ 1*x1
wx2 =~ 1*x2
wx3 =~ 1*x3 
wx4 =~ 1*x4
wx5 =~ 1*x5

wy1 =~ 1*y1
wy2 =~ 1*y2
wy3 =~ 1*y3
wy4 =~ 1*y4
wy5 =~ 1*y5

wx2 ~ ax*wx1 + by*wy1
wx3 ~ ax*wx2 + by*wy2
wx4 ~ ax*wx3 + by*wy3
wx5 ~ ax*wx4 + by*wy4

wy2 ~ bx*wx1 + ay*wy1
wy3 ~ bx*wx2 + ay*wy2
wy4 ~ bx*wx3 + ay*wy3
wy5 ~ bx*wx4 + ay*wy4

wx1 ~~ wy1 # Covariance

wx2 ~~ ur*wy2
wx3 ~~ ur*wy3
wx4 ~~ ur*wy4
wx5 ~~ ur*wy5

RIx ~~ varRIx*RIx
RIy ~~ varRIy*RIy 
RIx ~~ covRI*RIy

wx1 ~~ wx1 # Variances
wy1 ~~ wy1 
wx2 ~~ wx2 # Residual variances
wy2 ~~ wy2 
wx3 ~~ wx3 
wy3 ~~ wy3 
wx4 ~~ wx4 
wy4 ~~ wy4 
wx5 ~~ wx5 
wy5 ~~ wy5 


x1 ~~ 0*x1 
x2 ~~ 0*x2 
x3 ~~ 0*x3 
x4 ~~ 0*x4 
x5 ~~ 0*x5 

y1 ~~ 0*y1 
y2 ~~ 0*y2 
y3 ~~ 0*y3 
y4 ~~ 0*y4 
y5 ~~ 0*y5 

wx1 ~~ 0*RIx
wx1 ~~ 0*RIy
wy1 ~~ 0*RIx
wy1 ~~ 0*RIy

'


# DPM

DPM <- '
Zx =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
Zy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5

Zx~~x1
Zx~~y1
Zy~~x1
Zy~~y1
x1~~y1


Zx ~~ varRIx*Zx
Zy ~~ varRIy*Zy 
Zx ~~ covRI*Zy



x2~~ur*y2
x3~~ur*y3
x4~~ur*y4
x5~~ur*y5


x2~ax*x1
x3~ax*x2
x4~ax*x3
x5~ax*x4

y2~ay*y1
y3~ay*y2
y4~ay*y3
y5~ay*y4


x2~by*y1
x3~by*y2
x4~by*y3
x5~by*y4

y2~bx*x1
y3~bx*x2
y4~bx*x3
y5~bx*x4

'





# Estimate models
CLPMout <- sem(model = CLPM, data = df)
RICLPMout <- sem(model = RICLPM,   data  = df)
DPMout <- sem(model = DPM,   data  = df)

CLPMparam <- summary(CLPMout,standardized = TRUE, fit.measures=TRUE)
RICLPMparam <- summary(RICLPMout,standardized = TRUE, fit.measures=TRUE)
DPMparam <- summary(DPMout,standardized = TRUE, fit.measures=TRUE)


# Create tables to summarize the results
## Fit table
fintablefit <- cbind(Model = c("CLPM", "RICLPM", "DPM"),as.data.frame(rbind(round(CLPMparam$fit[c(15, 1, 3, 4, 5, 9, 10, 17, 21)], digits = 4),
                    round(RICLPMparam$fit[c(15, 1, 3, 4, 5, 9, 10, 17, 21)], digits = 4),
                    round(DPMparam$fit[c(15, 1, 3, 4, 5, 9, 10, 17, 21)], digits = 4))))



## Parameter estimates table
fintable <- cbind(Model = c(rep("CLPM", 4), rep("RICLPM", 7), rep("DPM", 7)), rbind(
  CLPMparam$pe[CLPMparam$pe$label == "ax", ][1, ],
  CLPMparam$pe[CLPMparam$pe$label == "ay", ][1, ],
  CLPMparam$pe[CLPMparam$pe$label == "bx", ][1, ],
  CLPMparam$pe[CLPMparam$pe$label == "by", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "ax", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "ay", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "bx", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "by", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "varRIx", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "varRIy", ][1, ],
RICLPMparam$pe[RICLPMparam$pe$label == "covRI", ][1, ],
DPMparam$pe[DPMparam$pe$label == "ax", ][1, ],
DPMparam$pe[DPMparam$pe$label == "ay", ][1, ],
DPMparam$pe[DPMparam$pe$label == "bx", ][1, ],
DPMparam$pe[DPMparam$pe$label == "by", ][1, ],
DPMparam$pe[DPMparam$pe$label == "varRIx", ][1, ],
DPMparam$pe[DPMparam$pe$label == "varRIy", ][1, ],
DPMparam$pe[DPMparam$pe$label == "covRI", ][1, ]))

quad <- c(alpha, alpha, beta_focal, beta_focal)
fintable$True <- c(quad, quad, 0, 0, NA, quad, 0, 0, NA)



# Summary tables
## Fit results
fintablefit

## Parameter estimate results
fintable


# If you are interested in the outputs from each model, use the following.
summary(CLPMout,standardized = TRUE, fit.measures=T)
summary(RICLPMout,standardized = TRUE, fit.measures=T)
summary(DPMout,standardized = TRUE, fit.measures=T)




