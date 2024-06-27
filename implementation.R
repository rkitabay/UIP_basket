

library(R2jags)
library(entropy)


##-----------------Trial Data--------------------##

Nbask <- 6 # Number of cancer types
pi.null <- 0.15 # Null response rate
pi.alt <- 0.35 # Alternative response rate

# Trial data containing number of patients (n) and responders (x) for each cancer type in the vemurafenib trial
data <- data.frame(i = 1:Nbask, x = c(8,0,1,1,6,2), n = c(19,10,26,8,14,7))






##-----------------Parameter Specification--------------------##
#---------Common to 3 UPSiDe Methods-----------#

M <- 84 # Total amount parameter (introduced in Section 2.1 in the main manuscript)

pi.min <- 0.05
# pi_min is the minimum value for response rates when estimating Unit Information (UI) for each cancer type.
# If the observed response rate is 0, the UI approaches infinity.
# We recommend setting the value to 0.05 because this is a common value that represents the null response rate for cancers with no standard of care (i.e., the potential minimum in oncology).

alp.min <- 0.5
bet.min <- 0.5
# alp_min and bet_min are the minimum values for the parameters of the beta prior for each cancer type.
# Alpha and beta in a beta distribution must be greater than 0; however, very small values (e.g., 0.001) can cause JAGS to produce infinite density errors in this model.
# We recommend setting both to 0.5 to represent weak prior information and to avoid errors in various trial settings.


#---------For UPSiDe-D-----------#

z <- 1 # Parameter of the Dirichlet prior

#---------For UPSiDe-JS2-----------#

s <- 0.01 # Parameter of the Gamma prior






##-----------------Analysis--------------------##
#---------For UPSiDe-D-----------#

Ncomb <- Nbask * (Nbask - 1) / 2 # Number of combinations between two cancer types
ref <- ref_matrix(Nbask) # Matrix to refer the order of each element when transposing the matrix to vector
jlist <- j_matrix(Nbask) # Matrix to refer the cancer type number other than i (i is a row number)

# Observed response rate with adjustments for extreme values to prevent the division by zero
Means <- data$x / data$n
Means[Means < 0.001] <- 0.001
Means[Means > 0.999] <- 0.999

# Prepare data for JAGS
jagsData <- list(x = data$x, n = data$n, I = Nbask, jlist = jlist, zvct = rep(z, Ncomb),
                 ref = ref, upM = M, Means = Means, pi.min = pi.min, alp.min = alp.min, bet.min = bet.min)

# Run JAGS model ("UPSiDe-D.txt" contains the model specification)
fit <- jags(model.file = "UPSiDe-D.txt", data = jagsData, parameters.to.save = c("postM", "pi", "wmtx", "Mwmtx"),
            n.chains = 4, n.iter = 50000, n.thin = 3, n.burnin = 5000, quiet = T, DIC = F)


# Summary of the posterior distribution for pi and M
sum.fit <- fit[["BUGSoutput"]][["summary"]]
summary.pi <- sum.fit[Nbask * (Nbask - 1) + 1:(Nbask + 1), ]
summary.pi

# Posterior probabilities that response rates exceed the null response rate
pi.success <- pi_success(fit = fit[["BUGSoutput"]][["sims.list"]]$pi, pn = pi.null)
pi.success

# Posterior mean of Mw (i.e., sample size of information borrowed between two cancer types)
Mwlist <- matrix(sum.fit[1:(Nbask * (Nbask - 1)), "mean"], nrow = Nbask)
Mw <- matrix(0, nrow = Nbask, ncol = Nbask)
for(i in 1:Nbask){
  for(j in 1:Nbask){
    if(i > j) Mw[i, j] <- Mwlist[i, j]
    if(i < j) Mw[i, j] <- Mwlist[i, j - 1]
  }
}
Mw

# Posterior mean of w (i.e., weight parameters between two cancer types)
wlist <- matrix(sum.fit[Nbask * (Nbask - 1) + Nbask + 1 + 1:(Nbask * (Nbask - 1)), "mean"], nrow = Nbask)
w <- matrix(0, nrow = Nbask, ncol = Nbask)
for(i in 1:Nbask){
  for(j in 1:Nbask){
    if(i > j) w[i, j] <- wlist[i, j]
    if(i < j) w[i, j] <- wlist[i, j - 1]
  }
}
w





#---------For UPSiDe-JS1-----------#

Ncomb <- Nbask * (Nbask - 1) / 2 # Number of combinations between two cancer types
ref <- ref_matrix(Nbask) # Matrix to refer the order of each element when transposing the matrix to vector
jlist <- j_matrix(Nbask) # Matrix to refer the cancer type number other than i (i is a row number)

# Observed response rate with adjustments for extreme values to prevent the division by zero
Means <- data$x / data$n
Means[Means < 0.001] <- 0.001
Means[Means > 0.999] <- 0.999

# Calculate JS distance between response rates of two cancer types
distance <- JS_distance(I = Nbask, n = data$n, x = data$x, q1 = pi.alt, q0 = pi.null, a = 1, b = 1)

# Calculate of matrix of w and Mw
expd.mtx <- exp(- distance)
diag(expd.mtx) <- 0
sum.expd <- sum(expd.mtx[upper.tri(expd.mtx)])
w <- expd.mtx / 2 / sum.expd
wmtx <- matrix(w[row(w) != col(w)], nrow = Nbask)
Mwmtx <- M * wmtx

# Estimate posterior response rates
avg <- c()
prec <- c()
alp <- c()
bet <- c()
MwUImtx <- matrix(0, nrow = Nbask, ncol = Nbask - 1)
for(i in 1:Nbask){
  for(j in 1:(Nbask - 1)){
    # Calculate total of Unit Information (UI) borrowed from cancer types other than i
    MwUImtx[i,j] <- Mwmtx[i,j] * min(1 / (pi.min * (1 - pi.min)), 1 / (Means[jlist[i,j]] * (1 - Means[jlist[i,j]])))
  }
  # Weighted average of observed response rates
  avg[i] <- sum(wmtx[i, 1:(Nbask - 1)] * Means[jlist[i, 1:(Nbask - 1)]]) / sum(wmtx[i, 1:(Nbask - 1)])
  # Precision of the informative prior
  prec[i] <- sum(MwUImtx[i, 1:(Nbask - 1)])
  # Parameters alpha and beta of the beta prior for each cancer type
  alp[i] <- max(avg[i] * (avg[i] * (1 - avg[i]) * prec[i] - 1), alp.min)
  bet[i] <- max((1 - avg[i]) * (avg[i] * (1 - avg[i]) * prec[i] - 1), bet.min)
}
# Posterior alpha and beta
post.alp <- alp + data$x
post.bet <- bet + (data$n - data$x)

# Summary of the posterior distribution for pi
mn <- c(post.alp / (post.alp + post.bet))
sd <- c(sqrt(post.alp * post.bet / ((post.alp + post.bet) ^ 2 * (post.alp + post.bet + 1))))
lc <- c(qbeta(0.025, shape1 = post.alp, shape2 = post.bet))
md <- c(qbeta(0.500, shape1 = post.alp, shape2 = post.bet))
uc <- c(qbeta(0.975, shape1 = post.alp, shape2 = post.bet))
summary.pi <- data.frame(mean = mn, sd = sd, Lower.95pCL = lc, median = md, Upper.95pCL = uc)
rownames(summary.pi) <- sapply(1:Nbask, function(x){paste(c("pi[", x, "]"), collapse = "")})
summary.pi

# Posterior probabilities that response rates exceed the null response rate
pi.success <- c(1 - pbeta(pi.null, shape1 = post.alp, shape2 = post.bet))
names(pi.success) <- sapply(1:Nbask, function(x){paste(c("pi[", x, "]"), collapse = "")})
pi.success

# Posterior mean of Mw (i.e., sample size of information borrowed between two cancer types)
Mw <- M * w
Mw

# Posterior mean of w (i.e., weight parameters between two cancer types)
w






#---------For UPSiDe-JS2-----------#

Ncomb <- Nbask * (Nbask - 1) / 2 # Number of combinations between two cancer types
ref <- ref_matrix(Nbask) # Matrix to refer the order of each element when transposing the matrix to vector
jlist <- j_matrix(Nbask) # Matrix to refer the cancer type number other than i (i is a row number)

# Observed response rate with adjustments for extreme values to prevent the division by zero
Means <- data$x / data$n
Means[Means < 0.001] <- 0.001
Means[Means > 0.999] <- 0.999

# Calculate JS distance between response rates of two cancer types
distance<-JS_distance(I = Nbask, n = data$n, x = data$x, q1 = pi.alt, q0 = pi.null, a = 1, b = 1)
dvct <- c(distance[upper.tri(distance)])


# Prepare data for JAGS
jagsData <- list(x = data$x, n = data$n, I = Nbask, jlist = jlist, pris = s, dvct = dvct, Ncomb = Ncomb,
                 ref = ref, upM = M, Means = Means, pi.min = pi.min, alp.min = alp.min, bet.min = bet.min)

# Run JAGS model ("UPSiDe-JS.txt" contains the model specification)
fit <- jags(model.file = "UPSiDe-JS.txt", data = jagsData, parameters.to.save = c("postM", "pi", "wmtx", "Mwmtx"),
            n.chains = 4, n.iter = 50000, n.thin = 3, n.burnin = 5000, quiet = T, DIC = F)

# Summary of the posterior distibution for pi and M
sum.fit <- fit[["BUGSoutput"]][["summary"]]
summary.pi <- sum.fit[Nbask * (Nbask - 1) + 1:(Nbask + 1), ]
summary.pi

# Posterior probabilities that response rates exceed the null response rate
pi.success <- pi_success(fit = fit[["BUGSoutput"]][["sims.list"]]$pi, pn = pi.null)
pi.success

# Posterior mean of Mw (i.e., sample size of information borrowed between two cancer types)
Mwlist <- matrix(sum.fit[1:(Nbask * (Nbask - 1)), "mean"], nrow = Nbask)
Mw <- matrix(0, nrow = Nbask, ncol = Nbask)
for(i in 1:Nbask){
  for(j in 1:Nbask){
    if(i > j) Mw[i, j] <- Mwlist[i, j]
    if(i < j) Mw[i, j] <- Mwlist[i, j - 1]
  }
}
Mw

# Posterior mean of w (i.e., weight parameters between two cancer types)
wlist <- matrix(sum.fit[Nbask * (Nbask - 1) + Nbask + 1 + 1:(Nbask * (Nbask - 1)), "mean"], nrow = Nbask)
w <- matrix(0, nrow = Nbask, ncol = Nbask)
for(i in 1:Nbask){
  for(j in 1:Nbask){
    if(i > j) w[i, j] <- wlist[i, j]
    if(i < j) w[i, j] <- wlist[i, j - 1]
  }
}
w



