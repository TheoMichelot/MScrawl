
## This script runs a simulation study similar to Section 5.1 of the paper
## "State-switching continuous-time correlated random walks". This code uses the 
## R package MScrawl, which implements the MCMC algorithm presented in the paper.

# Install MScrawl from Github
devtools::install_github("TheoMichelot/MScrawl")
library(MScrawl)
library(ggplot2)

# Colour palette (for plots)
pal <- c("firebrick", "royalblue")

###################
## Simulate data ##
###################
set.seed(2732)

# Number of simulated points
nbSim <- 1e4

# Times of simulation
simTimes <- 0.1*(1:nbSim)

# Simulation parameters
truebeta <- c(1, 0.3)
truesigma <- c(1, 3)
Q <- matrix(c(-0.03, 0.03,
              0.03, -0.03),
            nrow = 2, ncol = 2)

# Simulate locations
sim <- simIOU(obsTimes = simTimes, beta = truebeta, 
              sigma = truesigma, Q = Q)

# Thin to emulate real data
keep <- sample(1:nbSim, size=nbSim/10, replace=FALSE)
sim <- sim[sort(keep),]
nbObs <- nrow(sim)

# Plot thinned trajectory
p0 <- ggplot(sim, aes(x, y, col = factor(state), group = NA)) + 
    geom_path(size = 0.5) + 
    geom_point(size=0.5) + 
    scale_color_manual(values = pal, name = "State") + 
    coord_equal()

####################
## Setup for MCMC ##
####################
nbStates <- 2

# Proposal parameters
props <- list(betaSD = 0.05,
              sigmaSD = 0.05,
              updateLim = c(3,70),
              updateProbs = rep(1/68,68))

# Initial parameters
Q0 <-  matrix(c(-0.1, 0.1, 0.1, -0.1), nrow = 2, ncol = 2)
state0 <- sample(1:nbStates, size = nbObs, replace = TRUE)
state0[1] <- sim$state[1]
state0[nbObs] <- sim$state[nbObs]
inits <- list(beta = c(0.7, 0.5), 
              sigma = c(1.5, 2.5),
              Q = Q0, 
              state = state0)

# Prior parameters (centred on true parameters for beta and sigma)
priors <- list(mean = log(c(truebeta, truesigma)),
               sd = rep(1, 2*nbStates),
               shape = 2,
               rate = 10,
               con = 0)

# Tuning parameters
tunes <- list(thinStates = 10)

# Parameters of the Kalman filter
kalmanpars <- list(Hmat = matrix(0, nbObs, 4),
                   a0 = c(sim$x[1], 0, sim$y[1], 0),
                   P0 = diag(c(0, 1, 0, 1)))

# Run 1e5 MCMC iterations
nbIter <- 1e5
m <- runMCMC(track=sim[,1:3], nbStates=nbStates, nbIter=nbIter, inits=inits, 
             priors=priors, props=props, tunes=tunes, kalmanpars=kalmanpars)

# Unpack posterior samples
allparam <- m$allparam
allrates <- m$allrates
allstates <- m$allstates

#################################
## Plot summaries of estimates ##
#################################
# Define posterior samples to keep (remove first half as burn-in)
burnin <- nbIter/2
ind <- seq(burnin+1, nbIter, by=10) # thin by 10 for visualisation
pardf <- as.data.frame(allparam[ind,])

# beta vs sigma (state 1)
p1 <- ggplot(pardf, aes(V1, V3)) + 
    geom_point(col=pal[1], size=0.1) +
    coord_equal() +
    xlab(expression(beta[1])) + 
    ylab(expression(sigma[1])) +
    theme(axis.title = element_text(size=12), 
          axis.text = element_text(size=12)) +
    geom_point(aes(x,y), data.frame(x = truebeta[1], y = truesigma[1]), size = 3)

# beta vs sigma (state 2)
p2 <- ggplot(pardf, aes(V2, V4)) + 
    geom_point(col=pal[2], size=0.1) +
    coord_equal() +
    xlab(expression(beta[2])) + 
    ylab(expression(sigma[2])) +
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 12)) +
    geom_point(aes(x,y), data.frame(x = truebeta[2], y = truesigma[2]), size = 3)

# Plot posterior state probabilities
sp <- colMeans(allstates[-(1:nrow(allstates)/2),]) - 1
spdf <- data.frame(t = sim$time, sp = sp)
p3 <- ggplot(spdf, aes(t,sp)) + xlab(expression(t[i])) + 
    ylab(expression(Pr(S[t[i]]==2))) +
    geom_path(size=0.3) + geom_point(size=0.7, col=pal[sim$state]) +
    geom_hline(yintercept=0.5, size=0.3, linetype=3)

# Mean posterior state sequence
states <- round(colMeans(allstates[-(1:nrow(allstates)/2),]))

# Count how many locations were classified with posterior prob < 0.9
spmat <- cbind(1-sp, sp)
count <- 0
for(i in 1:nrow(spmat))
    if(spmat[i, sim$state[i]]<0.9)
        count <- count + 1
