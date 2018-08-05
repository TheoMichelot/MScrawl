
#' Make covariance matrix
#' 
#' @param beta Parameter beta of the OU process
#' @param sigma Parameter sigma of the OU process
#' @param dt Time interval
makeSigma <- function(beta,sigma,dt) {
    Q <- matrix(0,4,4)
    Q[1,1] <- sigma^2/(2*beta) * (1-exp(-2*beta*dt))
    Q[2,2] <- (sigma/beta)^2 * (dt + (1-exp(-2*beta*dt))/(2*beta) - 
                                    2*(1-exp(-beta*dt))/beta)
    Q[1,2] <- sigma^2/(2*beta^2) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt))
    Q[2,1] <- Q[1,2]
    Q[3:4,3:4] <- Q[1:2,1:2]
    return(Q)
}

#' Simulate from 2D IOU process
#' 
#' @param obsTimes Vector of times of observations
#' @param beta Parameter beta of OU velocity process
#' @param sigma Parameter sigma of OU velocity process
#' @param Q Infinitesimal generator
#' 
#' @return Simulated track
#' 
#' @export
simIOU <- function(obsTimes, beta, sigma, Q)
{
    nbStates <- length(beta)
    nbObs <- length(obsTimes)
    obsTimes <- obsTimes - obsTimes[1]
    
    ############################
    ## Simulate state process ##
    ############################
    S <- NULL
    switchTimes <- NULL
    st <- 1 # start in state 1
    t <- rexp(1,-Q[st,st])
    while(t<obsTimes[nbObs]) {
        if(nbStates>2) {
            probs <- Q[st,-st]/sum(Q[st,-st])
            st <- sample((1:nbStates)[-st],size=1,prob=probs)            
        } else {
            st <- (1:2)[-st]
        }
        
        S <- c(S,st)
        switchTimes <- c(switchTimes,t)
        t <- t + rexp(1,-Q[st,st])
    }
    
    times <- c(obsTimes,switchTimes)
    nbData <- length(times)
    
    S <- c(rep(NA,nbObs),S)
    S <- S[order(times)]
    S[1] <- 1
    
    for(t in 1:nbData)
        if(is.na(S[t]))
            S[t] <- S[t-1]
    
    times <- sort(times)
    isObs <- (times %in% obsTimes)
    dt <- diff(times)
    nbData <- length(times)
    
    ##############################################
    ## Simulate velocity and position processes ##
    ##############################################
    data <- matrix(0,nrow=nbData,ncol=4)
    for(t in 2:nbData) {
        s <- S[t]
        mu <- c(exp(-beta[s]*dt[t-1])*data[t-1,1], 
                data[t-1,2]+data[t-1,1]/beta[s]*(1-exp(-beta[s]*dt[t-1])),
                exp(-beta[s]*dt[t-1])*data[t-1,3], 
                data[t-1,4]+data[t-1,3]/beta[s]*(1-exp(-beta[s]*dt[t-1])))
        Sigma <- makeSigma(beta[s],sigma[s],dt[t-1])
        data[t,] <- mvrnorm(1,mu,Sigma)
    }
    
    return(data.frame(x=data[,2], y=data[,4], time=times, state=S, 
                      vx=data[,1], vy=data[,3])[isObs,])
}
