
#' Update transition rates
#' 
#' @param nbStates Number of states
#' @param data Matrix of all data (observations and transitions)
#' @param switch Matrix of transitions
#' @param priorShape Shape of prior gamma distribution for the transition rates
#' @param priorRate Rate of prior gamma distribution for the transition rates
#' @param priorCon Concentration of prior Dirichlet distribution for
#' transition probabilities
#' 
#' @return Updated transition matrix
#' 
#' @importFrom gtools rdirichlet
#' @export
updateQ <- function(nbStates, data, switch, priorShape, priorRate, priorCon)
{
    # time spent in each state
    timeInbStates <- sapply(1:nbStates, function(i) 
        sum(diff(switch[,"time"])[which(switch[,"state"]==i)],na.rm=TRUE))
    
    # sample rates out of each state
    shape <- priorShape + table(switch[,"state"])
    rate <- priorRate + timeInbStates
    r <- rgamma(nbStates,shape,rate)
    
    # sample transition probabilities
    allCounts <- table(data[-nrow(data),"state"],data[-1,"state"])
    nonDiagCounts <- matrix(t(allCounts)[!diag(nbStates)],nrow=nbStates,byrow=TRUE)
    trProbs <- t(apply(nonDiagCounts, 1, function(x)
        rdirichlet(n=1, alpha=x+priorCon)))
    
    # update generator matrix from rates and transition probs
    Q <- -diag(nbStates)
    Q[!Q] <- t(trProbs)
    Q <- t(Q*rep(r,each=nbStates))
    
    return(Q)
}
