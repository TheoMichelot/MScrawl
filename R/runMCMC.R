
#' Run MCMC iterations
#' 
#' @param track Dataframe of data, with columns "x", "y", and "time"
#' @param nbStates Number of states
#' @param nbIter Number of iterations
#' @param inits List of initial parameters 
#' (beta, sigma, Q, state)
#' @param priors List of parameters of prior distributions, with components:
#' \itemize{
#'   \item{"mean":} Vector of means for normal priors on movement parameters, of length 2*nbStates
#'   \item{"sd":} Vector of standard deviations for normal priors on movement parameters, of length
#' 2*nbStates
#'   \item{"shape":} Vector of shapes of gamma priors for the transition rates
#'   \item{"rate":} Vector of rates of gamma priors for the transition rates
#'   \item{"con":} Vector of concentrations of Dirichlet priors for transition probabilities
#' }
#' @param props List of parameters of proposal distributions, with components:
#' \itemize{
#'   \item{"betaSD":} Scalar standard deviation for normal proposal distribution of beta
#'   \item{"sigmaSD":} Scalar standard deviation for normal proposal distribution of sigma
#'   \item{"updateLim":} Vector of two values: min and max length of updated state sequence
#'   \item{"updateProbs":} Probability for each element of updateLim[1]:updateLim[2] (if NULL,
#'   all values are equiprobable)
#' }
#' @param tunes List of tuning parameters, with components:
#' \itemize{
#'   \item{"thinStates":} Thinning factor for the posterior state sequences (needed because
#' of memory limitations)
#' }
#' @param kalmanpars List of parameters of the Kalman filter, with components:
#' \itemize{
#'   \item{"Hmat":} Matrix of observation error variance (four columns, and one row 
#' for each row of data)
#'   \item{"a0":} Initial state estimate vector
#'   \item{"P0":} Initial estimate covariance matrix
#' }
#' @param updateState Logical. If FALSE, the state process is not updated
#' (for exploratory analysis only)
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats dnorm runif rnorm rexp rgamma
#' @export
#' 
#' @useDynLib MScrawl
runMCMC <- function(track, nbStates, nbIter, inits, priors, props, tunes, kalmanpars,
                    updateState=TRUE)
{
    nbObs <- nrow(track)
    
    ######################
    ## Unpack arguments ##
    ######################
    # initial parameters
    beta <- inits$beta
    sigma <- inits$sigma
    param <- c(beta,sigma)
    Q <- inits$Q
    state0 <- inits$state
    if(is.null(beta) | is.null(sigma) | is.null(Q) | is.null(state0))
        stop("'inits' should have components beta, sigma, Q, and state")
    
    # Kalman filter parameters
    Hmat <- kalmanpars$Hmat
    a0 <- kalmanpars$a0
    P0 <- kalmanpars$P0
    if(is.null(Hmat) | is.null(a0) | is.null(P0))
        stop("'kalmanpars' should have components Hmat, a0, and P0")
    
    # unpack prior parameters
    priorMean <- priors$mean
    priorSD <- priors$sd
    priorShape <- priors$shape
    priorRate <- priors$rate
    priorCon <- priors$con
    if(is.null(priorMean) | is.null(priorSD) | is.null(priorShape) | 
       is.null(priorRate) | is.null(priorCon)) {
        stop("'priors' should have components mean, sd, shape, rate, and con")
    }
    
    # unpack proposal parameters
    betaPropSD <- props$betaSD
    sigmaPropSD <- props$sigmaSD
    updateLim <- props$updateLim
    updateProbs <- props$updateProbs
    if(is.null(betaPropSD) | is.null(sigmaPropSD) | is.null(updateLim))
        stop("'props' should have components betaSD, sigmaSD, and updateLim")
    
    if(is.null(updateProbs))
        updateProbs <- rep(1, length(updateLim[1]:updateLim[2]))/length(updateLim[1]:updateLim[2])
    
    if(length(updateLim[1]:updateLim[2])!=length(updateProbs))
        stop("'updateProbs' has the wrong length")
    
    # unpack tuning parameters
    thinStates <- tunes$thinStates
    
    ####################
    ## Initialisation ##
    ####################
    # Prepare data structures
    obs <- matrix(c(track[,"x"],track[,"y"],track[,"time"],state0),ncol=4)
    colnames(obs) <- c("x","y","time","state")
    
    indSwitch <- which(obs[-1,"state"]!=obs[-nbObs,"state"])+1
    switch <- matrix(c(obs[indSwitch,"time"]-0.001,rle(obs[,"state"])$values[-1]),ncol=2)
    colnames(switch) <- c("time","state")
    
    data <- rbind(obs,cbind(NA,NA,switch))
    data <- data[order(data[,"time"]),]
    
    # initialise Hmat (rows of 0s for transitions)
    HmatAll <- matrix(0, nrow(data), 4)
    HmatAll[which(!is.na(data[,"x"])),] <- Hmat
    
    # initial likelihood
    oldllk <- kalman_rcpp(data=data, param=param, Hmat=HmatAll, a0=a0, P0=P0)
    # initial log-prior
    oldlogprior <- sum(dnorm(log(param),priorMean,priorSD,log=TRUE))
    
    ###############################
    ## Loop over MCMC iterations ##
    ###############################
    allparam <- matrix(NA,nrow=nbIter,ncol=2*nbStates)
    allrates <- matrix(NA,nrow=nbIter,ncol=nbStates*(nbStates-1))
    allstates <- matrix(NA,nrow=nbIter/thinStates,ncol=nbObs) # uses a lot of memory!
    accSwitch <- rep(0,nbIter)
    accParam <- rep(0,nbIter)
    allLen <- rep(NA,nbIter)
    
    t0 <- Sys.time()
    for(iter in 1:nbIter) {
        if(iter%%100==0)
            cat("\rIteration ", iter, "/", nbIter, "... ", round(Sys.time()-t0,2),
                " -- accSwitch = ", round(sum(accSwitch)/iter*100), "%",
                " -- accParam = ", round(sum(accParam)/iter*100),"%", sep="")
        
        ######################################
        ## 1. Update discrete state process ##
        ######################################
        if(updateState) {
            upState <- updateState(obs=obs, switch=switch, updateLim=updateLim, 
                                   updateProbs=updateProbs, Q=Q)
            newData <- upState$newData
            newSwitch <- upState$newSwitch
            allLen[iter] <- upState$len
            
            # update Hmat (rows of 0s for transitions)
            newHmatAll <- matrix(0, nrow(newData), 4)
            newHmatAll[which(!is.na(newData[,"x"])),] <- Hmat 
            
            # Calculate acceptance ratio
            newllk <- kalman_rcpp(data=newData, param=param, Hmat=newHmatAll, a0=a0, P0=P0)
            logHR <- newllk - oldllk
            
            if(log(runif(1))<logHR) {
                # Accept new state sequence
                accSwitch[iter] <- 1
                switch <- newSwitch
                data <- newData
                obs <- data[!is.na(data[,"x"]),]
                oldllk <- newllk
                HmatAll <- newHmatAll
            }
        }
        
        ###################################
        ## 2. Update movement parameters ##
        ###################################
        # On working scale
        betaprimeW <- rnorm(nbStates,log(beta),betaPropSD)
        sigmaprimeW <- rnorm(nbStates,log(sigma),sigmaPropSD)
        newlogprior <- sum(dnorm(c(betaprimeW,sigmaprimeW),priorMean,priorSD,log=TRUE))
        
        # On natural scale
        betaprime <- exp(betaprimeW)
        sigmaprime <- exp(sigmaprimeW)
        
        # Calculate acceptance ratio
        newllk <- kalman_rcpp(data=data, param=c(betaprime,sigmaprime), Hmat=HmatAll, a0=a0, P0=P0)
        logHR <- newllk + newlogprior - oldllk - oldlogprior

        if(log(runif(1))<logHR) {
            # Accept new parameter values
            accParam[iter] <- 1
            beta <- betaprime
            sigma <- sigmaprime
            param <- c(beta,sigma)
            oldllk <- newllk
            oldlogprior <- newlogprior
        }
        
        ###############################
        ## 3. Update switching rates ##
        ###############################
        Q <- updateQ(nbStates=nbStates, data=data, switch=switch, 
                     priorShape=priorShape, priorRate=priorRate,
                     priorCon=priorCon)
        
        #########################
        ## Save posterior draw ##
        #########################
        allparam[iter,] <- param
        allrates[iter,] <- Q[!diag(nbStates)]    
        if(iter%%thinStates==0)
            allstates[iter/thinStates,] <- obs[,"state"]
    }
    cat("\n")
    
    return(list(allparam = allparam,
                allrates = allrates,
                allstates = allstates,
                accSwitch = accSwitch,
                allLen = allLen))
}
