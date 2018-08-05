
#' Update state sequence
#' 
#' @param obs Matrix of observations, with columns "x", "y", "time", and "state"
#' @param switch Matrix of state switches, with columns "time" and "state"
#' @param updateLim Vector of two elements: min and max length of updated interval
#' (each an integer, corresponding to a number of observations)
#' @param updateProbs Vector of probabilities of picking each number from
#' updateLim[1] to updateLim[2], for the length of the update interval
#' @param Q Square matrix of transition rates
#' 
#' @return List of two elements:
#' \item{newSwitch}{Updated matrix of transitions}
#' \item{newData}{Updated matrix of all data (observations and transitions)}
#' 
#' @details The update is done with the function \code{sample_path} from the
#' package ECctmc (Fintzi, 2017).
#' 
#' @references
#' Jon Fintzi (2017). ECctmc: Simulation from Endpoint-Conditioned Continuous 
#' Time Markov Chains. R package version 0.2.4. 
#' https://CRAN.R-project.org/package=ECctmc
#' 
#' @importFrom ECctmc sample_path
#' @export
updateState <- function(obs, switch, updateLim, updateProbs=NULL, Q)
{
    nbObs <- nrow(obs)
    
    if(is.null(updateProbs))
        updateProbs <- rep(1, length(updateLim[1]:updateLim[2]))/length(updateLim[1]:updateLim[2])
    
    # select interval to update
    if(updateLim[1]<updateLim[2]) {
        len <- sample(updateLim[1]:updateLim[2], size=1, prob=updateProbs)
    } else {
        len <- updateLim[1]
    }
    begin <- sample(1:(nbObs-len),size=1)
    end <- begin + len
    Tbeg <- obs[begin,"time"]
    Tend <- obs[end,"time"]
    
    # sample state sequence conditional on start and end state
    path <- sample_path(a=obs[begin,"state"], b=obs[end,"state"], 
                        t0=Tbeg, t1=Tend, Q=Q)
    path <- path[-c(1,nrow(path)),] # remove 1st and last rows (observations)
    
    # update state sequence
    newSwitch <- rbind(switch[switch[,"time"]<Tbeg,], 
                       path, 
                       switch[switch[,"time"]>Tend,])
    
    # remove switches into same state
    fakeSwitch <- which(newSwitch[-1,"state"]==newSwitch[-nrow(newSwitch),"state"]) + 1
    if(length(fakeSwitch)>0)
        newSwitch <- newSwitch[-fakeSwitch,]
    
    newData <- rbind(obs,cbind(NA,NA,newSwitch))
    newData <- newData[order(newData[,"time"]),]
    
    # update state sequence for new switches
    ind <- which(newData[,"time"]>Tbeg & newData[,"time"]<Tend)
    for(t in ind) {
        if(!is.na(newData[t,"x"])) {
            newData[t,"state"] <- newData[t-1,"state"]
        }
    }
    
    return(list(newSwitch=newSwitch, newData=newData, len=len))
}
