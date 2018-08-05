
#' Plot posterior samples of movement parameters
#' 
#' @param allpar Matrix of posterior samples
#' @param ind Vector of indices that should be included (e.g. to remove burn-in)
#' @param log Logical. If TRUE, plots posterior samples on log scale
#' 
#' @importFrom graphics par plot
#' 
#' @export
plotPars <- function(allpar, ind=NULL, log=FALSE) 
{
    nbStates <- ncol(allpar)/2
    if(is.null(ind))
        ind <- 1:nrow(allpar)
    
    pal <- c("firebrick", "royalblue", "seagreen")
    par(mfrow=c(1,nbStates), mar=c(5,5,1,1))
    
    for(s in 1:nbStates) {
        if(log) {
            x <- log(allpar[ind,s])
            y <- log(allpar[ind,nbStates+s])
            xlab <- bquote(paste("log(",beta[.(s)],")"))
            ylab <- bquote(paste("log(",sigma[.(s)],")"))
        } else {
            x <- allpar[ind,s]
            y <- allpar[ind,nbStates+s]
            xlab <- bquote(beta[.(s)])
            ylab <- bquote(sigma[.(s)])
        }
        xlim <- range(x)
        ylim <- range(y)
        plot(x, y, pch=19, cex=0.2, col=pal[s], xlim=xlim, ylim=ylim, asp=1, 
             xlab=xlab, ylab=ylab)
    }
    
    par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
}
