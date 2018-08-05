
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "mat.hpp"

//' Kalman filter
//' 
//' This code is adapted from the package crawl (Johnson et al., 2008).
//' 
//' @param data Matrix of data, including (in that order) columns "x", "y", 
//' "time", and "state".
//' @param param Vector of movement parameters (beta and sigma)
//' @param Hmat Matrix of observation error variance (four columns, and one row 
//' for each row of data)
//' @param a0 Initial state estimate vector
//' @param P0 Initial estimate covariance matrix
//' 
//' @return Log-likelihood
//' 
//' @references
//' Johnson, D. S., London, J. M., Lea, M. A., & Durban, J. W. (2008).
//' Continuous-time correlated random walk model for animal telemetry data
//' Ecology, 89(5), 1208-1215.
//' 
//' @export
// [[Rcpp::export]]
double kalman_rcpp(arma::mat& data, arma::vec param, arma::mat& Hmat, 
                   arma::vec a0, arma::mat P0)
{
    int nstate = param.size()/2;
    int nbData = data.n_rows;
    
    // unpack data
    arma::mat X = data.cols(0,1); // (x,y)
    arma::vec time = data.col(2); // time
    arma::vec S = data.col(3); // state
    // time intervals
    arma::vec dt(nbData, fill::ones);
    dt.subvec(0,nbData-2) = diff(time);
    
    // unpack parameters
    arma::vec beta = param.subvec(0,nstate-1);
    arma::vec sigma = param.subvec(nstate,2*nstate-1);
    
    // define all matrices and vectors needed below
    arma::mat Z(2, 4, fill::zeros);
    Z(0,0) = 1;
    Z(1,2) = 1;
    arma::mat H(2,2, fill::zeros);
    arma::mat T(4,4);
    arma::mat Q(4,4);
    arma::mat F(2,2, fill::zeros);
    arma::mat K(4,2, fill::zeros);
    arma::mat L(4,4, fill::zeros);
    arma::colvec u(2, fill::zeros);
    double detF;
    
    // initial state mean
    arma::vec aest(4, fill::zeros); 
    aest = a0;
    // initial state covariance matrix
    arma::mat Pest(4,4, fill::zeros);
    Pest = P0;
    
    // Kalman filter iterations
    double llk = 0;
    for(int i=0; i<nbData; i++) {
        arma::mat T = makeT(beta(S(i)-1),dt(i));
        arma::mat Q = makeQ(beta(S(i)-1),sigma(S(i)-1),dt(i));

        if(R_IsNA(X(i,0))) {
            // if missing observation
            aest = T*aest;
            Pest = T*Pest*T.t() + Q;
        } else {
            H(0,0) = Hmat(i,0);
            H(1,1) = Hmat(i,1);
            H(0,1) = Hmat(i,2);
            H(1,0) = Hmat(i,3);

            // measurement residual
            u = X.row(i).t() - Z*aest;
            // residual covariance
            F = Z*Pest*Z.t() + H; 
            detF = F(0,0)*F(1,1) - F(1,0)*F(0,1);

            if(detF<=0) {
                aest = T*aest;
                Pest = T*Pest*T.t() + Q;
            } else {
                // update log-likelihood
                llk = llk - (log(detF) + dot(u,solve(F,u)))/2;
                // Kalman gain
                K = T*Pest*Z.t()*F.i();
                // update state estimate
                aest = T*aest + K*u;
                // update estimate covariance
                L = T - K*Z;
                Pest = T*Pest*L.t() + Q;
            }
        }
    }
    
    return llk;
}
