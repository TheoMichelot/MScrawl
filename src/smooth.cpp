
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "mat.hpp"

//' Kalman filter and smoother
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
//' @return List of: log-likelihood, state estimates, and state estimate covariance
//' 
//' @references
//' Johnson, D. S., London, J. M., Lea, M. A., & Durban, J. W. (2008).
//' Continuous-time correlated random walk model for animal telemetry data
//' Ecology, 89(5), 1208-1215.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List smooth_rcpp(arma::mat& data, arma::vec param, arma::mat& Hmat, 
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
    
    arma::mat pred(4,nbData, fill::zeros);
    arma::cube predVar(4,4,nbData, fill::zeros);
    
    arma::mat Z(2,4, fill::zeros); 
    Z(0,0) = 1;
    Z(1,2) = 1;
    
    arma::mat T(4,4, fill::zeros); 
    arma::mat Q(4,4, fill::zeros);
    
    arma::cube F(2,2,nbData, fill::zeros);
    arma::mat H(2,2, fill::zeros);
    arma::cube K(4,2,nbData, fill::zeros);
    arma::cube L(4,4,nbData, fill::zeros);
    arma::mat u(2,nbData, fill::zeros);
    
    // initial state mean
    arma::mat aest(4, nbData+1, fill::zeros);
    aest.col(0) = a0;

    // initial state covariance matrix
    arma::cube Pest(4,4,nbData+1, fill::zeros);
    Pest.slice(0) = P0;
    
    arma::colvec r(4, fill::zeros);
    arma::mat N(4,4, fill::zeros);
    
    double llk = 0;
    
    arma::mat res(nbData, 2, fill::zeros);
    
    // Forward filter
    for(int i=0; i<nbData; i++) {
        T = makeT(beta(S(i)-1),dt(i));
        Q = makeQ(beta(S(i)-1),sigma(S(i)-1),dt(i));
        
        if(R_IsNA(X(i,0))) {
            // missing observation
            aest.col(i+1) = T*aest.col(i);
            Pest.slice(i+1) = T*Pest.slice(i)*T.t() + Q;
            L.slice(i) = T;
        } else {
            H(0,0) = Hmat(i,0);
            H(1,1) = Hmat(i,1); 
            H(0,1) = Hmat(i,2); 
            H(1,0) = Hmat(i,2);
            
            u.col(i) = X.row(i).t() - Z*aest.col(i);
            F.slice(i) = Z*Pest.slice(i)*Z.t() + H;
            
            // standardized residuals (if F diagonal)
            res.row(i) = u.col(i).t()/sqrt(F.slice(i).diag().t());
            
            llk = llk - (log(det(F.slice(i))) + dot(u.col(i),solve(F.slice(i),u.col(i))))/2;
            
            K.slice(i) = T*Pest.slice(i)*Z.t()*F.slice(i).i();     
            L.slice(i) = T - K.slice(i)*Z;
            aest.col(i+1) = T*aest.col(i) + K.slice(i)*u.col(i);
            Pest.slice(i+1) = T*Pest.slice(i)*L.slice(i).t() + Q;
        }
    }
    
    // Backward smooth
    for(int j=nbData; j>0; j--) {
        if(R_IsNA(X(j-1,0)) || F.slice(j-1)(0,0)*F.slice(j-1)(1,1)==0) {
            r = L.slice(j-1).t() * r;
            N = L.slice(j-1).t() * N * L.slice(j-1);
        } else {
            r = Z.t()*solve(F.slice(j-1),u.col(j-1)) + L.slice(j-1).t() * r;
            N = Z.t() * solve(F.slice(j-1),Z) + L.slice(j-1).t()*N*L.slice(j-1); 
        }
        
        pred.col(j-1) = aest.col(j-1) + Pest.slice(j-1)*r;
        predVar.slice(j-1) = Pest.slice(j-1) - Pest.slice(j-1)*N*Pest.slice(j-1); 
    }
    
    return Rcpp::List::create(
        Rcpp::Named("llk") = llk, 
        Rcpp::Named("pred") = pred,
        Rcpp::Named("predVar") = predVar,
        Rcpp::Named("res") = res);
}
