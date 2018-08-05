
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Make T matrix
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
arma::mat makeT(double beta, double dt) {
    arma::mat T(4,4, fill::zeros);
    
    T(0,0) = 1;
    T(2,2) = 1;
    T(0,1) = (1-exp(-beta*dt))/beta;
    T(1,1) = exp(-beta*dt);
    T(2,3) = (1-exp(-beta*dt))/beta;
    T(3,3) = exp(-beta*dt);
    
    return T;
}

//' Make Q matrix
//' 
//' @param beta Parameter beta of OU velocity process
//' @param sigma Parameter sigma of OU velocity process
//' @param dt Length of time interval
arma::mat makeQ(double beta, double sigma, double dt) {
    arma::mat Q(4,4, fill::zeros);
    
    Q(0,0) = (sigma/beta)*(sigma/beta)*(dt - 2/beta*(1-exp(-beta*dt)) + 
        1/(2*beta)*(1-exp(-2*beta*dt)));
    Q(1,1) = sigma*sigma/(2*beta) * (1-exp(-2*beta*dt));
    Q(0,1) = sigma*sigma/(2*beta*beta) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt));
    Q(1,0) = Q(0,1);
    Q.submat(2,2,3,3) = Q.submat(0,0,1,1);
    
    return Q;
}
