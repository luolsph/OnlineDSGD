// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
arma::vec lasso_SGD(arma::vec beta_k, arma::mat X_new, arma::mat Y_new, int N_new,
    int maxit, double tol, double eta, double lambda_s){
    
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    int p = beta_k.n_elem;

    arma::vec beta_new;

    while (!stop_flag){
        niter += 1;
        // cout << niter << endl;
        vec g = X_new.t() * (X_new * beta_k - Y_new) / N_new;
        beta_k = beta_k - eta * g;

        vec s1 = beta_k - eta * lambda_s * ones<vec>(p);
        vec s2 = - beta_k - eta * lambda_s * ones<vec>(p);
    
        uvec c1 = (zeros<vec>(p) <= s1);
        uvec c2 = (zeros<vec>(p) <= s2);

        beta_k = c1 % s1 - c2 % s2; 
        
        double gnorm = sqrt(as_scalar(g.t() * g));

        if(gnorm < tol) {converged = TRUE; stop_flag = TRUE;}
        if(niter >= maxit) {stop_flag = TRUE;}
    }
    beta_new = beta_k;
    return beta_new;
}


//[[Rcpp::export]]
List online_Lasso_full_ASGD(arma::mat X, arma::vec y, arma::mat beta_lambda, arma::vec subset_index, int N_new,
    arma::uvec index1, arma::uvec index2, arma::mat gamma_new,
	arma::vec zz_r, arma::vec ztx_r, arma::vec zty_r, arma::mat ztX_r,
	arma::vec lambda_seq, double eta, int b, int maxit, double tol, arma::vec beta_tilde, arma::mat gamma_tilde){
    
    int p = beta_lambda.n_rows;
    int s = beta_lambda.n_cols;
    int n = y.n_elem;
    int sub_length = subset_index.n_elem;
    
    arma::vec pred_error = zeros<vec>(s);
    arma::mat beta_lambda_new = zeros<mat>(p, s);


    arma::vec beta_de = zeros<vec>(sub_length);
    arma::vec sd_de = zeros<vec>(sub_length);
    double sigma_ols;

    double lambda_s;


    // Calculate beta mat on every possible lambda
    for (int i = 0; i < s; i++){
        double lambda_i = lambda_seq(i); 
        vec beta_i = beta_lambda.col(i);
        vec beta_lambda_new_i = lasso_SGD(beta_i, X, y, N_new, maxit, tol, eta, lambda_i);
        beta_lambda_new.col(i) = beta_lambda_new_i;
    }

    if (b == 1){  // cross-validation within the first batch
        mat d1 = join_rows(X, y);
        mat train = d1.rows(index1);
        mat test = d1.rows(index2);
        
        mat X1 = train.cols(0, (p - 1));
        mat y1 = train.col(p);

        int N1 = y1.n_elem;

        mat X2 = test.cols(0, (p - 1));
        mat y2 = test.col(p);
        
        arma::mat beta_lambda_init = zeros<mat>(p, s);
        for (int i = 0; i < s; i++){
            vec beta_i = beta_lambda.col(i);
            double lambda_i = lambda_seq(i);
            vec beta_lambda_init_i = lasso_SGD(beta_i, X1, y1, N1, maxit, tol, eta, lambda_i);
            beta_lambda_init.col(i) = beta_lambda_init_i;
            double PE = as_scalar((y2 - X2 * beta_lambda_init_i).t() * 
               (y2 - X2 * beta_lambda_init_i)) / (n - N1);
            pred_error(i) = PE;
        }

    } else{
        for (int i = 0; i < s; i++){
            vec beta_i = beta_lambda.col(i);
            double PE = as_scalar((y - X * beta_i).t() * (y - X * beta_i)) / n;
            pred_error(i) = PE;
        }
    }
    
    uword min_index = pred_error.index_min();
    
    lambda_s = lambda_seq(min_index);
    
    // Lasso estimator
    arma::vec beta_new = beta_lambda_new.col(min_index);
    if (b == 1){beta_tilde = beta_new;} 
    else{beta_tilde = (beta_new / b) + (b - 1) * (beta_tilde / b);}
    
    // debias on selected predictor
    for (arma::uword l = 0; l < sub_length; l++){

        arma::uword r = subset_index(l) - 1;
   
        arma::vec col_Range = arma::regspace<arma::vec>(0, p - 1);
        arma::mat X_r = X.cols(find(col_Range != r));
        arma::vec x_r = X.col(r);
   
        
        gamma_new.col(l) = lasso_SGD(gamma_new.col(l), X_r, x_r, N_new, 
            maxit, tol, eta, lambda_s);
        
        if (b == 1){gamma_tilde.col(l) = gamma_new.col(l);} 
        else{gamma_tilde.col(l) = (gamma_new.col(l) / b) + (b - 1) * (gamma_tilde.col(l) / b);}
        
        vec z_r = x_r - X_r * gamma_tilde.col(l);
        zz_r(r) = zz_r(r) + as_scalar(z_r.t() * z_r);
        ztx_r(r) = ztx_r(r) + as_scalar(z_r.t() * x_r);
        zty_r(r) = zty_r(r) + as_scalar(z_r.t() * y);
        ztX_r.col(l) = ztX_r.col(l) + X.t() * z_r;

        beta_de(l) = beta_tilde(r) + (zty_r(r) - as_scalar(ztX_r.col(l).t() * beta_tilde)) / ztx_r(r);

        sd_de(l) = sqrt(zz_r(r)) / ztx_r(r);
    }
    sigma_ols = as_scalar(sqrt((y - X * beta_tilde).t() * (y - X * beta_tilde) / n));
  

    
	return List::create(
        Named("beta_de") = beta_de,
	 	Named("sd_de") = sd_de,
        Named("sigma_ols") = sigma_ols,
        Named("lambda_s") = lambda_s,
        Named("beta_lambda_new") = beta_lambda_new,
        Named("beta_new") = beta_new,
        Named("zz_r") = zz_r,
        Named("ztx_r") = ztx_r,
        Named("zty_r") = zty_r,
        Named("ztX_r") = ztX_r,
        Named("gamma_new") = gamma_new,
        Named("beta_tilde") = beta_tilde,
        Named("gamma_tilde") = gamma_tilde
        );
}



