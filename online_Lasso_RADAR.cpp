// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
arma::mat lasso_RADAR(arma::vec mu, arma::vec theta, arma::vec beta_k, arma::mat X_new, arma::mat Y_new, int N_new,
    int maxit, double tol, double eta, double lambda_s, double R_k){
    
    int d = mu.n_elem;
    double p = 2 * log(d) / (2 * log(d) - 1);
    double q = 2 * log(d);

    arma::vec beta_new;

    mu = mu + (X_new.t() * (X_new * theta - Y_new) / N_new)+ (lambda_s * sign(theta));

    double mu_q_norm = pow(sum(pow(abs(mu), q)), 1 / q);
    double xi_2 = (eta * mu_q_norm * R_k * (p - 1))  - 1;
    double xi = 0;
    if (xi_2 > 0){
          xi = xi_2;}
    
    theta = beta_k + ((pow(mu_q_norm, 2 - q) * R_k * R_k * eta * (p - 1) * pow(abs(mu), q - 1) % sign(mu)) /  (xi + 1) );

    mat A(d, 2, fill::zeros);

    A.col(0) = mu;
    A.col(1) = theta;

    return A;
}


//[[Rcpp::export]]
List online_Lasso_full(arma::mat X, arma::vec y, arma::mat beta_lambda, arma::vec subset_index, int N_new,
    arma::uvec index1, arma::uvec index2, arma::mat gamma_new,
	arma::vec zz_r, arma::vec ztx_r, arma::vec zty_r, arma::mat ztX_r,
	arma::vec lambda_seq, double eta, int b, int maxit, double tol, arma::vec beta_tilde, arma::mat gamma_tilde,
	double R_k, arma::mat mu, arma::mat theta, arma::vec tilde_theta, int count_int_k,
	arma::mat mu_gamma, arma::mat theta_gamma, arma::mat tilde_theta_gamma, int k){
    
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
        arma::vec beta_i = beta_lambda.col(i);
        arma::vec mu_i = mu.col(i);
        arma::vec theta_i = theta.col(i);
        arma::mat A = lasso_RADAR(mu_i, theta_i, beta_i, X, y, N_new, maxit, tol, eta, lambda_i, R_k);

        mu.col(i) = A.col(0);
        theta.col(i) = A.col(1);

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
            arma::vec mu_i = mu.col(i);
            arma::vec theta_i = theta.col(i);
            arma::mat A = lasso_RADAR(mu_i, theta_i, beta_i, X1, y1, N_new, maxit, tol, eta, lambda_i, R_k);

            mu.col(i) = A.col(0);
            theta.col(i) = A.col(1);

            double PE = as_scalar((y2 - X2 * theta.col(i)).t() * 
               (y2 - X2 * theta.col(i))) / (n - N1);
            pred_error(i) = PE;
        }

    } else{
        for (int i = 0; i < s; i++){
            vec beta_i = beta_lambda.col(i);
            double PE = as_scalar((y - X * theta.col(i)).t() * (y - X * theta.col(i))) / n;
            pred_error(i) = PE;
        }
    }
    
    uword min_index = pred_error.index_min();
    
    lambda_s = lambda_seq(min_index);

    // Lasso estimator
    arma::vec theta_new = theta.col(min_index);

    tilde_theta = (theta_new / count_int_k) + (count_int_k - 1) * (tilde_theta / count_int_k);
    // debias on every predictor

    arma::mat used_tilde_theta_gamma = zeros<mat>(p-1, sub_length);
    
    for (arma::uword l = 0; l < sub_length; l++){
        
        arma::uword r = subset_index(l) - 1;

        arma::vec col_Range = arma::regspace<arma::vec>(0, p - 1);
        arma::mat X_r = X.cols(find(col_Range != r));
        arma::vec x_r = X.col(r);


        arma::mat A_gamma = lasso_RADAR(mu_gamma.col(l), theta_gamma.col(l), gamma_new.col(l), 
            X_r, x_r, N_new, maxit, tol, eta, lambda_s, R_k);

        mu_gamma.col(l) = A_gamma.col(0);
        theta_gamma.col(l) = A_gamma.col(1);
        
        tilde_theta_gamma.col(l) = (theta_gamma.col(l) / count_int_k) + (count_int_k-1) * (tilde_theta_gamma.col(l) / count_int_k);
        
        
        
        if (k < 3){used_tilde_theta_gamma.col(l) = zeros<vec>(p-1);
        }  // We will not use the estimation when data size is small.
        else{used_tilde_theta_gamma.col(l) = gamma_new.col(l);}
        

        vec z_r = x_r - X_r * used_tilde_theta_gamma.col(l);               // We may use the gamma_new rather than tilde_theta_gamma
        
        
        if (k < 3){vec z_r = x_r*0 - X_r * used_tilde_theta_gamma.col(l);}  
        
        zz_r(r) = zz_r(r) + as_scalar(z_r.t() * z_r);
        ztx_r(r) = ztx_r(r) + as_scalar(z_r.t() * x_r);
        zty_r(r) = zty_r(r) + as_scalar(z_r.t() * y);
        ztX_r.col(l) = ztX_r.col(l) + X.t() * z_r;

       
        beta_de(l) = beta_lambda.col(1)(r) + (zty_r(r) - as_scalar(ztX_r.col(l).t() * beta_lambda.col(1))) / ztx_r(r);  
        sd_de(l) = sqrt(zz_r(r)) / ztx_r(r);
    }
    sigma_ols = as_scalar(sqrt((y - X * beta_lambda.col(1)).t() * (y - X * beta_lambda.col(1)) / n));

    
	return List::create(
        Named("beta_de") = beta_de,
	 	Named("sd_de") = sd_de,
        Named("sigma_ols") = sigma_ols,
        Named("lambda_s") = lambda_s,
        Named("beta_lambda_new") = beta_lambda, //beta_lambda_new,
 
        Named("zz_r") = zz_r,
        Named("ztx_r") = ztx_r,
        Named("zty_r") = zty_r,
        Named("ztX_r") = ztX_r,
        Named("gamma_new") = gamma_new,

        Named("mu") = mu,
        Named("theta") = theta,
        Named("tilde_theta") = tilde_theta,
        Named("mu_gamma") = mu_gamma,
        Named("theta_gamma") = theta_gamma,
        Named("tilde_theta_gamma") = tilde_theta_gamma
        );
}



