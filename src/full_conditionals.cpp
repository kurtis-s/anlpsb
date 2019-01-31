#include <Rcpp.h>
#include <RcppTN.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppTN)]]

/* Likelihood ############################################################### */
// [[Rcpp::export]]
double fnLogLikNoGamma(const double & y_tkj, const double & s_j,
                     const double & r_tk, const double & mu_tj) {
    // Does not include the y_tkj factorial part or the gamma part s_j
    const double s_j_inv = 1/s_j;
    const double param_product = r_tk * mu_tj * s_j;

    // For numerical stability when param_product is infinite
    if(y_tkj != 0) {
        return -y_tkj * log(1/param_product + 1) - s_j_inv * log(1 + param_product);
    }
    else { // The above fails when y_tkj is 0 and param_product is infinite: 0*log(Inf) = Nan
        return -s_j_inv * log(1 + param_product);
    }
}

// [[Rcpp::export]]
double fnLogLikNoFac(const double & y_tkj, const double & s_j,
                     const double & r_tk, const double & mu_tj) {
    // Does not include the y_tkj factorial part
    const double s_j_inv = 1/s_j;

    return lgamma(y_tkj + s_j_inv) - lgamma(s_j_inv) +
        fnLogLikNoGamma(y_tkj, s_j, r_tk, mu_tj);
}

// [[Rcpp::export]]
double fnOtuLogLikNoGamma(const double & s_j,
                        NumericVector r_vec, NumericVector mu_j_vec,
                        NumericVector rep_K, NumericVector y_j_vec) {
    // Log-likelihood for the jth OTU over time and replicates
    const int & n = rep_K.size();

    double total = 0;
    int tk=0;
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            total += fnLogLikNoGamma(y_j_vec(tk), s_j, r_vec(tk), mu_j_vec(i));
            tk++;
        }
    }

    return total;
}

// [[Rcpp::export]]
double fnOtuLogLikNoFac(const double & s_j,
                        NumericVector r_vec, NumericVector mu_j_vec,
                        NumericVector rep_K, NumericVector y_j_vec) {
    // Log-likelihood for the jth OTU over time and replicates
    const int & n = rep_K.size();

    double total = 0;
    int tk=0;
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            total += fnLogLikNoFac(y_j_vec(tk), s_j, r_vec(tk), mu_j_vec(i));
            tk++;
        }
    }

    return total;
}

// [[Rcpp::export]]
double fnTimepointLogLikNoGamma(NumericVector mu_t_vec, NumericVector s_vec,
                             const double & r_tk, NumericVector y_tk_vec) {
    // Log-likelihood for sample tk over OTUs
    const int & J = mu_t_vec.size();

    double log_lik_total = 0;
    for(int j=0; j<J; j++) {
        log_lik_total += fnLogLikNoGamma(y_tk_vec(j), s_vec(j), r_tk, mu_t_vec(j));
    }

    return log_lik_total;
}

/* Functions ################################################################ */
// [[Rcpp::export]]
bool acceptProposal(const double & curr_log_lik, const double & prop_log_lik, const double log_curr_to_prop_prob = 0, const double log_prop_to_curr_prob = 0) {
    double u = runif(1)(0);
    bool accept = FALSE;
    if(log(u) <= (prop_log_lik - curr_log_lik + log_prop_to_curr_prob - log_curr_to_prop_prob)) {
        accept = TRUE;
    }

    return accept;
}

/* Overdispersion ########################################################### */
// ### s_j
// [[Rcpp::export]]
double fnSjTildeLogLik(const double & s_j,
                       NumericVector r_vec, NumericVector mu_j_vec,
                       NumericVector rep_K, NumericVector y_j_vec,
                       const double & h_scal, const double & var_sig_2) {
    return R::dnorm(log(s_j), h_scal, sqrt(var_sig_2), TRUE) +
        fnOtuLogLikNoFac(s_j, r_vec, mu_j_vec, rep_K, y_j_vec);
}

// [[Rcpp::export]]
void fnSampSVec(NumericVector s_vec, NumericVector s_tilde_vec,
                  NumericVector r_vec, NumericMatrix Mu_mat,
                  NumericVector rep_K, NumericMatrix Y_mat,
                  const double & h_scal, const double & var_sig_2,
                  const double & s_j_tilde_proposal_sd) {
    const int J = s_vec.size();

    double curr_log_lik = 0;
    double prop_log_lik = 0;
    double s_j_prop = 0;
    for(int j=0; j<J; j++) {
        s_j_prop = exp(s_tilde_vec(j) + R::rnorm(0, s_j_tilde_proposal_sd));

        curr_log_lik = fnSjTildeLogLik(s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j), h_scal, var_sig_2);
        prop_log_lik = fnSjTildeLogLik(s_j_prop, r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j), h_scal, var_sig_2);

        if(acceptProposal(curr_log_lik, prop_log_lik)) {
            s_vec(j) = s_j_prop;
            s_tilde_vec(j) = log(s_j_prop);
        }
    }
}
// ###

/* Likelihood mean params ################################################### */
// ### w_ell
// [[Rcpp::export]]
double fnWlLogLik(const int & ell_idx1, const double & w_ell,
                  NumericVector param_vec, NumericVector lambda_vec,
                  NumericVector c_vec, NumericVector eta_vec,
                  const double & upsilon, const double & a_w,
                  const double & b_w, const double & u_2) {
    // param_vec is either r_tilde_vec or alpha0_vec

    // Guard to make sure w_ell is in the support
    if(w_ell < 0 || w_ell > 1) {
        return R_NegInf;
    }

    double normal_mean = (upsilon - w_ell*eta_vec[ell_idx1-1])/(1 - w_ell);

    double w_part_total = 0; // The part of the full conditional from w
    double param_part_total = 0; // The part of the full conditional from r_tk or alpha0_j
    double lambda_cur;
    double c_cur;
    double param_cur;
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for w_\ell^r; or j for w_\ell^\theta
        lambda_cur = lambda_vec[loop_idx];
        c_cur = c_vec[loop_idx];
        param_cur = param_vec[loop_idx];

        if(c_cur == ell_idx1) {
            if(lambda_cur == 0) {
                w_part_total += log(1 - w_ell);
                param_part_total += R::dnorm(param_cur, normal_mean, sqrt(u_2), TRUE);
            }
            else { // lambda = 1
                w_part_total += log(w_ell);
            }
        }
    }

    double prior_part = (a_w - 1)*log(w_ell) + (b_w - 1)*log(1 - w_ell); // From the prior

    return prior_part + w_part_total + param_part_total;
}

// [[Rcpp::export]]
void fnSampWVec(NumericVector w_vec, NumericVector param_vec,
                  NumericVector lambda_vec, NumericVector c_vec,
                  NumericVector eta_vec, const double & upsilon,
                  const double & a_w, const double & b_w, const double & u_2,
                  const double & w_ell_proposal_sd) {
    // param_vec is either r_tilde_vec or alpha0_vec
    const int L = w_vec.size();

    double curr_log_lik = 0;
    double prop_log_lik = 0;
    double w_ell_prop = -1;
    for(int ell=0; ell<L; ell++) {
        w_ell_prop = w_vec(ell) + R::rnorm(0, w_ell_proposal_sd);

        curr_log_lik = fnWlLogLik(ell + 1, w_vec(ell), param_vec, lambda_vec, c_vec, eta_vec, upsilon, a_w, b_w, u_2);
        prop_log_lik = fnWlLogLik(ell + 1, w_ell_prop, param_vec, lambda_vec, c_vec, eta_vec, upsilon, a_w, b_w, u_2);

        if(acceptProposal(curr_log_lik, prop_log_lik)) {
            w_vec(ell) = w_ell_prop;
        }
    }
}
// ###

// ### eta_ell
// [[Rcpp::export]]
void fnSampEtaVec(NumericVector eta_vec, NumericVector param_vec, NumericVector w_vec,
                  NumericVector lambda_vec, NumericVector c_vec,
                  const double & u_2, const double & upsilon, const double & b_eta_2) {
    // param_vec is either r_tilde_vec or alpha0_vec
    // eta_vec can be eta^r or eta^\theta
    const int L = eta_vec.size();
    const double & m3 = upsilon;
    const double rho3 = 1/b_eta_2;

    int ell_idx1;
    double norm_mean;
    double norm_precision;
    double norm_variance;
    double rho1;
    double rho2;
    double m1;
    double m2;
    double m1_total = 0;
    double m2_total = 0;
    double set_A_size = 0;
    double set_B_size = 0;
    for(int ell=0; ell<L; ell++) {
        ell_idx1 = ell + 1;
        m1_total = 0;
        m2_total = 0;
        set_A_size = 0;
        set_B_size = 0;
        for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for r.tk; or j for alpha0.j
            if(c_vec(loop_idx)==ell_idx1) {
                if(lambda_vec(loop_idx) == 1) { // Set 'A' in the full conditionals
                    m1_total += param_vec(loop_idx);
                    set_A_size++;
                }
                else { // lambda_vec[loop_idx]==0; set 'B' in full conditionals
                    m2_total += ( upsilon - (1-w_vec(ell)) * param_vec(loop_idx) )/w_vec(ell);
                    set_B_size++;
                }
            }
        }
        if(set_A_size > 0) {
            m1 = m1_total/set_A_size;
        }
        else {
            m1 = 0;
        }
        if(set_B_size > 0) {
            m2 = m2_total/set_B_size;
        }
        else {
            m2 = 0;
        }
        rho1 = (1/u_2) * set_A_size;
        rho2 = (1/u_2) * (w_vec(ell)/(1-w_vec(ell))) * set_B_size;
        norm_precision = rho1 + rho2 + rho3;
        norm_variance = 1/norm_precision;
        norm_mean = (rho1*m1 + rho2*m2 + rho3*m3)/norm_precision;

        eta_vec(ell) = R::rnorm(norm_mean, sqrt(norm_variance));
    }
}
// ###

// ### c_vec (c_tk or c_j)
double log_sum_exp(double u, double v) {
    // See http://andrewgelman.com/2016/06/11/log-sum-of-exponentials/
    return std::max(u, v) + log(exp(u - std::max(u, v)) + exp(v - std::max(u, v)));
}

// [[Rcpp::export]]
int fnSampC(const double & param, NumericVector psi_vec,
               NumericVector w_vec, NumericVector eta_vec,
               const double & upsilon, const double & u_2) {
    // param is either r_tilde_tk or alpha0_j
    const double & L = psi_vec.size();
    NumericVector log_probs(L);

    double phi1;
    double phi2;
    double mean1;
    double mean2;
    for(int ell=0; ell<L; ell++) {
        mean1 = eta_vec(ell);
        mean2 = (upsilon - w_vec(ell)*eta_vec(ell))/(1-w_vec(ell));
        phi1 = R::dnorm(param, mean1, sqrt(u_2), TRUE);
        phi2 = R::dnorm(param, mean2, sqrt(u_2), TRUE);
        log_probs(ell) = log(psi_vec(ell)) + log_sum_exp(log(w_vec(ell)) + phi1, log(1-w_vec(ell)) + phi2);
    }
    NumericVector probs = exp(log_probs - max(log_probs));

    return sample(L, 1, FALSE, probs)(0);
}

// [[Rcpp::export]]
void fnSampCVec(IntegerVector c_vec, IntegerVector d_vec,
                NumericVector param_vec, NumericVector psi_vec,
               NumericVector w_vec, NumericVector eta_vec,
               const double & upsilon, const double & u_2) {
    // param_vec is either r_tilde_vec or alpha0_vec

    int new_c;
    int old_c;
    double param;
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for w_\ell^r; or j for w_\ell^\theta
        param = param_vec(loop_idx);
        old_c = c_vec(loop_idx);
        new_c = fnSampC(param, psi_vec, w_vec, eta_vec, upsilon, u_2);

        c_vec(loop_idx) = new_c;

        if(new_c != old_c) { // Need to update cluster counts
            d_vec(old_c-1) -= 1; // Don't forget about C's 0 indexing here
            d_vec(new_c-1) += 1;
        }
    }
}
// ###

// ### lambda_vec (lambda_tk or lambda_j)
// [[Rcpp::export]]
void fnSampLambdaVec(IntegerVector lambda_vec, NumericVector param_vec,
                     IntegerVector c_vec, NumericVector w_vec,
                     NumericVector eta_vec, const double & upsilon,
                     const double & u_2) {
    // param_vec is either r_tilde_vec or alpha0_vec

    int new_lambda;
    NumericVector log_probs(2);
    NumericVector probs(2);
    double ell_idx0;
    double w_ell;
    double eta_ell;
    double mean1;
    double mean2;
    double p; // Probability that lambda=1
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for w_\ell^r; or j for w_\ell^\theta
        ell_idx0 = c_vec(loop_idx) - 1; // Don't forget about C's 0 indexing here
        w_ell = w_vec(ell_idx0);
        eta_ell = eta_vec(ell_idx0);

        mean1 = eta_ell;
        mean2 = (upsilon - w_ell*eta_ell)/(1-w_ell);

        log_probs(1) = log(w_ell) + R::dnorm(param_vec(loop_idx), mean1, sqrt(u_2), TRUE);
        log_probs(0) = log(1-w_ell) + R::dnorm(param_vec(loop_idx), mean2, sqrt(u_2), TRUE);

        probs = exp(log_probs - max(log_probs));
        p = probs(1)/sum(probs);
        new_lambda = R::rbinom(1, p);

        lambda_vec(loop_idx) = new_lambda;
    }
}
// ###

// ### r_tk
// [[Rcpp::export]]
double fnRtkTildeLogLik(const double & r_tk,
                        NumericVector mu_t_vec, NumericVector s_vec,
                        const int & lambda_tk, const double & eta_ell,
                        const double & w_ell, const double & upsilon_r,
                        const double & u_2, NumericVector y_tk_vec) {
    double prior_part = 0;
    if(lambda_tk == 1) {
        prior_part = R::dnorm(log(r_tk), eta_ell, sqrt(u_2), TRUE);
    }
    else { // lambda_tk == 0
        prior_part = R::dnorm(log(r_tk), (upsilon_r - w_ell*eta_ell)/(1 - w_ell), sqrt(u_2), TRUE);
    }

    double likelihood_part = fnTimepointLogLikNoGamma(mu_t_vec, s_vec, r_tk, y_tk_vec);

    return prior_part + likelihood_part;
}

// [[Rcpp::export]]
void fnSampRVec(NumericVector r_vec, NumericVector r_tilde_vec, NumericMatrix Mu_mat,
                NumericVector s_vec, NumericVector lambda_vec, NumericVector c_vec,
                NumericVector eta_vec, NumericVector w_vec, NumericVector rep_K,
                const double & upsilon_r, const double & u_2, NumericMatrix Y_mat,
                const double & r_tk_tilde_proposal_sd) {
    const int n = rep_K.size();

    int ell_idx0;
    int tk=0;
    double curr_log_lik;
    double prop_log_lik;
    double r_tk_prop;
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            ell_idx0 = c_vec[tk] - 1; // Don't forget about C's 0 indexing here

            r_tk_prop = exp(r_tilde_vec(tk) + R::rnorm(0, r_tk_tilde_proposal_sd));

            curr_log_lik = fnRtkTildeLogLik(r_vec(tk), Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));
            prop_log_lik = fnRtkTildeLogLik(r_tk_prop, Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));

            if(acceptProposal(curr_log_lik, prop_log_lik)) {
                r_vec(tk) = r_tk_prop;
                r_tilde_vec(tk) = log(r_tk_prop);
            }

            tk++;
        }
    }
}

// [[Rcpp::export]]
void fnSampRVecJoint(NumericVector r_vec, NumericVector r_tilde_vec, NumericMatrix Mu_mat,
                     NumericVector s_vec, NumericVector lambda_vec, NumericVector c_vec,
                     NumericVector eta_vec, NumericVector w_vec, NumericVector rep_K,
                     const double & upsilon_r, const double & u_2, NumericMatrix Y_mat,
                     const double & r_tk_tilde_proposal_sd) {
    const int n = rep_K.size();
    const int N = r_vec.size();

    int ell_idx0;
    int tk=0;
    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    NumericVector r_vec_prop(N);
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            ell_idx0 = c_vec[tk] - 1; // Don't forget about C's 0 indexing here

            r_vec_prop(tk) = exp(r_tilde_vec(tk) + R::rnorm(0, r_tk_tilde_proposal_sd));

            curr_log_lik += fnRtkTildeLogLik(r_vec(tk), Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));
            prop_log_lik += fnRtkTildeLogLik(r_vec_prop(tk), Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));

            tk++;
        }
    }

    if(acceptProposal(curr_log_lik, prop_log_lik)) {
        tk=0;
        for(int i=0; i<n; i++) {
            for(int k=0; k<rep_K(i); k++) {
                r_vec(tk) = r_vec_prop(tk);
                r_tilde_vec(tk) = log(r_vec_prop(tk));
            }
        }
    }
}

// [[Rcpp::export]]
void fnSampRVecJointVertical(NumericVector r_vec, NumericVector r_tilde_vec, NumericMatrix Mu_mat,
                             NumericVector s_vec, NumericVector lambda_vec, NumericVector c_vec,
                             NumericVector eta_vec, NumericVector w_vec, NumericVector rep_K,
                             const double & upsilon_r, const double & u_2, NumericMatrix Y_mat,
                             const double & r_tk_tilde_proposal_sd) {
    const int n = rep_K.size();
    const int N = r_vec.size();

    int ell_idx0;
    int tk=0;
    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    NumericVector r_vec_prop(N);
    double vertical_shift = R::rnorm(0, r_tk_tilde_proposal_sd);
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            ell_idx0 = c_vec[tk] - 1; // Don't forget about C's 0 indexing here

            r_vec_prop(tk) = exp(r_tilde_vec(tk) + vertical_shift);

            curr_log_lik += fnRtkTildeLogLik(r_vec(tk), Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));
            prop_log_lik += fnRtkTildeLogLik(r_vec_prop(tk), Mu_mat(i, _), s_vec, lambda_vec(tk), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(tk, _));

            tk++;
        }
    }

    if(acceptProposal(curr_log_lik, prop_log_lik)) {
        tk=0;
        for(int i=0; i<n; i++) {
            for(int k=0; k<rep_K(i); k++) {
                r_vec(tk) = r_vec_prop(tk);
                r_tilde_vec(tk) = log(r_vec_prop(tk));
            }
        }
    }
}
// ###

// ### alpha0_j
// [[Rcpp::export]]
double fnAlpha0jLogLik(const double & alpha0_j, NumericVector r_vec,
                        NumericVector mu_j_vec, const double & s_j,
                        const int & lambda_j, const double & eta_ell,
                        const double & w_ell, const double & upsilon_r,
                        const double & u_2, NumericVector rep_K,
                        NumericVector y_j_vec) {
    double prior_part = 0;
    if(lambda_j == 1) {
        prior_part = R::dnorm(alpha0_j, eta_ell, sqrt(u_2), TRUE);
    }
    else { // lambda_j == 0
        prior_part = R::dnorm(alpha0_j, (upsilon_r - w_ell*eta_ell)/(1 - w_ell), sqrt(u_2), TRUE);
    }

    double likelihood_part = fnOtuLogLikNoGamma(s_j, r_vec, mu_j_vec, rep_K, y_j_vec);

    return prior_part + likelihood_part;
}

// [[Rcpp::export]]
void fnSampAlpha0Vec(NumericVector alpha0_vec, NumericVector r_vec,
                        NumericMatrix Mu_mat, NumericVector s_vec,
                        NumericVector lambda_vec, NumericVector c_vec,
                        NumericVector eta_vec, NumericVector w_vec,
                        const double & upsilon_r, const double & u_2,
                        NumericVector rep_K, NumericMatrix Y_mat,
                        const double & alpha0_j_proposal_sd) {
    const double & J = alpha0_vec.size();

    int ell_idx0;
    NumericVector mu_j_vec_curr(J);
    NumericVector mu_tilde_j_vec_curr(J);
    NumericVector mu_j_vec_prop(J);
    NumericVector mu_tilde_j_vec_prop(J);
    double curr_log_lik;
    double prop_log_lik;
    double alpha0_j_prop;
    for(int j=0; j<J; j++) {
        ell_idx0 = c_vec[j] - 1; // Don't forget about C's 0 indexing here
        alpha0_j_prop = alpha0_vec(j) + R::rnorm(0, alpha0_j_proposal_sd);

        mu_j_vec_curr = Mu_mat(_, j);
        mu_tilde_j_vec_curr = log(mu_j_vec_curr);
        mu_tilde_j_vec_prop = ( mu_tilde_j_vec_curr - alpha0_vec(j) ) + alpha0_j_prop; // Subtract the old theta_0j (scalar) and then add the proposal theta_0j (scalar)
        mu_j_vec_prop = exp(mu_tilde_j_vec_prop);

        curr_log_lik = fnAlpha0jLogLik(alpha0_vec(j), r_vec, mu_j_vec_curr, s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));
        prop_log_lik = fnAlpha0jLogLik(alpha0_j_prop, r_vec, mu_j_vec_prop, s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));

        if(acceptProposal(curr_log_lik, prop_log_lik)) {
            alpha0_vec(j) = alpha0_j_prop;
            Mu_mat(_, j) = mu_j_vec_prop;
        }
    }
}

// [[Rcpp::export]]
void fnSampAlpha0VecJoint(NumericVector alpha0_vec, NumericVector r_vec,
                        NumericMatrix Mu_mat, NumericVector s_vec,
                        NumericVector lambda_vec, NumericVector c_vec,
                        NumericVector eta_vec, NumericVector w_vec,
                        const double & upsilon_r, const double & u_2,
                        NumericVector rep_K, NumericMatrix Y_mat,
                        const double & alpha0_j_proposal_sd) {
    const double & J = alpha0_vec.size();
    const double & n = rep_K.size();

    int ell_idx0;
    NumericMatrix Mu_mat_prop(n, J);
    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    NumericVector alpha0_vec_prop(J);
    for(int j=0; j<J; j++) {
        ell_idx0 = c_vec[j] - 1; // Don't forget about C's 0 indexing here
        alpha0_vec_prop(j) = alpha0_vec(j) + R::rnorm(0, alpha0_j_proposal_sd);

        Mu_mat_prop(_, j) = exp( (log(Mu_mat(_, j)) - alpha0_vec(j)) + alpha0_vec_prop(j) ); // Subtract the old theta_0j (scalar) and then add the proposal theta_0j (scalar)

        curr_log_lik += fnAlpha0jLogLik(alpha0_vec(j), r_vec, Mu_mat(_, j), s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));
        prop_log_lik += fnAlpha0jLogLik(alpha0_vec_prop(j), r_vec, Mu_mat_prop(_, j), s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));
    }
    if(acceptProposal(curr_log_lik, prop_log_lik)) {
        for(int j=0; j<J; j++) {
            alpha0_vec(j) = alpha0_vec_prop(j);
        }
        for(int i=0; i<n; i++) {
            for(int j=0; j<J; j++) {
                Mu_mat(i, j) = Mu_mat_prop(i, j);
            }
        }
    }
}

// [[Rcpp::export]]
void fnSampAlpha0VecJointVertical(NumericVector alpha0_vec, NumericVector r_vec,
                        NumericMatrix Mu_mat, NumericVector s_vec,
                        NumericVector lambda_vec, NumericVector c_vec,
                        NumericVector eta_vec, NumericVector w_vec,
                        const double & upsilon_r, const double & u_2,
                        NumericVector rep_K, NumericMatrix Y_mat,
                        const double & alpha0_j_proposal_sd) {
    const double & J = alpha0_vec.size();
    const double & n = rep_K.size();

    int ell_idx0;
    NumericMatrix Mu_mat_prop(n, J);
    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    NumericVector alpha0_vec_prop(J);
    double vertical_shift = R::rnorm(0, alpha0_j_proposal_sd);
    for(int j=0; j<J; j++) {
        ell_idx0 = c_vec[j] - 1; // Don't forget about C's 0 indexing here
        alpha0_vec_prop(j) = alpha0_vec(j) + vertical_shift;

        Mu_mat_prop(_, j) = exp( (log(Mu_mat(_, j)) - alpha0_vec(j)) + alpha0_vec_prop(j) ); // Subtract the old theta_0j (scalar) and then add the proposal theta_0j (scalar)

        curr_log_lik += fnAlpha0jLogLik(alpha0_vec(j), r_vec, Mu_mat(_, j), s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));
        prop_log_lik += fnAlpha0jLogLik(alpha0_vec_prop(j), r_vec, Mu_mat_prop(_, j), s_vec(j), lambda_vec(j), eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, rep_K, Y_mat(_, j));
    }
    if(acceptProposal(curr_log_lik, prop_log_lik)) {
        for(int j=0; j<J; j++) {
            alpha0_vec(j) = alpha0_vec_prop(j);
        }
        for(int i=0; i<n; i++) {
            for(int j=0; j<J; j++) {
                Mu_mat(i, j) = Mu_mat_prop(i, j);
            }
        }
    }
}
// ###

/* Environmental factors #################################################### */
// [[Rcpp::export]]
double fnThetamjLogLik(const double & theta_mj,
                       const double & tau_j_2, const double & s_j,
                        NumericVector r_vec, NumericVector mu_j_vec,
                        NumericVector rep_K, NumericVector y_j_vec) {
    double prior_part = R::dnorm(theta_mj, 0, sqrt(tau_j_2), TRUE);
    double likelihood_part = fnOtuLogLikNoGamma(s_j, r_vec, mu_j_vec, rep_K, y_j_vec);

    return prior_part + likelihood_part;
}

// [[Rcpp::export]]
void fnSampThetaMat(NumericMatrix Theta_mat, NumericVector tau_2_vec,
                      NumericVector s_vec, NumericVector r_vec,
                      NumericMatrix Mu_mat, NumericVector rep_K, NumericMatrix K_mat,
                      NumericMatrix Y_mat, const double & theta_mj_proposal_sd) {
    const int & M = Theta_mat.nrow();
    const int & J = Theta_mat.ncol();

    NumericVector mu_j_vec_prop;
    double curr_log_lik;
    double prop_log_lik;
    double theta_mj_prop;
    for(int m=0; m<M; m++) {
        for(int j=0; j<J; j++) {
            theta_mj_prop = Theta_mat(m, j) + R::rnorm(0, theta_mj_proposal_sd);
            mu_j_vec_prop = exp( ( log(Mu_mat(_, j)) - K_mat(_ , m)*Theta_mat(m, j) ) + K_mat(_ , m)*theta_mj_prop ); // Subtract out the old kernel part and add the new kernel part

            curr_log_lik = fnThetamjLogLik(Theta_mat(m, j), tau_2_vec(j),  s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j));
            prop_log_lik = fnThetamjLogLik(theta_mj_prop, tau_2_vec(j),  s_vec(j), r_vec, mu_j_vec_prop, rep_K, Y_mat(_, j));

            if(acceptProposal(curr_log_lik, prop_log_lik)) {
                Theta_mat(m, j) = theta_mj_prop;
                Mu_mat(_, j) = mu_j_vec_prop;
            }
        }
    }
}
// ###

// ### sigma_p_2
// [[Rcpp::export]]
void fnSampSigma2Vec(NumericVector sigma_2_vec, NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericVector iota_vec, Function rtigamma, const double & a_sigma, const double & b_sigma) {
    const double & P = sigma_2_vec.size();
    const double & J = Gamma_mat.nrow();

    double sigma_p_2_ubound;
    double alpha_param_part = 0;
    double beta_param_part = 0;
    double tot_nonzero_gamma = 0;
    double tot_squared_beta = 0;
    int j=0;
    for(int p=0; p<P; p++) {
        tot_nonzero_gamma = 0;
        tot_squared_beta = 0;
        for(j=0; j<J; j++) {
            if(Gamma_mat(j, p) != 0) {
                tot_nonzero_gamma++;
            }
            tot_squared_beta += pow(Beta_mat(j, p), 2);
        }

        alpha_param_part = a_sigma + tot_nonzero_gamma/2;
        beta_param_part = b_sigma + tot_squared_beta/2;

        // Find the upper bound for sigma_p_2
        sigma_p_2_ubound = R_PosInf;
        for(j=0; j<J; j++) {
            if( (Gamma_mat(j, p) != 0) && ( pow(Beta_mat(j, p)/iota_vec(p), 2) < sigma_p_2_ubound) ) {
                sigma_p_2_ubound = pow(Beta_mat(j, p)/iota_vec(p), 2);
                // Rcpp::Rcout << "p=" << p << std::endl;
                // Rcpp::Rcout << "sigma_p_2_ubound=" << sigma_p_2_ubound << std::endl;
            }
        }

        // Environment anlpsb("package:anlpsb");
        // Function rtigamma = anlpsb["rtigamma"];
        // Function rtigamma("rtigamma");
        NumericVector draw = rtigamma(1, alpha_param_part, beta_param_part, sigma_p_2_ubound);
        sigma_2_vec(p) = draw(0);
    }
}
// ###

// ### beta_jp and gamma_jp
// [[Rcpp::export]]
NumericVector fnMujVecFromBetajp(const double & beta_jp, const double & beta_jp_curr, NumericVector mu_j_vec_curr, NumericVector z_p_vec) {
    const double & J = mu_j_vec_curr.size();
    NumericVector mu_j_vec(J);
    mu_j_vec = exp( log(mu_j_vec_curr) + z_p_vec*(beta_jp - beta_jp_curr) );

    return mu_j_vec;
}

// [[Rcpp::export]]
double fnBetajpPriorLik(const double & beta_jp, const double & sigma_p_2, const double & iota_p) {
    // Conditional on gamma_jp
    double prior_part;
    if( (beta_jp > iota_p*sqrt(sigma_p_2)) & (beta_jp != 0) ){
        prior_part = R::dnorm(beta_jp, 0.0, sqrt(sigma_p_2), TRUE) - log( 1 - R::pnorm5(iota_p, 0.0, 1.0, TRUE, FALSE) );
    }
    else if( (beta_jp < -iota_p*sqrt(sigma_p_2)) & (beta_jp !=0) ){
        prior_part = R::dnorm(beta_jp, 0.0, sqrt(sigma_p_2), TRUE) - log( R::pnorm5(-iota_p, 0.0, 1.0, TRUE, FALSE) );
    }
    else if(beta_jp == 0) {
        prior_part = log(1);
    }
    else {
        throw std::range_error("Inconsistent beta_jp, iota_p, and sigma_p values");
    }

    return prior_part;
}

// [[Rcpp::export]]
double fnBetajpLogLik(const double & beta_jp, NumericVector mu_j_vec, const double & s_j, NumericVector r_vec, const double & sigma_p_2, NumericVector rep_K, NumericVector y_j_vec, const double & iota_p) {
    double prior_part = fnBetajpPriorLik(beta_jp, sigma_p_2, iota_p);
    double lik_part = fnOtuLogLikNoGamma(s_j, r_vec, mu_j_vec, rep_K, y_j_vec);

    return prior_part + lik_part;
}

// [[Rcpp::export]]
double fnProposeBetajp(int gamma_jp, double sigma_p, double iota_p) {
    double beta_jp_prop = R_PosInf;
    if(gamma_jp == 1) {
        beta_jp_prop = RcppTN::rtn1(0.0, sigma_p, iota_p*sigma_p, R_PosInf);
    }
    else if(gamma_jp == -1) {
        beta_jp_prop = RcppTN::rtn1(0.0, sigma_p, R_NegInf, -iota_p*sigma_p);
    }
    else if(gamma_jp == 0) {
        beta_jp_prop = 0.0;
    }
    else {
        throw std::range_error("No beta_jp proposed because invalid value detected for gamma_jp");
    }

    return beta_jp_prop;
}

// [[Rcpp::export]]
void fnSampBetaMat(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericVector s_vec, NumericVector r_vec, NumericMatrix X_mat, NumericVector sigma_2_vec, NumericVector rep_K, NumericMatrix Y_mat, NumericVector iota_vec, const double & beta_jp_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();
    // const double & n = Mu_mat.nrow();

    double log_curr_to_prop;
    double log_prop_to_curr;
    double curr_log_lik;
    double prop_log_lik;
    double beta_jp_prop;
    NumericVector mu_j_prop;
    for(int p=0; p<P; p++) {
        for(int j=0; j<J; j++) {
            if(Gamma_mat(j, p) != 0) {
                if(Gamma_mat(j, p) == -1) {
                    beta_jp_prop = RcppTN::rtn1(Beta_mat(j, p), beta_jp_proposal_sd, R_NegInf, -iota_vec(p)*sqrt(sigma_2_vec(p)));
                    mu_j_prop = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), Mu_mat(_, j), X_mat(_, p));

                    curr_log_lik = fnBetajpLogLik(Beta_mat(j, p), Mu_mat(_, j), s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));
                    prop_log_lik = fnBetajpLogLik(beta_jp_prop, mu_j_prop, s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));

                    log_curr_to_prop = log(RcppTN::dtn1(beta_jp_prop, Beta_mat(j, p), beta_jp_proposal_sd, R_NegInf, -iota_vec(p)*sqrt(sigma_2_vec(p))));
                    log_prop_to_curr = log(RcppTN::dtn1(Beta_mat(j, p), beta_jp_prop, beta_jp_proposal_sd, R_NegInf, -iota_vec(p)*sqrt(sigma_2_vec(p))));

                    if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
                        Beta_mat(j, p) = beta_jp_prop;
                        Mu_mat(_, j) = mu_j_prop;
                    }
                }
                else if(Gamma_mat(j, p) == 1) {
                    beta_jp_prop = RcppTN::rtn1(Beta_mat(j, p), beta_jp_proposal_sd, iota_vec(p)*sqrt(sigma_2_vec(p)), R_PosInf);
                    mu_j_prop = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), Mu_mat(_, j), X_mat(_, p));

                    curr_log_lik = fnBetajpLogLik(Beta_mat(j, p), Mu_mat(_, j), s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));
                    prop_log_lik = fnBetajpLogLik(beta_jp_prop, mu_j_prop, s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));

                    log_curr_to_prop = log(RcppTN::dtn1(beta_jp_prop, Beta_mat(j, p), beta_jp_proposal_sd, iota_vec(p)*sqrt(sigma_2_vec(p)), R_PosInf));
                    log_prop_to_curr = log(RcppTN::dtn1(Beta_mat(j, p), beta_jp_prop, beta_jp_proposal_sd, iota_vec(p)*sqrt(sigma_2_vec(p)), R_PosInf));

                    if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
                        Beta_mat(j, p) = beta_jp_prop;
                        Mu_mat(_, j) = mu_j_prop;
                    }
                }
                else {
                    throw std::range_error("This should never happen; gamma_jp==0");
                }

            }
        }
    }
}

// [[Rcpp::export]]
double fnProposeBetajpHat(const double & beta_jp_curr, const double & iota_p, const double & sigma_p_2) {
    double beta_jp_hat_prop;
    if(beta_jp_curr == 0.0) {
        beta_jp_hat_prop = RcppTN::rtn1(0.0, sqrt(sigma_p_2), -iota_p*sqrt(sigma_p_2), iota_p*sqrt(sigma_p_2));
    }
    else {
        beta_jp_hat_prop = beta_jp_curr;
    }

    return beta_jp_hat_prop;
}

// [[Rcpp::export]]
double fnLogGammaPrior(int gamma_jp, double pi0_p, double pi1_p) {
    double gamma_prior = R_NegInf;
    if(gamma_jp == 1) {
        // gamma_prior = 1 - pi0_p;
        gamma_prior = log(pi0_p) + log(pi1_p);
    }
    else if(gamma_jp == -1) {
        gamma_prior = log(pi0_p) + log(1-pi1_p);
    }
    else if(gamma_jp == 0) {
        // gamma_prior = pi0_p * pi1_p;
        gamma_prior = log(1 - pi0_p);
    }
    else {
        throw std::range_error("Unable to set gamma prior, invalid value for gamma_jp");
    }

    return gamma_prior;
}

// [[Rcpp::export]]
double fnLogBetajpTransition(const double & beta_jp_hat, const double & gamma_jp, const double & iota_p, const double & sigma_p_2) {
    double log_trans_prob;
    if(gamma_jp != 0) {
        log_trans_prob = log(1);
    }
    else {
        log_trans_prob = log(R::dnorm(beta_jp_hat, 0.0, sqrt(sigma_p_2), FALSE)) - log( 1 - 2*R::pnorm5(-iota_p, 0.0, 1.0, TRUE, FALSE));
    }

    return log_trans_prob;
}

// [[Rcpp::export]]
void fnSampBetaMatGammaMatJoint(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector sigma_2_vec, NumericVector iota_vec, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, NumericMatrix Y_mat, const double & beta_jp_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();

    double curr_log_lik;
    double prop_log_lik;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double log_gamma_prior_curr;
    double log_gamma_prior_prop;

    int gamma_jp_prop;
    double beta_jp_prop;
    double beta_jp_hat;
    double beta_jp_hat_prop;

    NumericVector mu_j_prop;
    for(int p=0; p<P; p++) {
        for(int j=0; j<J; j++) {
            beta_jp_hat = fnProposeBetajpHat(Beta_mat(j, p), iota_vec(p), sigma_2_vec(p));
            beta_jp_hat_prop = R::rnorm(beta_jp_hat, beta_jp_proposal_sd);
            if( fabs(beta_jp_hat_prop) > (iota_vec(p)*sqrt(sigma_2_vec(p))) ) {
                beta_jp_prop = beta_jp_hat_prop;
                gamma_jp_prop = R::sign(beta_jp_hat_prop);
            }
            else {
                beta_jp_prop = 0.0;
                gamma_jp_prop = 0;
            }
            mu_j_prop = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), Mu_mat(_, j), X_mat(_, p));

            log_gamma_prior_curr = fnLogGammaPrior(Gamma_mat(j, p), pi0_vec(p), pi1_vec(p));
            log_gamma_prior_prop = fnLogGammaPrior(gamma_jp_prop, pi0_vec(p), pi1_vec(p));

            curr_log_lik = log_gamma_prior_curr + fnBetajpLogLik(Beta_mat(j, p), Mu_mat(_, j), s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));
            prop_log_lik = log_gamma_prior_prop +  fnBetajpLogLik(beta_jp_prop, mu_j_prop, s_vec(j), r_vec, sigma_2_vec(p), rep_K, Y_mat(_, j), iota_vec(p));

            log_curr_to_prop = fnLogBetajpTransition(beta_jp_hat, Gamma_mat(j, p), iota_vec(p), sigma_2_vec(p));
            log_prop_to_curr = fnLogBetajpTransition(beta_jp_hat, gamma_jp_prop, iota_vec(p), sigma_2_vec(p));

             if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
                Gamma_mat(j, p) = gamma_jp_prop;
                Beta_mat(j, p) = beta_jp_prop;
                Mu_mat(_, j) = mu_j_prop;
            }
        }
    }
}

double sample_int(double x) {
    // Note -- this is already built into Rcpp.  Probably should change to just use the Rcpp version: https://github.com/RcppCore/Rcpp/blob/master/inst/include/Rcpp/sugar/functions/sample.h
    // Sample an integer from {1, 2, ..., x} uniformly
    double u = R::runif(0,1);
    double i=1;
    for(i=1; i<=x; i++) {
        if(u < i/x) {
            break;
        }
    }

    return i;
}

// [[Rcpp::export]]
void fnBetaGammaAddSwpDel(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericVector sigma_2_vec, NumericVector iota_vec, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, NumericMatrix Y_mat) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();

    double p1_idx0; // 0-indexed, different from math notation!
    double gamma_jp1;
    double num_p2_choices;
    double p2_choice_idx0;
    double p2_idx0;
    double gamma_jp2;
    NumericVector p2_choices(P);
    double gamma_jp1_prime;
    double gamma_jp2_prime;
    double beta_jp1_prime;
    double beta_jp2_prime;
    NumericVector mu_j_prop;
    double curr_log_lik;
    double prop_log_lik;
    const int num_swaps=10;
    for(int j=0; j<J; j++) {
        for(int k=0; k<num_swaps; k++) { // Inner loop here to increase the number of swaps
            p1_idx0 = sample_int(P) - 1;
            gamma_jp1 = Gamma_mat(j, p1_idx0);

            std::fill(p2_choices.begin(), p2_choices.end(), -100);
            num_p2_choices = 0;
            for(int p=0; p<P; p++) {
                if(gamma_jp1 != Gamma_mat(j, p)) {
                    p2_choices(num_p2_choices) = p;
                    num_p2_choices++;
                }
            }

            if(num_p2_choices > 0) {
                p2_choice_idx0 = sample_int(num_p2_choices) - 1;
                p2_idx0 = p2_choices(p2_choice_idx0);
                gamma_jp2 = Gamma_mat(j, p2_idx0);

                gamma_jp1_prime = gamma_jp2;
                gamma_jp2_prime = gamma_jp1;

                beta_jp1_prime = fnProposeBetajp(gamma_jp1_prime, sqrt(sigma_2_vec(p1_idx0)), iota_vec(p1_idx0));
                beta_jp2_prime = fnProposeBetajp(gamma_jp2_prime, sqrt(sigma_2_vec(p2_idx0)), iota_vec(p2_idx0));

                // Need to do this twice since we change two beta_jp's
                mu_j_prop = fnMujVecFromBetajp(beta_jp1_prime, Beta_mat(j, p1_idx0), Mu_mat(_, j), X_mat(_, p1_idx0));
                mu_j_prop = fnMujVecFromBetajp(beta_jp2_prime, Beta_mat(j, p2_idx0), mu_j_prop, X_mat(_, p2_idx0));

                curr_log_lik = fnOtuLogLikNoGamma(s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j)) + fnLogGammaPrior(gamma_jp1, pi0_vec(p1_idx0), pi1_vec(p1_idx0)) + fnLogGammaPrior(gamma_jp2, pi0_vec(p2_idx0), pi1_vec(p2_idx0));
                prop_log_lik = fnOtuLogLikNoGamma(s_vec(j), r_vec, mu_j_prop, rep_K, Y_mat(_, j)) + fnLogGammaPrior(gamma_jp1_prime, pi0_vec(p1_idx0), pi1_vec(p1_idx0)) + fnLogGammaPrior(gamma_jp2_prime, pi0_vec(p2_idx0), pi1_vec(p2_idx0));

                if(acceptProposal(curr_log_lik, prop_log_lik)) {
                    Gamma_mat(j, p1_idx0) = gamma_jp1_prime;
                    Beta_mat(j, p1_idx0) = beta_jp1_prime;
                    Gamma_mat(j, p2_idx0) = gamma_jp2_prime;
                    Beta_mat(j, p2_idx0) = beta_jp2_prime;
                    Mu_mat(_, j) = mu_j_prop;
                }
            }
            else {
                // Rcpp::Rcout << "No swap options" << ", j=" << j << ", p1_idx0=" << p1_idx0 << std::endl;
            }
        }
    }
}
// ###

// ### iota_p
// [[Rcpp::export]]
double fnIotapLogLik(const double & iota_p, NumericVector beta_p_vec, NumericVector gamma_p_vec, const double & sigma_p_2, const double & smallest_beta, const double & a_iota, const double & b_iota) {
    if(smallest_beta < 0) throw std::range_error("Negative 'smallest_beta' not handled");
    if(iota_p > smallest_beta/sqrt(sigma_p_2)) return R_NegInf; // Check for proper support of iota_p
    const double & J = beta_p_vec.size();

    double prior_part = (a_iota - 1)*log(iota_p) - b_iota*iota_p;
    double lik_part = 0.0;
    for(int j=0; j<J; j++) {
        if(gamma_p_vec(j) == -1) {
            lik_part += log(1/R::pnorm5(-iota_p, 0.0, 1.0, TRUE, FALSE));
        }
        else if(gamma_p_vec(j) == 1) {
            lik_part += log(1/(1 - R::pnorm5(iota_p, 0.0, 1.0, TRUE, FALSE)));
        }
    }

    return prior_part + lik_part;
}

// [[Rcpp::export]]
void fnSampIotaVec(NumericVector iota_vec, NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericVector sigma_2_vec, const double & a_iota, const double & b_iota, const double & iota_p_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();

    double curr_log_lik;
    double prop_log_lik;
    double iota_p_prop;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double smallest_beta; // Should always be positive
    for(int p=0; p<P; p++) {
        // Determine the support of iota_p
        smallest_beta = R_PosInf;
        for(int j=0; j<J; j++) {
            if( (Gamma_mat(j, p) != 0) & (std::abs(Beta_mat(j, p)) < smallest_beta) ) {
                smallest_beta = std::abs(Beta_mat(j, p));
            }
        }
        iota_p_prop = RcppTN::rtn1(iota_vec(p), iota_p_proposal_sd, 0, smallest_beta/sqrt(sigma_2_vec(p)));

        curr_log_lik = fnIotapLogLik(iota_vec(p), Beta_mat(_, p), Gamma_mat(_, p), sigma_2_vec(p), smallest_beta, a_iota, b_iota);
        prop_log_lik = fnIotapLogLik(iota_p_prop, Beta_mat(_, p), Gamma_mat(_, p), sigma_2_vec(p), smallest_beta, a_iota, b_iota);
        log_curr_to_prop = log(RcppTN::dtn1(iota_p_prop, iota_vec(p), iota_p_proposal_sd, 0, smallest_beta/sqrt(sigma_2_vec(p))));
        log_prop_to_curr = log(RcppTN::dtn1(iota_vec(p), iota_p_prop, iota_p_proposal_sd, 0, smallest_beta/sqrt(sigma_2_vec(p))));

        if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
            iota_vec(p) = iota_p_prop;
        }
    }
}
// ###

// ### Z (impute missing covarites)
void fnGenMuProp(NumericMatrix Mu_prop_mat, NumericVector mu_t_vec, NumericMatrix Beta_mat, const double & curr_cat, const double & max_cat, const double & p_missing_idx0) {
    const double & J = Beta_mat.nrow();

    double curr_z_eff;
    double prop_z_eff;
    for(int j=0; j<J; j++) {
        // Covariate effect
        if(curr_cat == 0) {
            curr_z_eff = 0; // Handle the case where all the indicators are 0
        }
        else {
            curr_z_eff = Beta_mat(j, p_missing_idx0 + (curr_cat-1)); // Subtract 1 to handle the case where all the indicators are 0
        }

        for(int l=0; l<=max_cat; l++) {
            if(l==0) {
                prop_z_eff = 0; // Handle the case where all the indicators are 0
            }
            else {
                prop_z_eff = Beta_mat(j, p_missing_idx0 + (l-1)); // Subtract 1 to handle the case where all the indicators are 0
            }
            Mu_prop_mat(j, l) = exp( log(mu_t_vec(j)) - curr_z_eff + prop_z_eff);
        }
    }
}

// [[Rcpp::export]]
double fnGetCurrCat(NumericVector categorical_vector) {
    // Given a 1-hot-encoded vector, returns the category
    // e.g. c(0, 0, 0) returns 0
    //      c(0, 1, 0) returns 2
    //      c(1, 0, 0) returns 1
    for(int i=0; i<categorical_vector.size(); i++) {
        if(categorical_vector(i) == 1) return (i+1);
    }

    return 0;
}

// [[Rcpp::export]]
void fnImputeZ(NumericMatrix X_mat, NumericMatrix Miss_var_ind_mat, NumericVector max_miss_cat_vec, NumericVector first_p_missing_idx1_vec, NumericMatrix Mu_mat, NumericMatrix Beta_mat, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector rep_K_cumsum, NumericMatrix Y_mat) {
    const double & J = Mu_mat.ncol();
    const double & n_missing_vars = Miss_var_ind_mat.ncol(); // Number of missing covariates
    const double & n = Miss_var_ind_mat.nrow(); // Number of timepoints

    int tk;
    double max_cat;
    double p_missing_idx0;
    double curr_cat;
    double new_cat;
    for(int b=0; b<n_missing_vars; b++) {
        p_missing_idx0 = first_p_missing_idx1_vec(b) - 1;
        max_cat = max_miss_cat_vec(b); // The number of possible categories the missing covariate can take
        NumericVector curr_categorical_vec(max_cat); // One-hot encoded vector for the current category
        NumericMatrix Mu_prop_mat(J, max_cat+1); // Need to add 1 more for the category where all the indicators are 0
        for(int i=0; i<n; i++) {
            if(Miss_var_ind_mat(i, b) == 1) {
                // Find the current one-hot encoded category
                for(int l=0; l<max_cat; l++) {
                    curr_categorical_vec(l) = X_mat(i, p_missing_idx0 + l);
                }
                curr_cat = fnGetCurrCat(curr_categorical_vec);

                // Get Mu_j proposals depending on the imputed category
                fnGenMuProp(Mu_prop_mat, Mu_mat(i, _), Beta_mat, curr_cat, max_cat, p_missing_idx0); // Modifies Mu_prop_mat in place

                NumericVector log_probs(max_cat + 1);
                tk=rep_K_cumsum(i);
                for(int l=0; l<(max_cat+1); l++) {
                    for(int k=0; k<rep_K(i); k++) {
                        for(int j=0; j<J; j++) {
                            log_probs(l) += fnLogLikNoGamma(Y_mat(tk, j), s_vec(j), r_vec(tk), Mu_prop_mat(j, l));
                        }
                        tk++;
                    }
                }
                NumericVector probs = exp(log_probs - max(log_probs));
                probs = probs/sum(probs); // Normalize

                // Sample new category (new category may be 0)
                new_cat = sample(max_cat+1, 1, FALSE, probs)(0) - 1;

                for(int j=0; j<J; j++) {
                    Mu_mat(i, j) = Mu_prop_mat(j, new_cat);
                }

                // Update the covariate matrix
                for(int l=0; l<max_cat; l++) {
                    X_mat(i, p_missing_idx0 + l) = 0;
                    if( (new_cat-1)==l ) {
                        X_mat(i, p_missing_idx0 + l) = 1;
                    }
                }
            }
        }
    }
}

// [[Rcpp::export]]
double fnDLogInvGamma(const double & x, const double & a, const double & b) {
    return a*log(b) - lgamma(a) + (-a-1)*log(x) - b/x;
}

// [[Rcpp::export]]
double fnIotapSigmapJointLogLik(const double & iota_p, const double & sigma_p_2, NumericVector beta_p_vec, NumericVector gamma_p_vec, NumericMatrix Mu_mat, NumericVector s_vec, NumericVector r_vec, const double & pi0_p, const double & pi1_p, NumericVector rep_K, const double & a_iota, const double & b_iota, const double & a_sigma, const double & b_sigma, NumericMatrix Y_mat) {
    // Guard for proper support
    if( (iota_p < 0.0) || (sigma_p_2 < 0.0) ) {
        return R_NegInf;
    }

    const double & J = beta_p_vec.size();

    double log_iota_prior = R::dgamma(iota_p, a_iota, 1/b_iota, TRUE);
    double log_sigma_p_2_prior = fnDLogInvGamma(sigma_p_2, a_sigma, b_sigma);

    double log_gamma_prior = 0;
    double log_lik_beta = 0;
    for(int j=0; j<J; j++) {
        log_gamma_prior += fnLogGammaPrior(gamma_p_vec(j), pi0_p, pi1_p);
        log_lik_beta += fnBetajpLogLik(beta_p_vec(j), Mu_mat(_, j), s_vec(j), r_vec, sigma_p_2, rep_K, Y_mat(_, j), iota_p);
    }

    return log_iota_prior + log_sigma_p_2_prior + log_gamma_prior + log_lik_beta;
}

// [[Rcpp::export]]
double fnConvertBetaHat(double & beta_jp_hat, const double & iota_p, const double & sigma_p_2) {
    double beta_jp = R_PosInf;
    if(fabs(beta_jp_hat) <= iota_p * sqrt(sigma_p_2)) {
        beta_jp = 0.0;
    }
    else {
        beta_jp = beta_jp_hat;
    }

    return beta_jp;
}

// [[Rcpp::export]]
void fnSampIotaVecSigma2VecJoint(NumericVector sigma_2_vec, NumericVector iota_vec, NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, const double & a_sigma, const double & b_sigma, const double & a_iota, const double & b_iota, NumericMatrix Y_mat, const double & iota_p_joint_proposal_sd, const double & sigma_p_2_joint_proposal_sd) {
    const double & J = Beta_mat.nrow();
    const double & P = Beta_mat.ncol();
    const double & n = Mu_mat.nrow();

    NumericMatrix Mu_mat_prop(n, J);
    double curr_log_lik;
    double prop_log_lik;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double beta_jp_hat;
    NumericVector beta_p_vec_prop(J);
    NumericVector gamma_p_vec_prop(J);
    double iota_p_prop;
    double sigma_p_2_prop;
    for(int p=0; p<P; p++) {
        iota_p_prop = iota_vec(p) + R::rnorm(0, iota_p_joint_proposal_sd);
        sigma_p_2_prop = sigma_2_vec(p) + R::rnorm(0, sigma_p_2_joint_proposal_sd);

        // Copy Mu_mat into Mu_mat_prop
        for(int i=0; i<n; i++) {
            for(int j=0; j<J; j++) {
                Mu_mat_prop(i, j) = Mu_mat(i, j);
            }
        }

        log_curr_to_prop = 0.0;
        log_prop_to_curr = 0.0;
        for(int j=0; j<J; j++) {
            beta_jp_hat = fnProposeBetajpHat(Beta_mat(j, p), iota_vec(p), sigma_2_vec(p));
            beta_p_vec_prop(j) = fnConvertBetaHat(beta_jp_hat, iota_p_prop, sigma_p_2_prop); // Modified in place

            Mu_mat_prop(_, j) = fnMujVecFromBetajp(beta_p_vec_prop(j), Beta_mat(j, p), Mu_mat_prop(_, j), X_mat(_, p));
            gamma_p_vec_prop(j) = R::sign(beta_p_vec_prop(j));
            log_curr_to_prop += fnLogBetajpTransition(beta_jp_hat, Gamma_mat(j, p), iota_vec(p), sigma_2_vec(p));
            log_prop_to_curr += fnLogBetajpTransition(beta_jp_hat, gamma_p_vec_prop(j), iota_p_prop, sigma_p_2_prop);
        }
        curr_log_lik = fnIotapSigmapJointLogLik(iota_vec(p), sigma_2_vec(p), Beta_mat(_, p), Gamma_mat(_, p), Mu_mat, s_vec, r_vec, pi0_vec(p), pi1_vec(p), rep_K, a_iota, b_iota, a_sigma, b_sigma, Y_mat);
        prop_log_lik = fnIotapSigmapJointLogLik(iota_p_prop, sigma_p_2_prop, beta_p_vec_prop, gamma_p_vec_prop, Mu_mat_prop, s_vec, r_vec, pi0_vec(p), pi1_vec(p), rep_K, a_iota, b_iota, a_sigma, b_sigma, Y_mat);

        if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
            iota_vec(p) = iota_p_prop;
            sigma_2_vec(p) = sigma_p_2_prop;
            Beta_mat(_, p) = beta_p_vec_prop;
            Gamma_mat(_, p) = gamma_p_vec_prop;
            for(int i=0; i<n; i++) {
                for(int j=0; j<J; j++) {
                    Mu_mat(i, j) = Mu_mat_prop(i, j);
                }
            }
        }
    }
}

// [[Rcpp::export]]
void fnSampBetaMatGammaMatJointWholeOtu(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector sigma_2_vec, NumericVector iota_vec, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, NumericMatrix Y_mat, const double & beta_jp_joint_otu_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();
    const double & n = Mu_mat.nrow();

    double curr_log_lik;
    double prop_log_lik;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double log_gamma_prior_curr;
    double log_gamma_prior_prop;
    double log_beta_prior_curr;
    double log_beta_prior_prop;

    int gamma_jp_prop;
    double beta_jp_prop;
    double beta_jp_hat;
    double beta_jp_hat_prop;

    NumericVector mu_j_prop(n);
    NumericVector beta_j_vec_prop(P);
    NumericVector gamma_j_vec_prop(P);
    for(int j=0; j<J; j++) {
        for(int i=0; i<n; i++) mu_j_prop(i) = Mu_mat(i, j);
        log_curr_to_prop = 0.0;
        log_prop_to_curr = 0.0;
        log_gamma_prior_curr = 0.0;
        log_gamma_prior_prop = 0.0;
        log_beta_prior_curr = 0.0;
        log_beta_prior_prop = 0.0;
        for(int p=0; p<P; p++) {
            beta_jp_hat = fnProposeBetajpHat(Beta_mat(j, p), iota_vec(p), sigma_2_vec(p));
            beta_jp_hat_prop = R::rnorm(beta_jp_hat, beta_jp_joint_otu_proposal_sd);
            beta_jp_prop = fnConvertBetaHat(beta_jp_hat_prop, iota_vec(p), sigma_2_vec(p));
            gamma_jp_prop = R::sign(beta_jp_prop);
            beta_j_vec_prop(p) = beta_jp_prop;
            gamma_j_vec_prop(p) = gamma_jp_prop;

            log_gamma_prior_curr += fnLogGammaPrior(Gamma_mat(j, p), pi0_vec(p), pi1_vec(p));
            log_gamma_prior_prop += fnLogGammaPrior(gamma_jp_prop, pi0_vec(p), pi1_vec(p));
            log_beta_prior_curr += fnBetajpPriorLik(Beta_mat(j, p), sigma_2_vec(p), iota_vec(p));
            log_beta_prior_prop += fnBetajpPriorLik(beta_jp_prop, sigma_2_vec(p), iota_vec(p));

            log_curr_to_prop += fnLogBetajpTransition(beta_jp_hat, Gamma_mat(j, p), iota_vec(p), sigma_2_vec(p));
            log_prop_to_curr += fnLogBetajpTransition(beta_jp_hat, gamma_jp_prop, iota_vec(p), sigma_2_vec(p));

            mu_j_prop = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), mu_j_prop, X_mat(_, p));
        }
        curr_log_lik = log_gamma_prior_curr + log_beta_prior_curr + fnOtuLogLikNoGamma(s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j));
        prop_log_lik = log_gamma_prior_prop + log_beta_prior_prop + fnOtuLogLikNoGamma(s_vec(j), r_vec, mu_j_prop, rep_K, Y_mat(_, j));

        if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
            Beta_mat(j, _) = beta_j_vec_prop;
            Gamma_mat(j, _) = gamma_j_vec_prop;
            for(int i=0; i<n; i++) {
                Mu_mat(i, j) = mu_j_prop(i);
            }
        }
    }
}

// [[Rcpp::export]]
void fnSampBetaMatGammaMatJointWholeOtu2(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector sigma_2_vec, NumericVector iota_vec, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, NumericMatrix Y_mat, const double & beta_jp_joint_otu_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();
    const double & n = Mu_mat.nrow();

    double curr_log_lik;
    double prop_log_lik;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double log_gamma_prior_curr;
    double log_gamma_prior_prop;
    double log_beta_prior_curr;
    double log_beta_prior_prop;

    int gamma_jp_prop;
    double beta_jp_prop;
    double beta_jp_hat;
    double beta_jp_hat_prop;
    double vertical_shift;

    NumericVector mu_j_prop(n);
    NumericVector beta_j_vec_prop(P);
    NumericVector gamma_j_vec_prop(P);
    for(int j=0; j<J; j++) {
        for(int i=0; i<n; i++) mu_j_prop(i) = Mu_mat(i, j);
        log_curr_to_prop = 0.0;
        log_prop_to_curr = 0.0;
        log_gamma_prior_curr = 0.0;
        log_gamma_prior_prop = 0.0;
        log_beta_prior_curr = 0.0;
        log_beta_prior_prop = 0.0;
        vertical_shift =  R::rnorm(0.0, beta_jp_joint_otu_proposal_sd);
        for(int p=0; p<P; p++) {
            beta_jp_hat = fnProposeBetajpHat(Beta_mat(j, p), iota_vec(p), sigma_2_vec(p));
            beta_jp_hat_prop = beta_jp_hat + vertical_shift;
            beta_jp_prop = fnConvertBetaHat(beta_jp_hat_prop, iota_vec(p), sigma_2_vec(p));
            gamma_jp_prop = R::sign(beta_jp_prop);
            beta_j_vec_prop(p) = beta_jp_prop;
            gamma_j_vec_prop(p) = gamma_jp_prop;

            log_gamma_prior_curr += fnLogGammaPrior(Gamma_mat(j, p), pi0_vec(p), pi1_vec(p));
            log_gamma_prior_prop += fnLogGammaPrior(gamma_jp_prop, pi0_vec(p), pi1_vec(p));
            log_beta_prior_curr += fnBetajpPriorLik(Beta_mat(j, p), sigma_2_vec(p), iota_vec(p));
            log_beta_prior_prop += fnBetajpPriorLik(beta_jp_prop, sigma_2_vec(p), iota_vec(p));

            log_curr_to_prop += fnLogBetajpTransition(beta_jp_hat, Gamma_mat(j, p), iota_vec(p), sigma_2_vec(p));
            log_prop_to_curr += fnLogBetajpTransition(beta_jp_hat, gamma_jp_prop, iota_vec(p), sigma_2_vec(p));

            mu_j_prop = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), mu_j_prop, X_mat(_, p));
        }
        curr_log_lik = log_gamma_prior_curr + log_beta_prior_curr + fnOtuLogLikNoGamma(s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j));
        prop_log_lik = log_gamma_prior_prop + log_beta_prior_prop + fnOtuLogLikNoGamma(s_vec(j), r_vec, mu_j_prop, rep_K, Y_mat(_, j));

        if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
            Beta_mat(j, _) = beta_j_vec_prop;
            Gamma_mat(j, _) = gamma_j_vec_prop;
            for(int i=0; i<n; i++) {
                Mu_mat(i, j) = mu_j_prop(i);
            }
        }
    }
}

// [[Rcpp::export]]
void fnSampBetaMatGammaMatJointWholeOtu3(NumericMatrix Beta_mat, NumericMatrix Gamma_mat, NumericMatrix Mu_mat, NumericMatrix X_mat, NumericVector sigma_2_vec, NumericVector iota_vec, NumericVector s_vec, NumericVector r_vec, NumericVector rep_K, NumericVector pi0_vec, NumericVector pi1_vec, NumericMatrix Y_mat, const double & beta_jp_joint_otu_proposal_sd) {
    const double & P = Beta_mat.ncol();
    const double & J = Beta_mat.nrow();
    const double & n = Mu_mat.nrow();

    double curr_log_lik;
    double prop_log_lik;
    double log_curr_to_prop;
    double log_prop_to_curr;
    double log_gamma_prior_curr;
    double log_gamma_prior_prop;
    double log_beta_prior_curr;
    double log_beta_prior_prop;

    int gamma_jp_prop;
    double beta_jp_prop;
    double beta_jp_hat;
    double beta_jp_hat_prop;

    NumericMatrix Mu_prop_mat(n, J);
    NumericVector beta_p_vec_prop(J);
    NumericVector gamma_p_vec_prop(J);
    for(int p=0; p<P; p++) {
        log_curr_to_prop = 0.0;
        log_prop_to_curr = 0.0;
        log_gamma_prior_curr = 0.0;
        log_gamma_prior_prop = 0.0;
        log_beta_prior_curr = 0.0;
        log_beta_prior_prop = 0.0;
        curr_log_lik = 0.0;
        prop_log_lik = 0.0;
        double vertical_shift;

        for(int i=0; i<n; i++) {
            for(int j=0; j<J; j++) {
                Mu_prop_mat(i, j) = Mu_mat(i, j);
            }
        }

        vertical_shift = R::rnorm(0.0, beta_jp_joint_otu_proposal_sd);
        for(int j=0; j<J; j++) {
            beta_jp_hat = fnProposeBetajpHat(Beta_mat(j, p), iota_vec(p), sigma_2_vec(p));
            beta_jp_hat_prop = beta_jp_hat + vertical_shift;
            beta_jp_prop = fnConvertBetaHat(beta_jp_hat_prop, iota_vec(p), sigma_2_vec(p));
            gamma_jp_prop = R::sign(beta_jp_prop);
            beta_p_vec_prop(j) = beta_jp_prop;
            gamma_p_vec_prop(j) = gamma_jp_prop;

            log_gamma_prior_curr += fnLogGammaPrior(Gamma_mat(j, p), pi0_vec(p), pi1_vec(p));
            log_gamma_prior_prop += fnLogGammaPrior(gamma_jp_prop, pi0_vec(p), pi1_vec(p));
            log_beta_prior_curr += fnBetajpPriorLik(Beta_mat(j, p), sigma_2_vec(p), iota_vec(p));
            log_beta_prior_prop += fnBetajpPriorLik(beta_jp_prop, sigma_2_vec(p), iota_vec(p));

            log_curr_to_prop += fnLogBetajpTransition(beta_jp_hat, Gamma_mat(j, p), iota_vec(p), sigma_2_vec(p));
            log_prop_to_curr += fnLogBetajpTransition(beta_jp_hat, gamma_jp_prop, iota_vec(p), sigma_2_vec(p));

            Mu_prop_mat(_, j) = fnMujVecFromBetajp(beta_jp_prop, Beta_mat(j, p), Mu_prop_mat(_, j), X_mat(_, p));

            curr_log_lik += log_gamma_prior_curr + log_beta_prior_curr + fnOtuLogLikNoGamma(s_vec(j), r_vec, Mu_mat(_, j), rep_K, Y_mat(_, j));
            prop_log_lik += log_gamma_prior_prop + log_beta_prior_prop + fnOtuLogLikNoGamma(s_vec(j), r_vec, Mu_prop_mat(_, j), rep_K, Y_mat(_, j));
        }

        if(acceptProposal(curr_log_lik, prop_log_lik, log_curr_to_prop, log_prop_to_curr)) {
            Beta_mat(_, p) = beta_p_vec_prop;
            Gamma_mat(_, p) = gamma_p_vec_prop;
            for(int i=0; i<n; i++) {
                for(int j=0; j<J; j++) {
                    Mu_mat(i, j) = Mu_prop_mat(i, j);
                }
            }
        }
    }
}
// ###
