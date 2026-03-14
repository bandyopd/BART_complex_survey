# ==============================================================================
# ComplexSurvey.R
#
# Soft Bayesian Additive Regression Trees (SBART) for clustered complex survey
# data with non-Gaussian errors.
#
# Implements three models described in:
#   Mandal, Linero, Bandyopadhyay & Sinha (2025).
#   "Soft Bayesian Additive Regression Trees (SBART) for correlated survey
#    response with non-Gaussian error."
#   Journal of Nonparametric Statistics.
#   DOI: 10.1080/10485252.2025.2574708
# ==============================================================================


# ==============================================================================
# 1. Internal Utility Functions (not exported)
# ==============================================================================

# Normalised rank transformation mapping a vector to [0, 1].
trank <- function(x) {
  x_unique <- unique(x)
  x_ranks  <- rank(x_unique, ties.method = "max")
  tx <- x_ranks[match(x, x_unique)] - 1
  tx <- tx / length(unique(tx))
  tx <- tx / max(tx)
  return(tx)
}

# Applies trank() column-wise to a predictor matrix X.
quantile_normalize_bart <- function(X) {
  apply(X = X, MARGIN = 2, trank)
}

# Jointly rank-normalises training and test predictors.
quantile_normalize <- function(X, test_X) {
  X_trans   <- quantile_normalize_bart(rbind(X, test_X))
  idx_train <- 1:nrow(X)
  X         <- X_trans[ idx_train, , drop = FALSE]
  test_X    <- X_trans[-idx_train, , drop = FALSE]
  return(list(X = X, test_X = test_X))
}

# Computes the sum of a vector within consecutive subgroups of given sizes.
sum_subgroups <- function(vector, subgroup_sizes) {
  sums        <- numeric(length(subgroup_sizes))
  start_index <- 1
  for (i in seq_along(subgroup_sizes)) {
    end_index   <- start_index + subgroup_sizes[i] - 1
    sums[i]     <- sum(vector[start_index:end_index])
    start_index <- end_index + 1
  }
  return(sums)
}

# Convenience wrapper around predict.glmnet.
predict_glmnet <- function(object, newx, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s      <- match.arg(s)
    lambda <- object[[s]]
    names(lambda) <- s
  } else {
    stop("Invalid form for s")
  }
  predict(object$glmnet.fit, newx, s = lambda, ...)
}

# Expands subject-level predictor matrix to observation level.
create_X <- function(x, ni) {
  tmp <- c(rep(0, ncol(x)))
  for (i in 1:nrow(x)) {
    for (j in 1:ni[i]) {
      tmp <- rbind(tmp, x[i, ])
    }
  }
  tmp <- tmp[-1, ]
  return(tmp)
}


# ==============================================================================
# 2. Exported Functions
# ==============================================================================

#' Estimate residual SD using cross-validated LASSO
#'
#' Fits a cross-validated LASSO model and returns the root mean squared
#' prediction error as an estimate of the residual standard deviation.
#' Useful for providing a data-driven starting value or prior scale for
#' sigma_0 in the SBART models.
#'
#' @param X numeric matrix or data.frame of predictors (n x p)
#' @param Y numeric response vector of length n
#' @param weights optional numeric weight vector of length n (default: all ones)
#'
#' @return Scalar estimate of the residual standard deviation
#' @export
GetSigma <- function(X, Y, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(Y))
  stopifnot(is.matrix(X) | is.data.frame(X))
  if (is.data.frame(X)) X <- model.matrix(~ . - 1, data = X)
  fit       <- glmnet::cv.glmnet(x = X, y = Y, weights = weights)
  fitted    <- predict_glmnet(fit, X)
  sigma_hat <- sqrt(mean((fitted - Y)^2))
  return(sigma_hat)
}


#' SBART with Skew-Normal Error for Complex Survey Data
#'
#' Fits a nonparametric Bayesian regression model with Skew-Normal errors
#' for clustered complex survey data with subject-specific survey weights.
#' The nonparametric regression function is modelled using Soft Bayesian
#' Additive Regression Trees (SBART; Linero & Yang 2018). Subject-specific
#' random intercepts account for within-cluster correlation.
#'
#' The hierarchical model is (Section 2.3 of Mandal et al. 2025):
#' \deqn{y_{ij} | b, r_i, Z_{ij}^+ \sim N(b(x_i) + r_i + \lambda Z_{ij}^+,
#'       \gamma^2 / w_i)}
#' \deqn{Z_{ij}^+ \sim N^+(0, 1/w_i), \quad r_i \sim N(0, \sigma_1^2/w_i)}
#' \deqn{\gamma \sim \text{Cauchy}^+(0, \gamma_0), \quad
#'       \lambda \sim N(0, d^2), \quad \sigma_1^2 \sim IG(a_1, c_1)}
#'
#' @param X numeric matrix of training predictors (n_obs x p)
#' @param Y numeric response vector of length n_obs
#' @param test_X numeric matrix of test predictors
#' @param hypers SoftBart hyperparameters object from \code{SoftBart::Hypers()};
#'   uses default if NULL
#' @param opts SoftBart options object from \code{SoftBart::Opts()};
#'   uses default if NULL
#' @param w numeric vector of survey weights of length n_obs
#' @param ni integer vector of cluster sizes of length n_subjects;
#'   must satisfy \code{sum(ni) == n_obs}
#'
#' @return A named list containing:
#' \describe{
#'   \item{y_hat_train}{posterior samples of conditional mean for training data
#'     (matrix: n_save x n_obs)}
#'   \item{y_hat_test}{posterior samples of conditional mean for test data}
#'   \item{y_hat_train_mean}{posterior mean of conditional mean, training}
#'   \item{y_hat_test_mean}{posterior mean of conditional mean, test}
#'   \item{f_hat_train}{posterior samples of b(x)+r (location), training}
#'   \item{f_hat_test}{posterior samples of b(x) (location), test}
#'   \item{f_hat_train_mean}{posterior mean of b(x)+r, training}
#'   \item{f_hat_test_mean}{posterior mean of b(x), test}
#'   \item{sigma}{posterior samples of total scale sigma}
#'   \item{tau}{posterior samples of symmetric scale gamma}
#'   \item{alpha}{posterior samples of skewness parameter alpha}
#'   \item{lambda}{posterior samples of lambda}
#'   \item{varcounts}{variable inclusion counts per iteration}
#'   \item{likelihood_mat}{log-likelihood matrix (n_save x n_obs) for LOO/LPML}
#'   \item{r}{final MCMC draw of random effects (expanded to obs. level)}
#'   \item{rsaved}{all MCMC draws of r (matrix: n_iter x n_obs)}
#' }
#'
#' @references
#' Mandal A, Linero AR, Bandyopadhyay D, Sinha D (2025).
#' Soft Bayesian Additive Regression Trees (SBART) for correlated survey
#' response with non-Gaussian error.
#' \emph{Journal of Nonparametric Statistics}.
#' \doi{10.1080/10485252.2025.2574708}
#'
#' @examples
#' \dontrun{
#' library(SoftBart)
#' dat    <- sim_fried_SN(N = 100, P = 5, alpha = 5, sigma = 1, sigma1 = 0.2)
#' hypers <- SoftBart::Hypers(dat$X, dat$Y, num_tree = 50)
#' opts   <- SoftBart::Opts(num_burn = 500, num_save = 500)
#' fit    <- SNSBARTW(X = dat$X, Y = dat$Y, test_X = dat$X,
#'                    w = dat$w, ni = dat$ni,
#'                    hypers = hypers, opts = opts)
#' sqrt(mean(dat$w * (dat$Y - fit$y_hat_train_mean)^2))
#' }
#'
#' @export
SNSBARTW <- function(X, Y, test_X, hypers = NULL, opts = NULL, w, ni) {

  if (is.null(opts))   opts   <- SoftBart::Opts()
  if (is.null(hypers)) hypers <- SoftBart::Hypers(X, Y)

  iter  <- opts$num_burn + opts$num_save
  burn  <- opts$num_burn
  opts  <- SoftBart::Opts()

  n_train     <- nrow(X)
  mu_hat      <- matrix(NA, iter, n_train)
  mu_hat_test <- matrix(NA, iter, nrow(test_X))
  varcounts   <- matrix(NA, nrow = iter, ncol = ncol(X))
  CP          <- matrix(NA, nrow = iter, ncol = n_train)

  forest <- SoftBart::MakeForest(hypers, opts)

  # Prior for lambda: lambda ~ N(0, d^2); large d = weakly informative
  d      <- 1000
  Z      <- rep(1, n_train)
  Lambda <- 1

  scale_Y  <- sd(Y)
  center_Y <- mean(Y)
  Y        <- (Y - center_Y) / scale_Y

  XXtest <- quantile_normalize(X, test_X)
  X      <- XXtest$X
  test_X <- XXtest$test_X

  s1      <- 0.1
  index_w <- cumsum(c(1, ni))[-(length(ni) + 1)]
  wu      <- w[index_w]

  ru <- rep(0, length(wu))
  for (j in seq_along(wu)) ru[j] <- rnorm(1, 0, sd = s1 / sqrt(wu[j]))
  r      <- rep(ru, times = ni)
  rsaved <- matrix(nrow = iter, ncol = length(r))

  EST_Tau    <- NULL
  EST_Lambda <- NULL

  for (i in 1:iter) {

    # Step 1: Update SBART trees b(x)
    R                <- Y - Lambda * Z - r
    mu_hat[i, ]      <- forest$do_gibbs_weighted(X, R, w, X, 1)
    mu_hat_test[i, ] <- forest$do_predict(test_X)
    varcounts[i, ]   <- as.numeric(forest$get_counts())
    Tau              <- forest$get_sigma()^2
    delta            <- Y - mu_hat[i, ] - r

    # Step 2: Update Z_ij+ (half-normal latent variables)
    Z <- truncnorm::rtruncnorm(n    = n_train,
                               a    = 0,
                               b    = Inf,
                               mean = delta * Lambda / (Lambda^2 + Tau),
                               sd   = sqrt((Tau / w) / (Lambda^2 + Tau)))

    # Step 3: Update lambda
    Lambda <- rnorm(1,
                    mean = d * sum(Z * delta * w) / (d * sum(w * Z^2) + Tau),
                    sd   = sqrt((Tau * d) / (d * sum(w * Z^2) + Tau)))

    # Step 4: Update subject random effects r_i
    mv    <- delta - Lambda * Z
    q1    <- sum_subgroups(mv, ni)
    denom <- ni * (s1^2) + Tau
    for (j in seq_along(wu)) {
      ru[j] <- rnorm(1,
                     mean = q1[j] * s1^2 / denom[j],
                     sd   = sqrt(Tau * s1^2 / (denom[j] * wu[j])))
    }

    # Step 5: Update sigma1
    s1 <- sqrt(1 / rgamma(1,
                           shape = 1 + length(ni) / 2,
                           rate  = 0.1 + 0.5 * sum(ru^2 * wu)))

    r           <- rep(ru, times = ni)
    rsaved[i, ] <- r
    CP[i, ]     <- dnorm(Y,
                         mean = mu_hat[i, ] + (Lambda / sqrt(w)) * sqrt(2 / pi),
                         sd   = sqrt(Tau / w))
    EST_Tau    <- c(EST_Tau,    Tau)
    EST_Lambda <- c(EST_Lambda, Lambda)

    if (i %% 100 == 0) cat("\rFinishing iteration", i, "of", iter, "\t")
  }

  f_hat_train      <- mu_hat[-c(1:burn), ]      * scale_Y + center_Y +
                      rsaved[-c(1:burn), ]       * scale_Y
  f_hat_test       <- mu_hat_test[-c(1:burn), ]  * scale_Y + center_Y
  f_hat_train_mean <- colMeans(f_hat_train)
  f_hat_test_mean  <- colMeans(f_hat_test)

  lambda <- EST_Lambda[-c(1:burn)] * scale_Y
  tau    <- sqrt(EST_Tau[-c(1:burn)] * scale_Y^2)
  alpha  <- lambda / tau
  sigma  <- tau * sqrt(1 + alpha^2)

  y_hat_train <- f_hat_train +
    matrix(lambda * sqrt(2 / pi), ncol = 1) %*% (1 / matrix(w, nrow = 1))
  y_hat_test  <- f_hat_test +
    matrix(lambda * sqrt(2 / pi), ncol = 1) %*%
    matrix(1, nrow = 1, ncol = nrow(test_X))

  y_hat_train_mean <- colMeans(y_hat_train)
  y_hat_test_mean  <- colMeans(y_hat_test)

  Y <- Y * scale_Y + center_Y
  like_iter      <- function(t) {
    sn::dsn(Y, xi = f_hat_train[t, ], omega = sigma[t] / sqrt(w),
            alpha = alpha[t], log = TRUE)
  }
  likelihood_mat <- t(sapply(1:length(alpha), like_iter))

  return(list(
    y_hat_train      = y_hat_train,      y_hat_test       = y_hat_test,
    y_hat_train_mean = y_hat_train_mean, y_hat_test_mean  = y_hat_test_mean,
    f_hat_train      = f_hat_train,      f_hat_test       = f_hat_test,
    f_hat_train_mean = f_hat_train_mean, f_hat_test_mean  = f_hat_test_mean,
    sigma            = sigma,            tau              = tau,
    alpha            = alpha,            lambda           = lambda,
    varcounts        = varcounts,        CP               = CP,
    likelihood_mat   = likelihood_mat,   r                = r,
    rsaved           = rsaved
  ))
}


#' SBART with Skew-t Error for Complex Survey Data
#'
#' Extends \code{SNSBARTW} by introducing a Gamma mixing variable to produce
#' heavier-tailed Skew-t errors. As the degrees of freedom nu -> Inf the model
#' reduces to \code{SNSBARTW}. See Section 2.4 of Mandal et al. (2025).
#'
#' The hierarchical model is:
#' \deqn{y_{ij} | b, r_i, Z_{ij}^+, G_{ij} \sim
#'       N(b(x_i) + r_i + \lambda Z_{ij}^+ / \sqrt{G_{ij}},
#'         \gamma^2 / (G_{ij} w_i))}
#' \deqn{Z_{ij}^+ \sim N^+(0, 1/w_i), \quad G_{ij} \sim \text{Gamma}(\nu/2, \nu/2)}
#' \deqn{r_i \sim N(0, \sigma_1^2/w_i), \quad \nu \sim \text{Gamma}(a_2, c_2)}
#'
#' @inheritParams SNSBARTW
#'
#' @return A named list containing the same elements as \code{SNSBARTW},
#'   plus \code{nu} (posterior samples of degrees of freedom).
#'
#' @references
#' Mandal A, Linero AR, Bandyopadhyay D, Sinha D (2025).
#' Soft Bayesian Additive Regression Trees (SBART) for correlated survey
#' response with non-Gaussian error.
#' \emph{Journal of Nonparametric Statistics}.
#' \doi{10.1080/10485252.2025.2574708}
#'
#' @examples
#' \dontrun{
#' library(SoftBart)
#' dat    <- sim_fried_ST(N = 100, P = 5, alpha = 5, sigma = 1,
#'                         sigma1 = 0.5, nu = 6)
#' hypers <- SoftBart::Hypers(dat$X, dat$Y, num_tree = 50)
#' opts   <- SoftBart::Opts(num_burn = 500, num_save = 500)
#' fit    <- TSBARTW(X = dat$X, Y = dat$Y, test_X = dat$X,
#'                   w = dat$w, ni = dat$ni,
#'                   hypers = hypers, opts = opts)
#' sqrt(mean(dat$w * (dat$Y - fit$y_hat_train_mean)^2))
#' }
#'
#' @export
TSBARTW <- function(X, Y, test_X, hypers = NULL, opts = NULL, w, ni) {

  if (is.null(opts))   opts   <- SoftBart::Opts()
  if (is.null(hypers)) hypers <- SoftBart::Hypers(X, Y)

  iter  <- opts$num_burn + opts$num_save
  burn  <- opts$num_burn
  opts  <- SoftBart::Opts()

  n_train     <- nrow(X)
  mu_hat      <- matrix(NA, iter, n_train)
  mu_hat_test <- matrix(NA, iter, nrow(test_X))
  forest      <- SoftBart::MakeForest(hypers, opts)

  d      <- 1000
  Z      <- rep(1, n_train)
  G      <- rep(1, n_train)
  Lambda <- 0.98
  nu     <- 6

  scale_Y  <- sd(Y)
  center_Y <- mean(Y)
  Y        <- (Y - center_Y) / scale_Y

  XXtest <- quantile_normalize(X, test_X)
  X      <- XXtest$X
  test_X <- XXtest$test_X

  s1      <- 0.1
  s1_a    <- 3
  s1_b    <- 0.1
  index_w <- cumsum(c(1, ni))[-(length(ni) + 1)]
  wu      <- w[index_w]

  ru <- rep(0, length(wu))
  for (j in seq_along(wu)) ru[j] <- rnorm(1, 0, sd = s1 / sqrt(wu[j]))
  r      <- rep(ru, times = ni)
  rsaved <- matrix(nrow = iter, ncol = length(r))

  EST_Tau <- EST_Lambda <- EST_nu <- NULL

  for (i in 1:iter) {

    # Step 1: Update SBART trees b(x)
    R                <- Y - (Lambda * Z) / sqrt(G) - r
    mu_hat[i, ]      <- forest$do_gibbs_weighted(X, R, w * G, X, 1)
    mu_hat_test[i, ] <- forest$do_predict(test_X)
    Tau              <- forest$get_sigma()^2
    delta            <- Y - mu_hat[i, ] - r

    # Step 2: Update Z_ij+
    for (j in 1:n_train) {
      Z[j] <- truncnorm::rtruncnorm(
        1, a = 0, b = Inf,
        mean = delta[j] * Lambda * sqrt(G[j]) / (Lambda^2 + Tau),
        sd   = sqrt((Tau / w[j]) / (Lambda^2 + Tau)))
    }

    # Step 3: Update lambda
    Lambda <- rnorm(1,
                    mean = d * sum(Z * delta * w * sqrt(G)) /
                             (d * sum(w * Z^2) + Tau),
                    sd   = sqrt((Tau * d) / (d * sum(w * Z^2) + Tau)))

    # Step 4: Update subject random effects r_i
    tmp1 <- G * w * (Y - mu_hat[i, ] - Lambda * Z / sqrt(G))
    q1   <- sum_subgroups(tmp1, ni)
    tmp2 <- G * w * (s1^2)
    q2   <- sum_subgroups(tmp2, ni)
    for (j in seq_along(wu)) {
      ru[j] <- rnorm(1,
                     mean = q1[j] * s1^2 / (Tau * wu[j] + q2[j]),
                     sd   = sqrt(Tau * s1^2 / (Tau * wu[j] + q2[j])))
    }

    # Step 5: Update sigma1
    s1 <- sqrt(invgamma::rinvgamma(1,
                                    shape = s1_a + length(ni) / 2,
                                    rate  = s1_b + 0.5 * sum(ru^2 * wu)))

    r           <- rep(ru, times = ni)
    rsaved[i, ] <- r

    # Step 6: Update G_ij via Metropolis-Hastings
    for (j in 1:n_train) {
      t1 <- delta[j]; t2 <- Z[j]; t3 <- w[j]
      loglik_g <- function(gval) {
        if (gval <= 0) return(-Inf)
        log(sqrt(gval)) -
          (1 / (2 * Tau)) * (t1 - Lambda * t2 / sqrt(gval))^2 * t3 * gval
      }
      tmp <- NULL
      invisible(capture.output(
        tmp <- diversitree::mcmc(
          lik = loglik_g, nsteps = 1, w = 1, x.init = G[j],
          prior = function(x) dgamma(x, shape = nu/2, rate = nu/2, log = TRUE),
          lower = 0, upper = Inf)
      ))
      G[j] <- tmp$pars
    }

    # Step 7: Update nu via Metropolis-Hastings
    loglik_nu <- function(eta) {
      if (eta <= 2) return(-Inf)
      N <- length(G); sumG <- sum(G); suml <- sum(log(G))
      ll <- N * (eta/2) * log(eta/2) - N * lgamma(eta/2) +
            (eta/2 - 1) * suml - (eta/2) * sumG
      return(ll)
    }
    nu <- local({
      tmp <- NULL
      capture.output(
        tmp <- diversitree::mcmc(
          lik = loglik_nu, nsteps = 1, w = 1, x.init = nu,
          prior = function(nu) {
            if (nu <= 2) return(-Inf)
            dgamma(nu, shape = 6, rate = 1, log = TRUE)
          },
          lower = 2, upper = Inf)
      )
      tmp$pars
    })

    EST_nu     <- c(EST_nu,     nu)
    EST_Tau    <- c(EST_Tau,    Tau)
    EST_Lambda <- c(EST_Lambda, Lambda)

    if (i %% 100 == 0) cat("\rFinishing iteration", i, "of", iter, "\t")
  }

  f_hat_train      <- mu_hat[-c(1:burn), ]     * scale_Y + center_Y +
                      rsaved[-c(1:burn), ]      * scale_Y
  f_hat_test       <- mu_hat_test[-c(1:burn), ] * scale_Y + center_Y
  f_hat_train_mean <- colMeans(f_hat_train)
  f_hat_test_mean  <- colMeans(f_hat_test)

  lambda <- EST_Lambda[-c(1:burn)] * scale_Y
  tau    <- sqrt(EST_Tau[-c(1:burn)] * scale_Y^2)
  alpha  <- lambda / tau
  sigma  <- tau * sqrt(1 + alpha^2)
  NU     <- EST_nu[-c(1:burn)]
  Enu    <- sqrt(NU / 2) * gamma((NU - 1) / 2) / gamma(NU / 2)

  y_hat_train <- f_hat_train +
    matrix(lambda * Enu * sqrt(2 / pi), ncol = 1) %*% (1 / matrix(w, nrow = 1))
  y_hat_test  <- f_hat_test +
    matrix(lambda * Enu * sqrt(2 / pi), ncol = 1) %*%
    matrix(1, nrow = 1, ncol = nrow(test_X))

  y_hat_train_mean <- colMeans(y_hat_train)
  y_hat_test_mean  <- colMeans(y_hat_test)

  Y              <- Y * scale_Y + center_Y
  likelihood_mat <- matrix(NA, nrow = length(sigma), ncol = length(Y))
  for (t in 1:length(sigma)) {
    for (j in 1:length(Y)) {
      likelihood_mat[t, j] <- sn::dst(Y[j],
                                       xi    = f_hat_train[t, j],
                                       omega = sqrt(sigma[t]^2 / w[j]),
                                       alpha = alpha[t],
                                       nu    = NU[t],
                                       log   = TRUE)
    }
  }

  return(list(
    y_hat_train      = y_hat_train,      y_hat_test       = y_hat_test,
    y_hat_train_mean = y_hat_train_mean, y_hat_test_mean  = y_hat_test_mean,
    f_hat_train      = f_hat_train,      f_hat_test       = f_hat_test,
    f_hat_train_mean = f_hat_train_mean, f_hat_test_mean  = f_hat_test_mean,
    sigma            = sigma,            tau              = tau,
    alpha            = alpha,            lambda           = lambda,
    nu               = NU,
    likelihood_mat   = likelihood_mat,   r                = r,
    rsaved           = rsaved
  ))
}


#' SBART Quantile Regression for Complex Survey Data
#'
#' Fits a nonparametric Bayesian quantile regression at a target quantile q
#' for clustered complex survey data. Uses the Asymmetric Laplace Distribution
#' (ALD) likelihood via its location-scale mixture of normals representation
#' (Kozumi & Kobayashi 2011). See Section 2.2 of Mandal et al. (2025).
#'
#' The hierarchical model is:
#' \deqn{y_{ij} | b, r_i, E_{ij} \sim N(b(x_i) + u_1 E_{ij} +
#'       u_2 \sqrt{E_{ij}} r_i, u_2^2 E_{ij} \sigma_0^2 / w_i)}
#' \deqn{E_{ij} \sim \text{Exp}((\sigma_0^2 + \sigma_1^2)/w_i), \quad
#'       r_i \sim N(0, \sigma_1^2/w_i)}
#' where \eqn{u_1 = (1-2q)/(q(1-q))} and \eqn{u_2 = \sqrt{2/(q(1-q))}}.
#'
#' @inheritParams SNSBARTW
#' @param q numeric scalar; target quantile in (0, 1)
#'
#' @return A named list containing:
#' \describe{
#'   \item{y_hat_train}{posterior samples of b(x) for training data}
#'   \item{y_hat_test}{posterior samples of b(x) for test data}
#'   \item{y_hat_train_mean}{posterior mean of b(x), training}
#'   \item{y_hat_test_mean}{posterior mean of b(x), test}
#'   \item{f_hat_train}{same as y_hat_train}
#'   \item{f_hat_test}{same as y_hat_test}
#'   \item{f_hat_train_mean}{same as y_hat_train_mean}
#'   \item{f_hat_test_mean}{same as y_hat_test_mean}
#'   \item{sigma}{posterior samples of sigma0}
#'   \item{varcounts}{variable inclusion counts per iteration}
#'   \item{likelihood_mat}{log ALD likelihood matrix (n_save x n_obs)}
#'   \item{ru}{final MCMC draw of random effects at subject level}
#'   \item{rsaved}{all MCMC draws of r (matrix: n_iter x n_obs)}
#' }
#'
#' @references
#' Mandal A, Linero AR, Bandyopadhyay D, Sinha D (2025).
#' Soft Bayesian Additive Regression Trees (SBART) for correlated survey
#' response with non-Gaussian error.
#' \emph{Journal of Nonparametric Statistics}.
#' \doi{10.1080/10485252.2025.2574708}
#'
#' @examples
#' \dontrun{
#' library(SoftBart)
#' dat    <- sim_fried_QR(N = 100, P = 5, sigma = 1, sigma1 = 0.5, q = 0.5)
#' hypers <- SoftBart::Hypers(dat$X, dat$Y, num_tree = 50)
#' opts   <- SoftBart::Opts(num_burn = 500, num_save = 500)
#' fit    <- QRSBARTW(X = dat$X, Y = dat$Y, test_X = dat$X,
#'                    w = dat$w, ni = dat$ni, q = 0.5,
#'                    hypers = hypers, opts = opts)
#' sqrt(mean(dat$w * (dat$Y - fit$f_hat_train_mean)^2))
#' }
#'
#' @export
QRSBARTW <- function(X, Y, test_X, hypers = NULL, opts = NULL, w, ni, q) {

  if (is.null(opts))   opts   <- SoftBart::Opts()
  if (is.null(hypers)) hypers <- SoftBart::Hypers(X, Y)

  iter <- opts$num_burn + opts$num_save
  burn <- opts$num_burn
  opts <- SoftBart::Opts()

  n_train     <- nrow(X)
  mu_hat      <- matrix(NA, iter, n_train)
  mu_hat_test <- matrix(NA, iter, nrow(test_X))
  varcounts   <- matrix(NA, nrow = iter, ncol = ncol(X))
  forest      <- SoftBart::MakeForest(hypers, opts)

  # ALD mixture constants
  U1 <- (1 - 2 * q) / (q - q^2)
  U2 <- sqrt(2 / (q - q^2))

  E <- rep(1, n_train)
  for (j in 1:n_train) E[j] <- rexp(1, 1 / w[j])

  scale_Y  <- sd(Y)
  center_Y <- mean(Y)
  Y        <- (Y - center_Y) / scale_Y

  XXtest <- quantile_normalize(X, test_X)
  X      <- XXtest$X
  test_X <- XXtest$test_X

  s1      <- 0.1
  index_w <- cumsum(c(1, ni))[-(length(ni) + 1)]
  wu      <- w[index_w]

  ru <- rep(0, length(wu))
  for (j in seq_along(wu)) ru[j] <- rnorm(1, 0, sd = s1 / sqrt(wu[j]))
  r        <- rep(ru, times = ni)
  rsaved   <- matrix(nrow = iter, ncol = length(r))
  s1_saved <- numeric(iter)
  EST_Tau  <- NULL

  for (i in 1:iter) {

    # Step 1: Update SBART trees b(x)
    R                <- Y - r * U2 * sqrt(E) - U1 * E
    mu_hat[i, ]      <- forest$do_gibbs_weighted(X, R, w / (E * U2^2), X, 1)
    mu_hat_test[i, ] <- forest$do_predict(test_X)
    varcounts[i, ]   <- as.numeric(forest$get_counts())

    # Step 2: Update sigma0^2 via Inverse-Gamma full conditional
    tb  <- 0.001 + sum((Y - r * U2 * sqrt(E) - mu_hat[i, ] - U1 * E)^2 *
                         w / (2 * U2^2 * E))
    ta  <- 0.001 + length(Y) / 2
    tau <- invgamma::rinvgamma(1, ta, tb)
    forest$set_sigma(sqrt(tau))

    # Step 3: Update E_ij via Metropolis-Hastings
    tmp_E        <- E
    log_cond_Eij <- function(E_val, yij, mu_ij, r_i, w_i, sigma0_sq, u1, u2) {
      if (E_val <= 0) return(-Inf)
      delta <- yij - mu_ij
      term1 <-  delta^2 / (u2^2 * E_val)
      term2 <- -2 * delta * r_i / (u2 * sqrt(E_val))
      term3 <-  2 * u1 * r_i / u2 * sqrt(E_val)
      term4 <-  (u1^2) / u2^2 * E_val
      -0.5 * log(E_val) - (w_i / (2 * sigma0_sq)) * (term1 + term2 + term3 + term4)
    }

    for (j in seq_along(Y)) {
      loglik_fun <- function(E_val) {
        log_cond_Eij(E_val, yij = Y[j], mu_ij = mu_hat[i, j], r_i = r[j],
                     w_i = w[j], sigma0_sq = tau, u1 = U1, u2 = U2)
      }
      tmp <- NULL
      capture.output(
        tmp <- diversitree::mcmc(
          lik = loglik_fun, nsteps = 1, w = 1, x.init = tmp_E[j],
          prior = function(E_val) dexp(E_val, rate = w[j] / (tau + s1^2),
                                       log = TRUE),
          lower = 1e-8, upper = Inf)
      )
      E[j] <- tmp$pars
    }

    # Step 4: Update subject random effects r_i
    num1 <- ((Y - mu_hat[i, ] - U1 * E) * s1^2) / (U2^2 * sqrt(E))
    nu1  <- sum_subgroups(num1, ni)
    for (j in seq_along(wu)) {
      ru[j] <- rnorm(1,
                     mean = nu1[j] / (ni[j] * s1^2 + tau),
                     sd   = sqrt((tau * s1^2) / ((ni[j] * s1^2 + tau) * wu[j])))
    }

    # Step 5: Update sigma1^2 via Metropolis-Hastings
    R_sum <- sum(ru^2 * wu)
    S_sum <- sum(E * w)
    M     <- sum(ni)
    n_sub <- length(ni)

    logpost_s1sq <- function(s1sq) {
      if (s1sq <= 0) return(-Inf)
      -(0.001 + n_sub / 2 + 1) * log(s1sq) -
        M * log(tau + s1sq) -
        (0.001 + R_sum / 2) / s1sq -
        S_sum / (tau + s1sq)
    }

    s1sq_new <- local({
      out <- NULL
      invisible(capture.output(
        out <- diversitree::mcmc(
          lik = logpost_s1sq, nsteps = 1, w = 0.1, x.init = s1^2,
          prior = function(x) 0, lower = 1e-8, upper = Inf)
      ))
      out$pars
    })

    s1          <- sqrt(s1sq_new)
    s1_saved[i] <- s1
    r           <- rep(ru, times = ni)
    rsaved[i, ] <- r
    EST_Tau     <- c(EST_Tau, tau)

    if (i %% 100 == 0) cat("\rFinishing iteration", i, "of", iter, "\t")
  }

  f_hat_train      <- mu_hat[-c(1:burn), ]     * scale_Y + center_Y
  f_hat_test       <- mu_hat_test[-c(1:burn), ] * scale_Y + center_Y
  f_hat_train_mean <- colMeans(f_hat_train)
  f_hat_test_mean  <- colMeans(f_hat_test)

  sig  <- sqrt(EST_Tau[-c(1:burn)] * scale_Y^2)
  sig1 <- s1_saved[-c(1:burn)] * scale_Y

  Y              <- Y * scale_Y + center_Y
  likelihood_mat <- matrix(NA, nrow = length(sig), ncol = length(Y))
  for (t in 1:length(sig)) {
    for (j in 1:length(Y)) {
      # gamma^2 = sigma0^2 + sigma1^2 is the total ALD scale
      likelihood_mat[t, j] <- log(ald::dALD(Y[j],
                                              mu    = f_hat_train[t, j],
                                              sigma = (sig[t]^2 + sig1[t]^2) / w[j],
                                              p     = q))
    }
  }

  return(list(
    y_hat_train      = f_hat_train,      y_hat_test       = f_hat_test,
    y_hat_train_mean = f_hat_train_mean, y_hat_test_mean  = f_hat_test_mean,
    f_hat_train      = f_hat_train,      f_hat_test       = f_hat_test,
    f_hat_train_mean = f_hat_train_mean, f_hat_test_mean  = f_hat_test_mean,
    sigma            = sig,
    varcounts        = varcounts,
    likelihood_mat   = likelihood_mat,
    ru               = ru,
    rsaved           = rsaved
  ))
}


# ==============================================================================
# 3. Simulation Functions
# ==============================================================================

#' Simulate clustered survey data with Skew-Normal errors
#'
#' Generates synthetic clustered survey data from the Skew-Normal model
#' using the Friedman (1991) regression function. Used for simulation
#' studies in Mandal et al. (2025).
#'
#' @param N integer; number of subjects
#' @param P integer; number of predictors
#' @param alpha numeric; skewness parameter
#' @param sigma numeric; total scale of Skew-Normal error
#' @param sigma1 numeric; standard deviation of random intercept
#'
#' @return A list with components X, Y, w, ni, mu, EY, EYr, Z, tau, lambda
#'
#' @export
sim_fried_SN <- function(N, P, alpha, sigma, sigma1) {
  lambda <- alpha * sigma / sqrt(1 + alpha^2)
  tau    <- sigma / sqrt(1 + alpha^2)
  ni     <- extraDistr::rdunif(N, 3, 6)
  X      <- matrix(runif(N * P), nrow = N)
  X      <- create_X(X, ni)
  mu     <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
             10 * X[, 4] + 5 * X[, 5]
  w1     <- runif(N, 1, 10)
  w      <- w1 / mean(w1)
  r      <- rnorm(N, 0, sd = sigma1 / sqrt(w))
  r      <- rep(r, times = ni)
  w      <- rep(w, times = ni)
  Z      <- abs(rnorm(sum(ni), mean = 0, sd = 1 / sqrt(w)))
  Y      <- mu + r + lambda * Z + rnorm(sum(ni), mean = 0, sd = tau / sqrt(w))
  EY     <- mu + lambda * sqrt(2 / pi) * sqrt(1 / w)
  EYr    <- mu + r + lambda * sqrt(2 / pi) * sqrt(1 / w)
  return(list(X = X, Y = Y, w = w, ni = ni, mu = mu,
              EY = EY, EYr = EYr, Z = Z, tau = tau, lambda = lambda))
}


#' Simulate clustered survey data with Skew-t errors
#'
#' Generates synthetic clustered survey data from the Skew-t model
#' using the Friedman (1991) regression function. Used for simulation
#' studies in Mandal et al. (2025).
#'
#' @param N integer; number of subjects
#' @param P integer; number of predictors
#' @param alpha numeric; skewness parameter
#' @param sigma numeric; total scale of Skew-t error
#' @param sigma1 numeric; standard deviation of random intercept
#' @param nu numeric; degrees of freedom (must be > 2)
#'
#' @return A list with components X, Y, w, ni, mu, EY, EYr
#'
#' @export
sim_fried_ST <- function(N, P, alpha, sigma, sigma1, nu) {
  lambda <- alpha * sigma / sqrt(1 + alpha^2)
  tau    <- sigma / sqrt(1 + alpha^2)
  ni     <- extraDistr::rdunif(N, 3, 6)
  X      <- matrix(runif(N * P), nrow = N)
  X      <- create_X(X, ni)
  mu     <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
             10 * X[, 4] + 5 * X[, 5]
  w1     <- runif(N, 1, 10)
  w      <- w1 / mean(w1)
  r      <- rnorm(N, 0, sd = sigma1 / sqrt(w))
  r      <- rep(r, times = ni)
  w      <- rep(w, times = ni)
  Z      <- abs(rnorm(sum(ni), mean = 0, sd = 1 / sqrt(w)))
  G      <- rgamma(sum(ni), shape = nu / 2, rate = nu / 2)
  Y      <- mu + r + (lambda * Z + rnorm(sum(ni), mean = 0,
                                          sd = tau / sqrt(w))) / sqrt(G)
  Enu    <- sqrt(nu / 2) * gamma((nu - 1) / 2) / gamma(nu / 2)
  EY     <- mu + lambda * sqrt(2 / (pi * w)) * Enu
  EYr    <- mu + r + lambda * sqrt(2 / (pi * w)) * Enu
  return(list(X = X, Y = Y, w = w, ni = ni, mu = mu, EY = EY, EYr = EYr))
}


#' Simulate clustered survey data from an Asymmetric Laplace Distribution
#'
#' Generates synthetic clustered survey data from the ALD model using the
#' Friedman (1991) regression function. Used for simulation studies in
#' Mandal et al. (2025).
#'
#' @param N integer; number of subjects
#' @param P integer; number of predictors
#' @param sigma numeric; symmetric error SD (sigma0)
#' @param sigma1 numeric; random effect SD
#' @param q numeric; target quantile in (0, 1)
#'
#' @return A list with components X, Y, w, ni, mu
#'
#' @export
sim_fried_QR <- function(N, P, sigma, sigma1, q) {
  ni        <- extraDistr::rdunif(N, 3, 6)
  X         <- matrix(runif(N * P), nrow = N)
  X         <- create_X(X, ni)
  mu        <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
               10 * X[, 4] + 5 * X[, 5]
  w_subject <- runif(N, 1, 5)
  w_subject <- w_subject / mean(w_subject)
  w         <- rep(w_subject, times = ni)
  r_subject <- rnorm(N, mean = 0, sd = sigma1 / sqrt(w_subject))
  r         <- rep(r_subject, times = ni)
  u1        <- (1 - 2 * q) / (q * (1 - q))
  u2        <- sqrt(2 / (q * (1 - q)))
  gamma_sq  <- (sigma1^2 + sigma^2) / w
  E         <- rexp(length(mu), rate = 1 / gamma_sq)
  Z         <- rnorm(length(mu), mean = 0, sd = sqrt(gamma_sq))
  Y         <- mu + u1 * E + u2 * sqrt(E) * (Z + r)
  return(list(X = X, Y = Y, w = w, ni = ni, mu = mu))
}
