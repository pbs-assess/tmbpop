#' Fit the model to data
#'
#' @param obj A list of class `obj` coming from either [build()]
#' (user-supplied real data) or [sim()] (simulated data).
#' @param trace If `TRUE`, extra information will be printed.
#' @param control A list to pass to the `control` argument of
#'   [stats::nlminb()].
#'
#' @return A list with the following elements:
#'   * `theta`: vector of estimated fixed parameters
#'   * `se.theta`: vector of standard errors of fixed parameters
#'   * many other elements to be described later TODO
#' @export
#' @useDynLib tmbpop
#'
#' @examples
#' sim_dat <- sim(
#'   survey1 = pop_example$S1t,
#'   survey2 = pop_example$S2t,
#'   survey3 = pop_example$S3t,
#'   paa.catch.female = pop_example$patC1,
#'   paa.catch.male = pop_example$patC2,
#'   n.trips.paa.catch = pop_example$ntC,
#'   paa.survey1.female = pop_example$patS11,
#'   paa.survey1.male = pop_example$patS12,
#'   n.trips.paa.survey1 = pop_example$ntS1,
#'   catch = pop_example$Ct,
#'   paa.mature = pop_example$ma,
#'   weight.female = pop_example$wa1,
#'   weight.male = pop_example$wa2,
#'   misc.fixed.param = pop_example$MiscFixedParam,
#'   theta.ini = NULL, # if NULL, uses default values
#'   lkhd.paa = "binomial",
#'   var.paa.add = TRUE,
#'   enable.priors = TRUE
#' )
#' model <- fit(sim_dat)
#'
#' # Table of estimated fixed parameters, their standard errors and the true theta
#' table.theta <- cbind(sim_dat$true.theta, model$theta, model$se.theta)
#' dimnames(table.theta)[[2]] <- c("True", "Estimate", "s.e.")
#' round(table.theta, 3)
#'
#' # Plot predicted recruits with +- 2s.e. as envelope, with true recuits overlaid
#' lb.R <- model$R - 2 * model$se.R
#' ub.R <- model$R + 2 * model$se.R
#'
#' plot(sim_dat$years, model$R,
#'   type = "o", xlab = "years", ylab = "Recruits", pch = 19, cex = 0.5,
#'   ylim = c(2000, 7000), main = paste0("Predicted recruits, simulated data")
#' )
#' polygon(
#'   x = c(sim_dat$years, sim_dat$years[seq(sim_dat$TC, 1)]),
#'   y = c(lb.R, ub.R[seq(sim_dat$TC, 1)]), col = "#4f4f4f40", border = NA
#' )
#' lines(sim_dat$years, sim_dat$true.R, col = "red")
#'
#' # Plot predicted biomass with +- 2s.e. as envelope, with true biomass overlaid
#' lb.B <- model$B - 2 * model$se.B
#' ub.B <- model$B + 2 * model$se.B
#'
#' plot(sim_dat$years, model$B,
#'   type = "o", xlab = "years", ylab = "SSB", pch = 19, cex = 0.5,
#'   ylim = c(2000, 22000), main = paste0("Predicted biomass, simulated data")
#' )
#' polygon(
#'   x = c(sim_dat$years, sim_dat$years[seq(sim_dat$TC, 1)]),
#'   y = c(lb.B, ub.B[seq(sim_dat$TC, 1)]), col = "#4f4f4f40", border = NA
#' )
#' lines(sim_dat$years, sim_dat$true.B, col = "red")
#'
#' # Sample via NUTS MCMC with Stan
#' \donttest{
#' library("tmbstan")
#' # options(mc.cores = parallel::detectCores()) # for parallel processing
#' pop_mcmc <- tmbstan(
#'   model$obj,
#'   chains = 1, # using only 1 chain and...
#'   iter = 600, # only 600 iterations for a quick example
#'   init = list("last.par.best"),
#'   control = list(adapt_delta = 0.9, max_treedepth = 20L) # as needed, see ?stan
#' )
#'
#' pars <- c("logR0", "logM1", "logh", "logmuC", "deltaC", "logqS1") # a selection
#' bayesplot::mcmc_trace(as.array(pop_mcmc), pars = pars)
#' bayesplot::mcmc_hist(as.array(pop_mcmc), pars = pars)
#' bayesplot::mcmc_pairs(as.array(pop_mcmc), pars = pars)
#' }

fit <- function(obj,
  trace = TRUE,
  control = list(eval.max = 5000, iter.max = 5000)) {

  # Setup, minimization and sdreport
  if (!is(obj, 'obj'))
    stop('Please supply an object of class "obj".')

  .obj <- TMB::MakeADFun(
    data = obj$datalist,
    parameters = obj$parlist,
    random = 'logRt',
    DLL = "tmbpop",
    silent = !trace
  )

  wallclock <- proc.time()[3]
  opt <-
    nlminb(
      start = .obj$par,
      objective = .obj$fn,
      gradient = .obj$gr,
      control = control
    )

  rep <- TMB::sdreport(.obj)
  summary.rep <- summary(rep) # ADREPORT
  meanrep <- .obj$report(.obj$env$last.par.best) # REPORT
  elapsed.time <- proc.time()[3] - wallclock

  # Outputs
  if (trace) {
    message('Optimization and sdreport took ', round(elapsed.time, 1),
      ' seconds.')
  }

  est.theta <- summary.rep[(obj$length.theta + obj$TC + 1):(2 * obj$length.theta +
      obj$TC), 1]
  se.est.theta <- summary.rep[(obj$length.theta + obj$TC + 1):(2 *
      obj$length.theta + obj$TC), 2]

  pred.Rt <- summary.rep[dimnames(summary.rep)[[1]] == 'Rt', 1]
  se.pred.Rt <- summary.rep[dimnames(summary.rep)[[1]] == 'Rt', 2]
  pred.Bt <- summary.rep[dimnames(summary.rep)[[1]] == 'Bt', 1]
  se.pred.Bt <- summary.rep[dimnames(summary.rep)[[1]] == 'Bt', 2]
  pred.Vt <- summary.rep[dimnames(summary.rep)[[1]] == 'Vt', 1]
  se.pred.Vt <- summary.rep[dimnames(summary.rep)[[1]] == 'Vt', 2]
  pred.ut <- summary.rep[dimnames(summary.rep)[[1]] == 'ut', 1]
  se.pred.ut <- summary.rep[dimnames(summary.rep)[[1]] == 'ut', 2]

  pred.uat1 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'uat1', 1],
      obj$A, obj$TC)
  se.pred.uat1 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'uat1', 2],
      obj$A, obj$TC)
  pred.uat2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'uat2', 1],
      obj$A, obj$TC)
  se.pred.uat2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'uat2', 2],
      obj$A, obj$TC)

  pred.sag1 <- matrix(summary.rep[dimnames(summary.rep)[[1]] == 'sag1', 1],
    obj$A, 4) # G=4 surveys
  se.pred.sag1 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'sag1', 2],
      obj$A, 4) # G=4 surveys
  pred.sag2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'sag2', 1],
      obj$A, 4) # G=4 surveys
  se.pred.sag2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'sag2', 2],
      obj$A, 4) # G=4 surveys

  pred.Nat1 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'Nat1', 1],
      obj$A, obj$TC)
  se.pred.Nat1 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'Nat1', 2],
      obj$A, obj$TC)
  pred.Nat2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'Nat2', 1],
      obj$A, obj$TC)
  se.pred.Nat2 <-
    matrix(summary.rep[dimnames(summary.rep)[[1]] == 'Nat2', 2],
      obj$A, obj$TC)

  list(
    theta = est.theta,
    se.theta = se.est.theta,
    R = pred.Rt,
    se.R = se.pred.Rt,
    B = pred.Bt,
    se.B = se.pred.Bt,
    V = pred.Vt,
    se.V = se.pred.Vt,
    u = pred.ut,
    se.u = se.pred.ut,
    uaa.female = pred.uat1,
    se.uaa.female = se.pred.uat1,
    uaa.male = pred.uat2,
    se.uaa.male = se.pred.uat2,
    select.female = pred.sag1,
    se.select.female = se.pred.sag1,
    select.male = pred.sag2,
    se.select.male = se.pred.sag2,
    abund.female = pred.Nat1,
    se.abund.female = se.pred.Nat1,
    abund.male = pred.Nat2,
    se.abund.male = se.pred.Nat2,
    mean.survey1 = meanrep$meanS1t,
    mean.survey2 = meanrep$meanS2t,
    mean.survey3 = meanrep$meanS3t,
    mean.paa.catch.female = meanrep$meanpatC1,
    mean.paa.catch.male = meanrep$meanpatC2,
    mean.paa.survey1.female = meanrep$meanpatS11,
    mean.paa.survey1.male = meanrep$meanpatS12,
    obj = .obj,
    opt = opt
  )
}
