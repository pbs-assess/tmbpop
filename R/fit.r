#' Fits POP model to data supplied as POPobj
#'
#' @param obj A list of class `obj` coming from either `build()`
#' (user-supplied data) or `sim()` (simulated data).
#' @param trace TODO
#' @param optim.control TODO
#'
#' @return A list with the following objects:
#'   * theta: vector of estimated fixed parameters
#'   * se.theta: vector of standard errors of fixed parameters
#' @export

fit <- function(obj,
  trace = TRUE,
  optim.control = list(eval.max = 5000, iter.max = 5000)) {

  # Setup, minimization and sdreport
  if (!is(obj, 'obj'))
    stop('Please supply an object of class "obj".')

  .obj <- TMB::MakeADFun(
    data = obj$datalist,
    parameters = obj$parlist,
    random = 'logRt',
    DLL = "POP",
    silent = !trace
  )

  wallclock <- proc.time()[3]
  opt <-
    nlminb(
      start = .obj$par,
      objective = .obj$fn,
      gradient = .obj$gr,
      control = optim.control
    )

  rep <- TMB::sdreport(.obj)
  summary.rep <- summary(rep) # ADREPORT stuff
  meanrep <- .obj$report(.obj$env$last.par.best) # REPORT stuff
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
    'theta' = est.theta,
    'se.theta' = se.est.theta,
    'R' = pred.Rt,
    'se.R' = se.pred.Rt,
    'B' = pred.Bt,
    'se.B' = se.pred.Bt,
    'V' = pred.Vt,
    'se.V' = se.pred.Vt,
    'u' = pred.ut,
    'se.u' = se.pred.ut,
    'uaa.female' = pred.uat1,
    'se.uaa.female' = se.pred.uat1,
    'uaa.male' = pred.uat2,
    'se.uaa.male' = se.pred.uat2,
    'select.female' = pred.sag1,
    'se.select.female' = se.pred.sag1,
    'select.male' = pred.sag2,
    'se.select.male' = se.pred.sag2,
    'abund.female' = pred.Nat1,
    'se.abund.female' = se.pred.Nat1,
    'abund.male' = pred.Nat2,
    'se.abund.male' = se.pred.Nat2,
    'mean.survey1' = meanrep$meanS1t,
    'mean.survey2' = meanrep$meanS2t,
    'mean.survey3' = meanrep$meanS3t,
    'mean.paa.catch.female' = meanrep$meanpatC1,
    'mean.paa.catch.male' = meanrep$meanpatC2,
    'mean.paa.survey1.female' = meanrep$meanpatS11,
    'mean.paa.survey1.male' = meanrep$meanpatS12,
    obj = .obj,
    opt = opt
  )
}
