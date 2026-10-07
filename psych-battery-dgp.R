#' ============================================================
#' psych-battery-dgp.R
#'
#' Shared data-generating process (DGP) and diagnostic/sweep engine for a
#' simulated battery of four correlated psychological tests, grounded in
#' Vermeiren, Vandendaele & Brysbaert (2023, Behavior Research Methods):
#'
#'   - StuVoc1  specialized-word vocabulary test, 50 items, 4AFC
#'   - StuVoc2  general-word vocabulary test (parallel form), 50 items, 4AFC
#'   - GK       general-knowledge test, 65 items, 4AFC
#'   - ART3     author recognition test, 60 real-author + 30 foil names,
#'              yes/no recognition (not 4AFC)
#'
#' Every simulated participant completes every test (fully crossed, no
#' missingness). Three things are independently controllable, matching
#' the three "dials" this battery is meant to teach with:
#'
#'   (i)   participant sample size (`n`)
#'   (ii)  measurement error, via a per-test `reliability`
#'   (iii) the *true* (construct-level, disattenuated) correlation between
#'         tests' underlying abilities, which -- together with (ii) --
#'         determines how strong a correlation the *observed* test scores
#'         will show (classical test theory attenuation: observed
#'         correlation = true correlation * sqrt(reliability_x * reliability_y))
#'
#' Design (trial/item level up, as in sampling-distributions.R, so every
#' simulated score is *structurally* bounded -- no clipping anywhere):
#'
#'   Standardized true ability (correlated across tests):
#'     z_i ~ MVN(0, Sigma_true)     (one row per person, one column per test)
#'   Per test t, calibrated location/scale (see calibrate_ability_general()):
#'     theta_ti   = theta_mu_t + theta_sigma_t * z_it
#'     p_ti       = guessing_t + (1 - guessing_t) * plogis(theta_ti)
#'   4AFC tests (StuVoc1, StuVoc2, GK):
#'     score_ti   ~ Binomial(n_items_t, p_ti)
#'   ART3 (60 author + 30 foil items, no guessing floor):
#'     hits_i             ~ Binomial(60, p_ART3,i)   -- correct "yes" to an author
#'     correct_rejects_i  ~ Binomial(30, p_ART3,i)   -- correct "no" to a foil
#'     score_ART3,i       =  hits_i + correct_rejects_i
#'     false_negative_rate_i = 1 - hits_i / 60   (said "no" to a real author)
#'     false_positive_rate_i = 1 - correct_rejects_i / 30 (said "yes" to a foil)
#'
#' `reliability` and `Sigma_true` aren't set directly on the observed-score
#' scale -- `calibrate_ability_general()` solves each test's ability
#' location/scale numerically so its *emergent* reliability matches the
#' target, and `build_sigma_true()` converts a target *observed*-score
#' correlation matrix into the *true*-score (disattenuated) correlation
#' matrix that actually needs to be fed to the multivariate draw, via
#' r_true = r_observed / sqrt(reliability_x * reliability_y). This mirrors
#' real test theory: you cannot set reliability, true correlation, and
#' observed correlation all independently -- fixing any two determines the
#' third.
#'
#' A generic sweep/diagnostic engine (`run_sweep()` / `summarize_sweep()`)
#' repeatedly draws a dataset, runs a correlation test on a chosen pair of
#' tests, and classifies the result (significant? correct sign? CI covers
#' the truth? how inflated is the estimate, among significant results?) --
#' the same four-way diagnostic quartet (power/Type I, sign consistency/
#' Type S, Type M inflation ratio, coverage) used in the project's heavier
#' Aims 1-2 simulation work, scaled down to run inline rather than as a
#' background job.
#'
#' Requires: ggplot2, MASS (called via `::` only, so it never masks
#' dplyr::select()), dplyr, tidyr, purrr.
#'
#' Quarto usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   d <- simulate_battery(n = 500, seed = 1)
#'   cor(d[c("StuVoc1", "StuVoc2", "GK", "ART3")])
#'   sweep <- run_sweep(n_grid = c(40, 100, 200), true_corr_grid = c(0, 0.2, 0.67),
#'                       reps = 200, seed_base = 42)
#'   plot_power_heatmap(summarize_sweep(sweep))
#'   ```
#' ============================================================

library(ggplot2)

# ---- test specifications (Vermeiren et al. 2023, Study 5 defaults) --------

#' Item counts, target mean, and guessing floor for each test. `mu` and
#' `n_items` are taken from Vermeiren et al.'s Study 5 (Table 6 for the
#' vocabulary/GK tests, Fig. 16 for ART3); `guessing = 0.25` reflects the
#' 4-alternative-forced-choice format of StuVoc1/StuVoc2/GK, and
#' `guessing = 0` reflects that ART3 is a yes/no recognition task, not 4AFC.
default_test_specs <- function() {
  list(
    StuVoc1 = list(n_items = 50, mu = 32.3, guessing = 0.25),
    StuVoc2 = list(n_items = 50, mu = 34.0, guessing = 0.25),
    GK      = list(n_items = 65, mu = 42.2, guessing = 0.25),
    ART3    = list(n_items = 90, mu = 57.1, guessing = 0.00, split = c(authors = 60, foils = 30))
  )
}

#' Default per-test reliability (Vermeiren et al. 2023, Study 5: alpha for
#' StuVoc1/StuVoc2/GK, alpha for ART3 -- see Table 6 and Fig. 16).
default_reliability_vec <- function() {
  c(StuVoc1 = 0.88, StuVoc2 = 0.91, GK = 0.83, ART3 = 0.90)
}

#' Default target *observed*-score correlation matrix between tests
#' (Vermeiren et al. 2023, Study 5, Tables 7-8 / Fig. 12): this is what the
#' simulated *summed scores* should correlate at when reliability is left
#' at its default -- not the latent/true-score correlation, which is
#' recovered internally via `build_sigma_true()`.
default_r_obs_matrix <- function() {
  nm <- c("StuVoc1", "StuVoc2", "GK", "ART3")
  m <- matrix(c(
    1.00, 0.82, 0.67, 0.52,
    0.82, 1.00, 0.71, 0.40,
    0.67, 0.71, 1.00, 0.31,
    0.52, 0.40, 0.31, 1.00
  ), nrow = 4, dimnames = list(nm, nm))
  m
}

# ---- calibration (generalizes sampling-distributions.R's calibrate_ability()) --

#' Expectation of `f(theta)` under theta ~ N(theta_mu, theta_sigma), by
#' numerical integration (same device as sampling-distributions.R).
expect_theta <- function(f, theta_mu, theta_sigma) {
  integrate(function(theta) f(theta) * dnorm(theta, theta_mu, theta_sigma),
            lower = theta_mu - 10 * theta_sigma,
            upper = theta_mu + 10 * theta_sigma)$value
}

#' Moments of the observed score implied by (theta_mu, theta_sigma,
#' n_items, guessing). Generalizes sampling-distributions.R's
#' implied_moments() with a guessing floor: p = guessing + (1-guessing) *
#' plogis(theta), which collapses to the original model when guessing = 0.
implied_moments_general <- function(theta_mu, theta_sigma, n_items, guessing = 0) {
  p_fun  <- function(th) guessing + (1 - guessing) * plogis(th)
  p_bar  <- expect_theta(p_fun, theta_mu, theta_sigma)
  p2_bar <- expect_theta(function(th) p_fun(th)^2, theta_mu, theta_sigma)
  pq_bar <- expect_theta(function(th) { p <- p_fun(th); p * (1 - p) }, theta_mu, theta_sigma)

  mean_score  <- n_items * p_bar
  var_true    <- n_items^2 * (p2_bar - p_bar^2)
  var_error   <- n_items * pq_bar
  reliability <- var_true / (var_true + var_error)

  list(mean_score = mean_score, var_true = var_true, var_error = var_error,
       reliability = reliability, sd_true = sqrt(var_true),
       sd_total = sqrt(var_true + var_error))
}

#' Solve for (theta_mu, theta_sigma) so the trial-level model's emergent
#' mean score matches `mu` and its emergent reliability matches
#' `reliability`, for a test with `n_items` items and a `guessing` floor.
#' Same robust nested-uniroot pattern as sampling-distributions.R's
#' calibrate_ability(): an outer search over theta_sigma, an inner search
#' over theta_mu, with the outer bracket doubled until it contains a root
#' and an informative error if the target isn't achievable.
#'
#' `stop()` above catches targets that are *structurally* unreachable
#' (the search bracket never contains a root). `tolerance` is a second,
#' softer check: `uniroot()` is only guaranteed to land within its own
#' numerical tolerance of an exact root, and `implied_moments_general()`
#' itself relies on `integrate()`, so it's possible (if rare in practice)
#' to land on a (theta_mu, theta_sigma) pair whose *achieved* reliability
#' drifts from the requested target by more than is acceptable for a
#' sweep without the achieved value being checked. `warning()` -- not
#' `stop()` -- because the calibration still *succeeded* in the sense
#' that returned, just not as tightly as requested; downstream sweep
#' code should still see a result, with the warning as the signal to
#' look closer.
calibrate_ability_general <- function(mu, reliability = 0.80, n_items = 40, guessing = 0, tolerance = 0.01) {
  floor_score <- guessing * n_items
  stopifnot(
    "`mu` must be strictly between `guessing * n_items` and `n_items`" =
      mu > floor_score && mu < n_items,
    "`reliability` must be strictly between 0 and 1" =
      reliability > 0 && reliability < 1
  )

  theta_mu_for_sigma <- function(theta_sigma) {
    uniroot(function(theta_mu) implied_moments_general(theta_mu, theta_sigma, n_items, guessing)$mean_score - mu,
            lower = -15, upper = 15, extendInt = "yes")$root
  }

  reliability_gap <- function(theta_sigma) {
    theta_mu <- theta_mu_for_sigma(theta_sigma)
    implied_moments_general(theta_mu, theta_sigma, n_items, guessing)$reliability - reliability
  }
  gap_at <- function(theta_sigma) tryCatch(reliability_gap(theta_sigma), error = function(e) NA_real_)
  lower <- 0.02
  upper <- 8
  max_upper <- 128
  g <- gap_at(upper)
  while ((is.na(g) || g < 0) && upper < max_upper) {
    upper <- upper * 2
    g <- gap_at(upper)
  }
  fail_message <- sprintf(
    paste("reliability = %.3f is not achievable at mu = %.1f with n_items = %d",
          "and guessing = %.2f (it would require an implausibly large,",
          "numerically unstable ability spread). Try a lower reliability",
          "(<= ~0.99 is usually safe), a larger n_items, or a mu further",
          "from the guessing floor."),
    reliability, mu, n_items, guessing
  )
  if (is.na(g) || g < 0) stop(fail_message)
  theta_sigma <- tryCatch(uniroot(reliability_gap, lower = lower, upper = upper)$root,
                           error = function(e) stop(fail_message))
  theta_mu <- theta_mu_for_sigma(theta_sigma)
  moments <- implied_moments_general(theta_mu, theta_sigma, n_items, guessing)

  achieved_gap <- abs(moments$reliability - reliability)
  if (achieved_gap > tolerance) {
    warning(sprintf(
      paste("calibrate_ability_general(): achieved reliability (%.4f) differs",
            "from the requested target (%.4f) by %.4f, more than `tolerance`",
            "(%.3f). The calibration still returned a result, but treat",
            "diagnostics computed at this (mu, reliability, n_items, guessing)",
            "setting with caution -- this usually happens near the extremes of",
            "what's achievable (very high reliability, or mu close to the",
            "guessing floor or ceiling)."),
      moments$reliability, reliability, achieved_gap, tolerance
    ))
  }

  list(theta_mu = theta_mu, theta_sigma = theta_sigma, moments = moments)
}

# A sweep calls calibrate_ability_general() with the same (mu, reliability,
# n_items, guessing) many times over (reliability/n don't change within a
# cell's replicates); cache by rounded arguments so repeated calls in a
# sweep are near-instant instead of re-running uniroot()/integrate().
.calibration_cache <- new.env(parent = emptyenv())
get_calibration <- function(mu, reliability, n_items, guessing = 0) {
  key <- sprintf("mu=%.4f|rel=%.5f|n=%d|g=%.3f", mu, reliability, n_items, guessing)
  if (is.null(.calibration_cache[[key]])) {
    .calibration_cache[[key]] <- calibrate_ability_general(mu, reliability, n_items, guessing)
  }
  .calibration_cache[[key]]
}

# ---- true-score correlation matrix (attenuation) ---------------------------

#' Disattenuated (true-score) correlation implied by a target *observed*
#' correlation and the two tests' reliabilities: r_true = r_obs /
#' sqrt(rel_x * rel_y). This is the classical test theory "correction for
#' attenuation" formula, applied in reverse of its usual use (normally used
#' to estimate a disattenuated correlation from data; here used to find the
#' latent correlation that will *produce* a target observed correlation
#' once measurement error is added back in).
disattenuate <- function(r_obs, reliability_x, reliability_y) {
  denom <- sqrt(reliability_x * reliability_y)
  r_true <- r_obs / denom
  if (abs(r_true) > 1) {
    warning(sprintf(
      paste("Target observed correlation r = %.2f is not achievable at",
            "reliabilities %.2f and %.2f (it would require a true-score",
            "correlation above 1); capping the true-score correlation at",
            "%.3f."),
      r_obs, reliability_x, reliability_y, sign(r_true) * 0.999
    ))
    r_true <- sign(r_true) * 0.999
  }
  r_true
}

#' Build the 4x4 true-score (latent) correlation matrix to feed to
#' MASS::mvrnorm(), from a target *observed*-score correlation matrix and a
#' named vector of per-test reliabilities. If the resulting matrix isn't
#' positive-definite (possible after capping individual cells, or after
#' overriding a single pair of tests in a sweep), it's nudged to the
#' nearest valid correlation matrix by eigenvalue clipping, with a warning.
build_sigma_true <- function(r_obs_matrix, reliability) {
  nm <- colnames(r_obs_matrix)
  stopifnot(
    "`r_obs_matrix` must be named on both dimensions" = !is.null(nm),
    "`reliability` must have an entry for every test in `r_obs_matrix`" =
      all(nm %in% names(reliability))
  )
  Sigma <- diag(length(nm))
  dimnames(Sigma) <- list(nm, nm)
  for (i in seq_along(nm)) {
    for (j in seq_along(nm)) {
      if (i == j) next
      Sigma[i, j] <- disattenuate(r_obs_matrix[nm[i], nm[j]], reliability[[nm[i]]], reliability[[nm[j]]])
    }
  }
  ev <- eigen(Sigma, symmetric = TRUE)
  if (any(ev$values < 1e-8)) {
    warning("Sigma_true was not positive-definite after disattenuation; nudging to the nearest valid correlation matrix (eigenvalue clipping).")
    vals <- pmax(ev$values, 1e-6)
    adj <- ev$vectors %*% diag(vals) %*% t(ev$vectors)
    d <- sqrt(diag(adj))
    Sigma <- adj / outer(d, d)
    dimnames(Sigma) <- list(nm, nm)
  }
  Sigma
}

# ---- simulation -------------------------------------------------------------

#' Simulate one study: `n` participants' scores on all four tests.
#'
#' @param n              sample size.
#' @param reliability    named list/vector overriding one or more tests'
#'                        reliability (default: `default_reliability_vec()`).
#'                        Dial (ii) -- measurement error.
#' @param Sigma_true      a 4x4 named true-score correlation matrix (as
#'                        produced by `build_sigma_true()`). If NULL
#'                        (default), it's built from
#'                        `default_r_obs_matrix()` and whatever
#'                        `reliability` is in effect, so the *observed*
#'                        scores land near the Vermeiren-default
#'                        correlations unless overridden. Dial (iii) --
#'                        true/construct-level correlation.
#' @param seed           for reproducibility.
#' @param keep_latent    if TRUE, also return each test's `_theta`
#'                        (ability, logit scale) and `_p` (per-item
#'                        probability) columns -- useful for directly
#'                        showing measurement error (true theta vs.
#'                        observed score).
#' @return a data frame: `id`, one score column per test (`StuVoc1`,
#'   `StuVoc2`, `GK`, `ART3`), plus `ART3_fp_rate`/`ART3_fn_rate` (per-
#'   person false-positive/false-negative rates on the recognition task).
simulate_battery <- function(n, reliability = NULL, Sigma_true = NULL, seed = NULL, keep_latent = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  specs <- default_test_specs()
  rel <- as.list(default_reliability_vec())
  if (!is.null(reliability)) {
    stopifnot("`reliability` names must be a subset of the four test names" =
                all(names(reliability) %in% names(rel)))
    rel[names(reliability)] <- as.list(reliability)
  }
  if (is.null(Sigma_true)) {
    Sigma_true <- build_sigma_true(default_r_obs_matrix(), rel)
  }
  stopifnot(
    "`Sigma_true` must be a 4x4 matrix named StuVoc1/StuVoc2/GK/ART3" =
      all(dim(Sigma_true) == 4) && setequal(colnames(Sigma_true), names(specs)),
    "`Sigma_true` must have a unit diagonal (it's a correlation matrix)" =
      isTRUE(all.equal(unname(diag(Sigma_true)), rep(1, 4), tolerance = 1e-6))
  )
  Sigma_true <- Sigma_true[names(specs), names(specs)]

  z <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = Sigma_true)
  colnames(z) <- names(specs)

  out <- data.frame(id = seq_len(n))
  for (t in names(specs)) {
    sp  <- specs[[t]]
    cal <- get_calibration(sp$mu, rel[[t]], sp$n_items, sp$guessing)
    theta <- cal$theta_mu + cal$theta_sigma * z[, t]
    p <- sp$guessing + (1 - sp$guessing) * plogis(theta)

    if (identical(t, "ART3")) {
      hits    <- rbinom(n, sp$split[["authors"]], p)
      rejects <- rbinom(n, sp$split[["foils"]], p)
      out$ART3          <- hits + rejects
      out$ART3_fn_rate  <- 1 - hits / sp$split[["authors"]]
      out$ART3_fp_rate  <- 1 - rejects / sp$split[["foils"]]
    } else {
      out[[t]] <- rbinom(n, sp$n_items, p)
    }
    if (keep_latent) {
      out[[paste0(t, "_theta")]] <- theta
      out[[paste0(t, "_p")]]     <- p
    }
  }
  out
}

#' Simulate `n` participants on just *two* named tests, with a directly
#' specified true (observed-score-scale) correlation and reliabilities.
#' This is what the NHST/power/error-rate sweep engine below uses, rather
#' than `simulate_battery()`'s full four-test matrix: powering/diagnosing a
#' test only ever concerns one focal pair, and restricting the draw to
#' that pair avoids needing the *other* two tests' correlations with the
#' focal pair to stay mutually consistent (a full 4x4 correlation matrix
#' stops being positive-definite surprisingly easily once one cell is
#' pushed to an extreme value for a sweep, as the other cells were not
#' built with that in mind).
#'
#' @param n           sample size.
#' @param var_x, var_y  which two tests to simulate (must be among the
#'                      four names in `default_test_specs()`).
#' @param true_corr    target correlation between the two tests, on
#'                      whichever scale `corr_scale` selects.
#' @param corr_scale   `"observed"` (default): `true_corr` is the target
#'                      *observed*-score correlation, and the function
#'                      disattenuates internally to find the latent
#'                      correlation that will produce it once measurement
#'                      error is added -- use this for power/NHST work,
#'                      where the population parameter being tested is the
#'                      observed-score correlation. `"latent"`: `true_corr`
#'                      is used directly as the correlation between the two
#'                      tests' *true* abilities, with no compensation -- use
#'                      this to demonstrate attenuation itself (holding the
#'                      construct-level relationship fixed while varying
#'                      reliability, so the *observed* correlation visibly
#'                      drops as reliability drops).
#' @param reliability_x, reliability_y  reliability for each test
#'                      (defaults to each test's Vermeiren-default value).
#' @param seed         for reproducibility.
#' @return a data frame: `id`, and one score column per requested test
#'   (plus `ART3_fp_rate`/`ART3_fn_rate` if `var_x`/`var_y` includes ART3).
simulate_pair <- function(n, var_x, var_y, true_corr = 0, corr_scale = c("observed", "latent"),
                           reliability_x = NULL, reliability_y = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  corr_scale <- match.arg(corr_scale)
  specs <- default_test_specs()
  stopifnot(
    "`var_x`/`var_y` must each be one of StuVoc1, StuVoc2, GK, ART3" =
      var_x %in% names(specs) && var_y %in% names(specs),
    "`var_x` and `var_y` must be different tests" = var_x != var_y
  )
  rel_defaults <- default_reliability_vec()
  if (is.null(reliability_x)) reliability_x <- unname(rel_defaults[[var_x]])
  if (is.null(reliability_y)) reliability_y <- unname(rel_defaults[[var_y]])

  r_true <- if (corr_scale == "observed") disattenuate(true_corr, reliability_x, reliability_y) else true_corr
  Sigma <- matrix(c(1, r_true, r_true, 1), nrow = 2)
  z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)

  out <- data.frame(id = seq_len(n))
  for (slot in list(list(name = var_x, z = z[, 1], rel = reliability_x),
                     list(name = var_y, z = z[, 2], rel = reliability_y))) {
    sp  <- specs[[slot$name]]
    cal <- get_calibration(sp$mu, slot$rel, sp$n_items, sp$guessing)
    theta <- cal$theta_mu + cal$theta_sigma * slot$z
    p <- sp$guessing + (1 - sp$guessing) * plogis(theta)
    if (identical(slot$name, "ART3")) {
      hits    <- rbinom(n, sp$split[["authors"]], p)
      rejects <- rbinom(n, sp$split[["foils"]], p)
      out$ART3         <- hits + rejects
      out$ART3_fn_rate <- 1 - hits / sp$split[["authors"]]
      out$ART3_fp_rate <- 1 - rejects / sp$split[["foils"]]
    } else {
      out[[slot$name]] <- rbinom(n, sp$n_items, p)
    }
  }
  out
}

#' Simulate `n_studies` independent studies of size `n`, for sampling-
#' distribution work (distribution of a sample statistic across repeated
#' studies). Returns one long data frame with a `study` id column, in the
#' same spirit as sampling-distributions.R's simulate_studies().
simulate_studies_battery <- function(n, n_studies, ..., seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- lapply(seq_len(n_studies), function(s) {
    d <- simulate_battery(n, ..., seed = NULL)
    d$study <- paste("Study", s)
    d
  })
  do.call(rbind, out)
}

# ---- NHST / diagnostic / sweep engine ---------------------------------------

#' Run a two-sided Pearson correlation test between two tests' scores.
fit_focal_cor <- function(data, var_x, var_y) {
  ct <- stats::cor.test(data[[var_x]], data[[var_y]])
  data.frame(estimate = unname(ct$estimate), ci_low = ct$conf.int[1],
             ci_high = ct$conf.int[2], p_value = ct$p.value)
}

#' Classify one fitted test's result against the *true* value used to
#' simulate the data: significant? (CI/p-value excludes zero at `alpha`),
#' sign-correct? (only meaningful when true_value != 0 -- NA under the
#' null, since there's no "correct" nonzero sign to match), CI coverage?,
#' and `inflation_ratio` (the Type M/magnitude-error ratio,
#' `|estimate| / |true_value|` -- also NA under the null, for the same
#' reason sign-correctness is: "exaggeration relative to zero" isn't a
#' meaningful quantity, since any nonzero estimate divided by zero is
#' undefined). This mirrors the project's heavier Aims 1-2 `glmer`
#' simulation work, which computed the same ratio (there called
#' `inflation_ratio`) from an identical "true value is known because we
#' simulated it" setup.
diagnose_cor_result <- function(fit, true_value, alpha = 0.05) {
  significant     <- fit$p_value < alpha
  sign_correct    <- if (true_value == 0) NA else (sign(fit$estimate) == sign(true_value))
  covers          <- (true_value >= fit$ci_low) & (true_value <= fit$ci_high)
  inflation_ratio <- if (true_value == 0) NA_real_ else abs(fit$estimate / true_value)
  data.frame(estimate = fit$estimate, p_value = fit$p_value,
             significant = significant, sign_correct = sign_correct, covers = covers,
             inflation_ratio = inflation_ratio)
}

#' Sweep sample size, the focal pair's true (population) observed-score
#' correlation, and (optionally) both focal tests' shared reliability,
#' repeating `reps` independent simulated studies per grid cell. Mirrors
#' the project's Aims 1-2 simulation pattern (deterministic per-cell
#' seeding, one row per replicate) scaled down to run inline.
#'
#' @param n_grid           sample sizes to sweep.
#' @param true_corr_grid    population *observed*-score correlations
#'                          between `var_x` and `var_y` to sweep (include 0
#'                          to get the empirical false-positive rate).
#' @param reliability_grid  reliability values to sweep, applied to *both*
#'                          `var_x` and `var_y` (other tests keep their
#'                          default reliability); NULL (default) holds
#'                          reliability fixed at the Vermeiren-default
#'                          value for `var_x`.
#' @param reps             replicate studies per grid cell.
#' @param var_x, var_y     the focal pair of tests to test/diagnose.
#' @param alpha            significance threshold.
#' @param seed_base         base seed; each (n, true_corr, reliability, rep)
#'                          cell gets its own deterministic seed derived
#'                          from it, so the whole sweep is reproducible and
#'                          any single cell can be reproduced in isolation.
#' @return a tidy data frame, one row per replicate, with columns `n`,
#'   `true_corr`, `reliability`, `rep`, `estimate`, `p_value`,
#'   `significant`, `sign_correct`, `covers`.
run_sweep <- function(n_grid, true_corr_grid, reliability_grid = NULL, reps = 200,
                       var_x = "StuVoc1", var_y = "GK", alpha = 0.05, seed_base = 1234) {
  rel_base <- default_reliability_vec()
  if (is.null(reliability_grid)) reliability_grid <- unname(rel_base[[var_x]])

  grid <- expand.grid(n = n_grid, true_corr = true_corr_grid,
                       reliability = reliability_grid, rep = seq_len(reps))

  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    n <- grid$n[i]; true_corr <- grid$true_corr[i]
    reliability <- grid$reliability[i]; rep_i <- grid$rep[i]

    seed_cell <- seed_base + rep_i * 1e6 + n * 1e3 +
      round(reliability * 100) + round((true_corr + 1) * 1000)

    d <- simulate_pair(n, var_x, var_y, true_corr = true_corr,
                        reliability_x = reliability, reliability_y = reliability,
                        seed = seed_cell)
    fit <- fit_focal_cor(d, var_x, var_y)
    diag <- diagnose_cor_result(fit, true_corr, alpha)
    rows[[i]] <- cbind(data.frame(n = n, true_corr = true_corr, reliability = reliability, rep = rep_i), diag)
  }
  do.call(rbind, rows)
}

#' Summarize `run_sweep()` output into per-cell metrics: `power_or_type1`
#' (the significance rate -- power when `true_corr != 0`, the empirical
#' Type I/false-positive rate when `true_corr == 0`), `sign_consistency`
#' (proportion of *significant* results with the correct sign -- NA for
#' `true_corr == 0` cells, where there's no correct sign to match),
#' `type_m_ratio` (median `|estimate| / |true_corr|` among *significant*
#' results only -- NA for `true_corr == 0` cells, for the same reason
#' `sign_consistency` is -- the Type M/magnitude-error "exaggeration
#' ratio": how much the significance filter inflates a reported effect,
#' reported as the median rather than the mean because significance
#' filtering disproportionately selects the occasional large outlier
#' estimate, especially at small `n`, and the median is more robust to
#' that), and `coverage` (proportion of CIs containing the true
#' correlation).
summarize_sweep <- function(results) {
  agg <- aggregate(
    cbind(significant, covers) ~ n + true_corr + reliability,
    data = results, FUN = mean
  )
  names(agg)[names(agg) == "significant"] <- "power_or_type1"
  names(agg)[names(agg) == "covers"]      <- "coverage"

  sig_df <- results[results$significant, ]

  sign_df <- sig_df[!is.na(sig_df$sign_correct), ]
  if (nrow(sign_df) > 0) {
    sign_agg <- aggregate(sign_correct ~ n + true_corr + reliability, data = sign_df, FUN = mean)
    names(sign_agg)[names(sign_agg) == "sign_correct"] <- "sign_consistency"
  } else {
    sign_agg <- unique(results[c("n", "true_corr", "reliability")])
    sign_agg$sign_consistency <- NA_real_
  }

  type_m_df <- sig_df[!is.na(sig_df$inflation_ratio), ]
  if (nrow(type_m_df) > 0) {
    type_m_agg <- aggregate(inflation_ratio ~ n + true_corr + reliability, data = type_m_df,
                             FUN = median)
    names(type_m_agg)[names(type_m_agg) == "inflation_ratio"] <- "type_m_ratio"
  } else {
    type_m_agg <- unique(results[c("n", "true_corr", "reliability")])
    type_m_agg$type_m_ratio <- NA_real_
  }

  out <- merge(agg, sign_agg, by = c("n", "true_corr", "reliability"), all.x = TRUE)
  out <- merge(out, type_m_agg, by = c("n", "true_corr", "reliability"), all.x = TRUE)
  n_reps <- aggregate(rep ~ n + true_corr + reliability, data = results, FUN = length)
  names(n_reps)[names(n_reps) == "rep"] <- "n_reps"
  out <- merge(out, n_reps, by = c("n", "true_corr", "reliability"))
  out[order(out$true_corr, out$reliability, out$n), ]
}

# ---- figures -----------------------------------------------------------------

#' Faceted line/point plot of one or more summary metrics against `x`
#' (typically sample size), one colour per value of `color_by` (typically
#' the true correlation), one facet per value of `facet_by` (typically
#' reliability, if more than one value was swept).
plot_metrics_vs_x <- function(summary_df, x = "n", metrics = "power_or_type1",
                               color_by = "true_corr", facet_by = NULL) {
  long <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(x = summary_df[[x]], value = summary_df[[m]], metric = m,
               color_by = factor(summary_df[[color_by]]),
               facet_by = if (!is.null(facet_by)) summary_df[[facet_by]] else NA)
  }))
  p <- ggplot(long, aes(x = x, y = value, color = color_by, group = color_by)) +
    geom_line() + geom_point() +
    scale_x_continuous(name = x) +
    scale_y_continuous(name = "Value", limits = c(0, 1)) +
    labs(color = color_by) +
    theme_minimal()
  if (length(metrics) > 1) p <- p + facet_grid(metric ~ ., switch = "y")
  if (!is.null(facet_by)) p <- p + facet_wrap(~facet_by, labeller = label_both)
  p
}

#' Two-dial tile heatmap of one summary metric, e.g. sample size (x) by
#' reliability (y), at a fixed true correlation -- or sample size (x) by
#' true correlation (y), at a fixed reliability. Cells are labelled with
#' the metric's value (as a percentage) for readability.
plot_power_heatmap <- function(summary_df, x = "n", y = "true_corr", fill = "power_or_type1",
                                fill_name = "Rate") {
  df <- summary_df
  df$.x <- factor(df[[x]])
  df$.y <- factor(df[[y]])
  df$.fill <- df[[fill]]
  ggplot(df, aes(x = .x, y = .y, fill = .fill)) +
    geom_tile(color = "white") +
    geom_text(aes(label = scales::percent(.fill, accuracy = 1)), color = "black", size = 3) +
    scale_fill_viridis_c(name = fill_name, limits = c(0, 1)) +
    scale_x_discrete(name = x) +
    scale_y_discrete(name = y) +
    theme_minimal()
}
