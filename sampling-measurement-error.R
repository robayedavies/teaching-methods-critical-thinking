#' ============================================================
#' sampling-measurement-error.R
#'
#' Demonstrates the two distinct "error" dials exposed by the shared
#' battery simulator in psych-battery-dgp.R (StuVoc1, StuVoc2, GK, ART3;
#' grounded in Vermeiren, Vandendaele & Brysbaert, 2023):
#'
#'   MEASUREMENT error -- noise in a single person's score, caused by a
#'     test having only finitely many items (or author/foil names) to
#'     sample responses from. Controlled by `reliability` (dial ii, from
#'     psych-battery-dgp.R). Lowering it widens the scatter of a
#'     person's *observed* score around their *true* ability, and --
#'     crucially -- *attenuates* (shrinks toward zero) any correlation
#'     computed from a noisy observed score, with no change in sample
#'     size at all.
#'
#'   SAMPLING error -- noise in a *study's* estimate (a sample mean or a
#'     sample correlation) caused by drawing only `n` participants from
#'     the population rather than all of them. Controlled by
#'     participant sample size `n` (dial i). Shrinking it narrows the
#'     *spread* of the sampling distribution of a statistic across
#'     repeated studies, without shifting where that distribution is
#'     centred.
#'
#' The explicit teaching contrast this script is built around, and the
#' reason the two kinds of error get separate figures rather than one:
#'
#'   changing N changes the SPREAD of a statistic's sampling distribution
#'   (plot_sampling_distribution());
#'
#'   changing reliability changes its CENTRE, via attenuation
#'   (plot_attenuation_vs_reliability()), not (much) its spread --
#'   these are different error sources with different fixes (a bigger N
#'   does little to undo attenuation from an unreliable measure, and a
#'   more reliable measure does little to narrow a sampling distribution
#'   that's wide because N is small).
#'
#' Four simulate_*() data generators (each a thin, teaching-oriented
#' wrapper around psych-battery-dgp.R's calibration/simulation engine)
#' pair with four plot_*() figures:
#'
#'   simulate_single_test()            -> plot_true_vs_observed()
#'   simulate_parallel_forms_battery() -> plot_parallel_forms()
#'   simulate_sampling_distribution()  -> plot_sampling_distribution()
#'   simulate_attenuation_sweep()      -> plot_attenuation_vs_reliability()
#'
#' Requires: psych-battery-dgp.R sourced first (uses default_test_specs(),
#' default_reliability_vec(), default_r_obs_matrix(), get_calibration(),
#' simulate_pair(), disattenuate()); ggplot2; ggdist (for
#' plot_sampling_distribution() only, same guarded-require pattern as
#' sampling-distributions.R's plot_dotpile()).
#'
#' Quarto usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("sampling-measurement-error.R")
#'
#'   # measurement error: one test, one reliability
#'   plot_true_vs_observed("StuVoc1", n = 300, seed = 1)
#'   plot_parallel_forms("StuVoc1", n = 300, seed = 2)
#'
#'   # sampling error: spread shrinks as N grows, centre doesn't move
#'   plot_sampling_distribution(n_grid = c(20, 80, 320), n_studies = 300,
#'                               statistic = "cor",
#'                               var_x = "StuVoc1", var_y = "GK", seed = 3)
#'
#'   # measurement error again, now via attenuation of a correlation
#'   plot_attenuation_vs_reliability(reliability_grid = seq(0.2, 0.95, by = 0.15),
#'                                    n = 300, reps = 300,
#'                                    var_x = "StuVoc1", var_y = "GK", seed = 4)
#'   ```
#' ============================================================

library(ggplot2)

# ---- measurement error: one test, one person at a time ---------------------

#' Simulate `n` participants on a single test from the battery, keeping
#' each person's standardized true ability (`true_ability`, z-scale) --
#' the single-test building block behind plot_true_vs_observed().
#'
#' @param test         one of "StuVoc1", "StuVoc2", "GK", "ART3".
#' @param n            number of simulated participants.
#' @param reliability  this test's reliability (dial ii); NULL (default)
#'                      uses its Vermeiren-default value.
#' @param seed         for reproducibility.
#' @return a data frame: `id`, `test`, `true_ability` (standardized,
#'   mean 0 / sd 1), `theta` (ability, logit scale), `p` (per-item/
#'   per-name probability implied by theta), `score`, `n_items` (test
#'   length score is bounded by), `reliability` (the value actually
#'   used).
simulate_single_test <- function(test = "StuVoc1", n = 500, reliability = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  specs <- default_test_specs()
  stopifnot(
    "`test` must be one of StuVoc1, StuVoc2, GK, ART3" = test %in% names(specs),
    "`n` must be a positive integer" = n >= 1
  )
  rel <- if (is.null(reliability)) unname(default_reliability_vec()[[test]]) else reliability
  stopifnot("`reliability` must be strictly between 0 and 1" = rel > 0 && rel < 1)

  sp  <- specs[[test]]
  cal <- get_calibration(sp$mu, rel, sp$n_items, sp$guessing)
  z     <- rnorm(n)
  theta <- cal$theta_mu + cal$theta_sigma * z
  p     <- sp$guessing + (1 - sp$guessing) * plogis(theta)

  if (identical(test, "ART3")) {
    hits     <- rbinom(n, sp$split[["authors"]], p)
    rejects  <- rbinom(n, sp$split[["foils"]], p)
    score    <- hits + rejects
    n_items  <- sum(sp$split)
  } else {
    score    <- rbinom(n, sp$n_items, p)
    n_items  <- sp$n_items
  }

  data.frame(id = seq_len(n), test = test, true_ability = z, theta = theta, p = p,
             score = score, n_items = n_items, reliability = rel)
}

#' A person's observed score plotted against their (otherwise
#' unobservable) true ability -- the single clearest picture of
#' measurement error: vertical scatter of points around the expected-
#' score curve is score variation that has nothing to do with sample
#' size, only with how few items/names the test has to work with.
#'
#' @inheritParams simulate_single_test
#' @return a ggplot object.
plot_true_vs_observed <- function(test = "StuVoc1", n = 300, reliability = NULL, seed = NULL) {
  d   <- simulate_single_test(test, n, reliability, seed)
  sp  <- default_test_specs()[[test]]
  rel <- d$reliability[1]
  cal <- get_calibration(sp$mu, rel, sp$n_items, sp$guessing)

  z_grid     <- seq(-3.2, 3.2, length.out = 200)
  theta_grid <- cal$theta_mu + cal$theta_sigma * z_grid
  p_grid     <- sp$guessing + (1 - sp$guessing) * plogis(theta_grid)
  curve_df   <- data.frame(true_ability = z_grid, expected_score = d$n_items[1] * p_grid)

  ggplot(d, aes(true_ability, score)) +
    geom_point(color = "#3B6EA5", alpha = 0.6) +
    geom_line(data = curve_df, aes(true_ability, expected_score), color = "#B23A48", linewidth = 1) +
    labs(x = "True ability (standardized)", y = sprintf("%s score", test),
         subtitle = sprintf(
           "reliability = %.2f \u2014 vertical scatter around the red curve is measurement error",
           rel)) +
    theme_minimal()
}

#' Two independent "administrations" of the same test to the same `n`
#' simulated people (sharing the same true ability, independent item
#' draws) -- the parallel-forms/test-retest check used throughout this
#' project (same device as sampling-distributions.R's
#' simulate_parallel_forms()), generalized across all four battery
#' tests.
#'
#' @inheritParams simulate_single_test
#' @return a data frame: `id`, `test`, `true_ability`, `form1`, `form2`
#'   (two independently simulated scores per person), `reliability`.
simulate_parallel_forms_battery <- function(test = "StuVoc1", n = 500, reliability = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  specs <- default_test_specs()
  stopifnot("`test` must be one of StuVoc1, StuVoc2, GK, ART3" = test %in% names(specs))
  rel <- if (is.null(reliability)) unname(default_reliability_vec()[[test]]) else reliability
  stopifnot("`reliability` must be strictly between 0 and 1" = rel > 0 && rel < 1)

  sp  <- specs[[test]]
  cal <- get_calibration(sp$mu, rel, sp$n_items, sp$guessing)
  z     <- rnorm(n)
  theta <- cal$theta_mu + cal$theta_sigma * z
  p     <- sp$guessing + (1 - sp$guessing) * plogis(theta)

  draw_once <- function() {
    if (identical(test, "ART3")) {
      rbinom(n, sp$split[["authors"]], p) + rbinom(n, sp$split[["foils"]], p)
    } else {
      rbinom(n, sp$n_items, p)
    }
  }
  data.frame(id = seq_len(n), test = test, true_ability = z,
             form1 = draw_once(), form2 = draw_once(), reliability = rel)
}

#' Scatterplot of one test's two parallel-form scores, with the
#' observed correlation annotated against its reliability target --
#' reliability *is* parallel-forms correlation in this model.
#'
#' @inheritParams simulate_single_test
#' @return a ggplot object.
plot_parallel_forms <- function(test = "StuVoc1", n = 300, reliability = NULL, seed = NULL) {
  d     <- simulate_parallel_forms_battery(test, n, reliability, seed)
  r_obs <- cor(d$form1, d$form2)

  ggplot(d, aes(form1, form2)) +
    geom_point(color = "#3B6EA5", alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "#B23A48") +
    coord_equal() +
    labs(x = "Form 1 score", y = "Form 2 score",
         subtitle = sprintf("%s \u2014 target reliability = %.2f, observed r(form1, form2) = %.2f",
                             test, d$reliability[1], r_obs)) +
    theme_minimal()
}

# ---- sampling error: repeated studies, varying N ----------------------------

#' Simulate `n_studies` independent studies of size `n`, each
#' contributing one value of a chosen summary statistic -- the building
#' block for a sampling distribution.
#'
#' @param n            participants per study (dial i).
#' @param n_studies    number of independent studies to simulate.
#' @param statistic    `"mean"` (mean score on `var_x`) or `"cor"`
#'                      (Pearson correlation between `var_x` and
#'                      `var_y`).
#' @param var_x, var_y which test(s); `var_y` is required when
#'                      `statistic = "cor"`.
#' @param reliability  reliability applied to `var_x` (and `var_y`, for
#'                      `"cor"`); NULL (default) uses Vermeiren-default
#'                      values. Held fixed across studies here -- this
#'                      function is about dial (i) only; see
#'                      simulate_attenuation_sweep() for dial (ii).
#' @param true_corr    population *observed*-score correlation to
#'                      simulate from, for `statistic = "cor"`; NULL
#'                      (default) uses the Vermeiren-default value for
#'                      the chosen pair.
#' @param seed         for reproducibility.
#' @return a data frame: `study`, `n`, `statistic`, `value`.
simulate_sampling_distribution <- function(n, n_studies, statistic = c("mean", "cor"),
                                            var_x = "StuVoc1", var_y = NULL,
                                            reliability = NULL, true_corr = NULL, seed = NULL) {
  statistic <- match.arg(statistic)
  stopifnot(
    "`var_y` must be supplied when statistic = 'cor'" = statistic != "cor" || !is.null(var_y),
    "`n_studies` must be at least 2 to show a distribution" = n_studies >= 2
  )
  if (!is.null(seed)) set.seed(seed)
  rel_default <- default_reliability_vec()
  rel_x <- if (is.null(reliability)) unname(rel_default[[var_x]]) else reliability

  if (statistic == "mean") {
    values <- vapply(seq_len(n_studies), function(s) {
      mean(simulate_single_test(var_x, n, reliability = rel_x, seed = NULL)$score)
    }, numeric(1))
  } else {
    rel_y <- if (is.null(reliability)) unname(rel_default[[var_y]]) else reliability
    tc    <- if (is.null(true_corr)) default_r_obs_matrix()[var_x, var_y] else true_corr
    values <- vapply(seq_len(n_studies), function(s) {
      d <- simulate_pair(n, var_x, var_y, true_corr = tc, corr_scale = "observed",
                          reliability_x = rel_x, reliability_y = rel_y, seed = NULL)
      cor(d[[var_x]], d[[var_y]])
    }, numeric(1))
  }
  data.frame(study = seq_len(n_studies), n = n, statistic = statistic, value = values)
}

#' A "pile of dots" (ggdist::stat_dots, as in sampling-distributions.R)
#' for the sampling distribution of a chosen statistic, one panel per
#' sample size in `n_grid` -- the sampling-error figure: panels narrow
#' (spread shrinks) as N grows from top to bottom, but stay centred in
#' the same place.
#'
#' @inheritParams simulate_sampling_distribution
#' @param n_grid    sample sizes to compare, one facet panel each.
#' @param n_studies independent studies simulated *per* panel.
#' @return a ggplot object.
plot_sampling_distribution <- function(n_grid, n_studies = 300, statistic = c("mean", "cor"),
                                        var_x = "StuVoc1", var_y = NULL,
                                        reliability = NULL, true_corr = NULL, seed = NULL) {
  statistic <- match.arg(statistic)
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plot_sampling_distribution(). Install with install.packages('ggdist').")
  }
  if (!is.null(seed)) set.seed(seed)
  combined <- do.call(rbind, lapply(n_grid, function(n) {
    simulate_sampling_distribution(n, n_studies, statistic, var_x, var_y, reliability, true_corr, seed = NULL)
  }))
  n_levels <- sort(unique(combined$n))
  combined$n_label <- factor(combined$n, levels = n_levels, labels = paste("n =", n_levels))
  xlab <- if (statistic == "mean") sprintf("Sample mean (%s)", var_x) else sprintf("Sample correlation, r(%s, %s)", var_x, var_y)

  ggplot(combined, aes(x = value)) +
    ggdist::stat_dots(color = "#3B6EA5", fill = "#3B6EA5") +
    facet_wrap(~n_label, ncol = 1) +
    labs(x = xlab,
         subtitle = "Each dot is one simulated study's estimate \u2014 larger N narrows the spread (sampling error), not the centre") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
}

# ---- measurement error again: attenuation of a correlation ------------------

#' Hold the *construct-level* (true-score) correlation between two
#' tests fixed, and vary their (shared) reliability -- shows the
#' observed-score correlation shrinking toward zero as reliability
#' drops, purely from measurement error, with sample size untouched.
#'
#' @param reliability_grid  reliability values to sweep (applied to
#'                           both `var_x` and `var_y` equally).
#' @param n                 participants per simulated study (held
#'                           fixed across the sweep).
#' @param reps              independent studies simulated per
#'                           reliability value.
#' @param var_x, var_y      the focal pair of tests.
#' @param true_corr_latent  the *true* (disattenuated) correlation to
#'                           hold fixed; NULL (default) backs it out
#'                           from the Vermeiren-default *observed*
#'                           correlation and reliabilities for this
#'                           pair, via disattenuate().
#' @param seed              for reproducibility.
#' @return a data frame: `reliability`, `rep`, `r_observed`, with the
#'   fixed `true_corr_latent` value attached as an attribute.
simulate_attenuation_sweep <- function(reliability_grid, n = 300, reps = 200,
                                        var_x = "StuVoc1", var_y = "GK",
                                        true_corr_latent = NULL, seed = NULL) {
  stopifnot(
    "`var_x` and `var_y` must be different tests" = var_x != var_y,
    "`reliability_grid` values must all be strictly between 0 and 1" =
      all(reliability_grid > 0 & reliability_grid < 1),
    "`reps` must be at least 2" = reps >= 2
  )
  if (!is.null(seed)) set.seed(seed)
  if (is.null(true_corr_latent)) {
    rel_default      <- default_reliability_vec()
    r_obs_default    <- default_r_obs_matrix()[var_x, var_y]
    true_corr_latent <- disattenuate(r_obs_default, rel_default[[var_x]], rel_default[[var_y]])
  }
  grid <- expand.grid(reliability = reliability_grid, rep = seq_len(reps))
  grid$r_observed <- mapply(function(rel, rp) {
    d <- simulate_pair(n, var_x, var_y, true_corr = true_corr_latent, corr_scale = "latent",
                        reliability_x = rel, reliability_y = rel, seed = NULL)
    cor(d[[var_x]], d[[var_y]])
  }, grid$reliability, grid$rep)
  attr(grid, "true_corr_latent") <- true_corr_latent
  grid
}

#' The attenuation figure: observed r (y) against reliability (x), with
#' a ribbon showing the empirical spread across `reps` studies, the
#' theoretical attenuation line `r_obs = r_true * reliability`, and a
#' dotted reference line at the fixed true (construct-level)
#' correlation -- the ceiling you'd see with perfectly reliable
#' measurement.
#'
#' @inheritParams simulate_attenuation_sweep
#' @return a ggplot object.
plot_attenuation_vs_reliability <- function(reliability_grid, n = 300, reps = 200,
                                             var_x = "StuVoc1", var_y = "GK",
                                             true_corr_latent = NULL, seed = NULL) {
  raw <- simulate_attenuation_sweep(reliability_grid, n, reps, var_x, var_y, true_corr_latent, seed)
  tc  <- attr(raw, "true_corr_latent")

  summ_list <- lapply(split(raw$r_observed, raw$reliability), function(v) {
    data.frame(mean = mean(v), lo = unname(quantile(v, .025)), hi = unname(quantile(v, .975)))
  })
  summ_df <- do.call(rbind, summ_list)
  summ_df$reliability <- as.numeric(names(summ_list))
  summ_df <- summ_df[order(summ_df$reliability), ]

  curve_df <- data.frame(reliability = sort(unique(raw$reliability)))
  curve_df$theory <- tc * curve_df$reliability

  ggplot(summ_df, aes(reliability, mean)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#3B6EA5", alpha = 0.2) +
    geom_line(color = "#3B6EA5", linewidth = 1) +
    geom_point(color = "#3B6EA5", size = 2) +
    geom_line(data = curve_df, aes(reliability, theory), color = "#B23A48", linetype = "dashed") +
    geom_hline(yintercept = tc, color = "gray40", linetype = "dotted") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Reliability (both tests, held equal)",
         y = sprintf("Observed r(%s, %s)", var_x, var_y),
         subtitle = sprintf(
           paste("True (construct-level) correlation held fixed at %.2f (dotted line);",
                 "dashed line = attenuation formula, r_obs = r_true \u00d7 reliability"),
           tc)) +
    theme_minimal()
}
