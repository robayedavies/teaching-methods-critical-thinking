#' ============================================================
#' error-rates-process-control.R
#'
#' False-positive and false-negative rates as long-run *process-control*
#' error rates, in the Neyman-Pearson decision-procedure sense (Neyman &
#' Pearson, 1933) -- built in two parts on top of the battery simulated by
#' psych-battery-dgp.R (StuVoc1, StuVoc2, GK, ART3; grounded in Vermeiren,
#' Vandendaele & Brysbaert, 2023).
#'
#' PART A -- a concrete anchor. ART3's item-level task (recognize 60 real
#' authors, reject 30 foil names) already produces genuine false positives
#' (saying "yes" to a foil) and false negatives (saying "no" to a real
#' author) for every simulated participant -- no abstraction needed yet.
#' This script turns that into a person-level *decision rule*
#' ("classify someone as a frequent reader if they recognize at least k of
#' the 60 authors") and asks the Neyman-Pearson question directly: not
#' "was this one classification right?", but "if this exact rule were
#' applied over and over, how often would it be wrong, and in which
#' direction?" -- i.e. the rule's false-positive and false-negative rate
#' are *properties of the rule*, estimable from repetition, not
#' properties of any one person's test result.
#'
#' PART B -- the same logic, abstracted to NHST. A significance test at a
#' fixed alpha is exactly this kind of rule: "reject H0 if p < alpha."
#' Neyman & Pearson's own framing is explicit that alpha (their "size" of
#' the critical/rejection region) and power (1 - their "Type II" rate)
#' are long-run error rates of a fixed testing *procedure*, chosen in
#' advance of seeing the data and justified by their long-run behaviour
#' under repeated use -- not a measure of evidence for or against any
#' particular result, which was Fisher's (distinct, and often conflated
#' with Neyman-Pearson's) reading of a p-value. See Perezgonzalez (2015)
#' for an accessible tutorial contrasting the two traditions. This script
#' reuses psych-battery-dgp.R's run_sweep()/summarize_sweep() -- already
#' built to compute power, the empirical Type I/false-positive rate, sign
#' consistency (Type S), the Type M "exaggeration ratio" (how inflated a
#' significant estimate tends to be), and CI coverage -- varying all three
#' DGP dials (sample size, reliability, true correlation) to show that
#' same quartet behaving exactly like the ART3 classification rule's error
#' rates: stable long-run properties of a fixed procedure, not properties
#' of any one study.
#'
#' Requires: psych-battery-dgp.R sourced first (uses default_test_specs(),
#' default_reliability_vec(), get_calibration(), simulate_battery(),
#' run_sweep(), summarize_sweep(), plot_metrics_vs_x(),
#' plot_power_heatmap(), default_r_obs_matrix()); ggplot2.
#'
#' Quarto usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("error-rates-process-control.R")
#'
#'   # Part A: ART3 as a concrete false-positive/false-negative anchor
#'   plot_art3_rates(n = 500, seed = 1)
#'   plot_art3_roc(n = 300, k_grid = seq(30, 50, by = 5), reps = 150, seed_base = 2)
#'
#'   # Part B: the same error-rate quartet, abstracted to a fixed-alpha
#'   # NHST procedure, swept across all three DGP dials
#'   plot_reliability_heatmap("StuVoc1", "GK", n_grid = c(20, 50, 100, 200),
#'                             reliability_grid = seq(0.3, 0.9, by = 0.1), seed_base = 3)
#'   plot_true_corr_heatmap("StuVoc1", "GK", n_grid = c(20, 50, 100, 200),
#'                           true_corr_grid = c(0, 0.2, 0.67), seed_base = 4)
#'   plot_process_control_overview("StuVoc1", "GK", n_grid = c(20, 50, 100, 200, 400), seed_base = 5)
#'   plot_type_m_curve("StuVoc1", "GK", n_grid = c(10, 20, 50, 100, 200, 400), seed_base = 6)
#'   ```
#' ============================================================

library(ggplot2)

# ---- Part A: ART3 as a concrete false-positive/false-negative anchor -------

#' Simulate `n` participants on ART3 (60 author names, 30 foils), unpacking
#' each person's item-level hits/rejects and the false-positive/false-
#' negative rates they imply, plus their standardized true ability.
#'
#' @param n            number of simulated participants.
#' @param reliability  ART3's reliability; NULL (default) uses its
#'                      Vermeiren-default value.
#' @param seed         for reproducibility.
#' @return a data frame: `id`, `true_ability` (standardized, mean 0 / sd 1),
#'   `hits` (correct "yes" to an author, out of 60), `rejects` (correct
#'   "no" to a foil, out of 30), `score` (hits + rejects), `fn_rate`
#'   (proportion of authors missed), `fp_rate` (proportion of foils
#'   wrongly accepted), `reliability`.
simulate_art3_battery <- function(n, reliability = NULL, seed = NULL) {
  stopifnot("`n` must be a positive integer" = n >= 1)
  sp      <- default_test_specs()$ART3
  rel_val <- if (is.null(reliability)) unname(default_reliability_vec()[["ART3"]]) else reliability
  cal     <- get_calibration(sp$mu, rel_val, sp$n_items, sp$guessing)

  d <- simulate_battery(n, reliability = list(ART3 = rel_val), seed = seed, keep_latent = TRUE)
  hits         <- round((1 - d$ART3_fn_rate) * sp$split[["authors"]])
  rejects      <- round((1 - d$ART3_fp_rate) * sp$split[["foils"]])
  true_ability <- (d$ART3_theta - cal$theta_mu) / cal$theta_sigma

  data.frame(id = d$id, true_ability = true_ability, hits = hits, rejects = rejects,
             score = d$ART3, fn_rate = d$ART3_fn_rate, fp_rate = d$ART3_fp_rate,
             reliability = rel_val)
}

#' Each simulated participant's two item-level error rates, plotted
#' against each other -- the concrete starting point for this script:
#' these are real false positives (yes to a foil) and false negatives (no
#' to a real author), computed per person, with no NHST machinery
#' involved yet.
#'
#' @inheritParams simulate_art3_battery
#' @return a ggplot object.
plot_art3_rates <- function(n = 500, reliability = NULL, seed = NULL) {
  d <- simulate_art3_battery(n, reliability, seed)
  ggplot(d, aes(fp_rate, fn_rate)) +
    geom_point(color = "#3B6EA5", alpha = 0.5) +
    labs(x = "False-positive rate (said \"yes\" to a foil name)",
         y = "False-negative rate (said \"no\" to a real author)",
         subtitle = sprintf(
           "ART3, n = %d, reliability = %.2f -- each point is one simulated participant",
           n, d$reliability[1])) +
    theme_minimal()
}

#' One application of a person-level classification rule: simulate `n`
#' participants on ART3, classify each as a "frequent reader" if they
#' recognize at least `k` of the 60 authors, and compare that call
#' against each person's *true* ability (above or below the population
#' average true recognition ability, i.e. `ART3_theta > theta_mu`) -- a
#' real, if simplified, ground truth, independent of the noisy
#' classification itself.
#'
#' @param n            number of simulated participants (one application
#'                       of the rule).
#' @param k            the rule's threshold: classify as positive if
#'                       `hits >= k` (0-60).
#' @param reliability  ART3's reliability; NULL (default) uses its
#'                      Vermeiren-default value.
#' @param seed         for reproducibility.
#' @return a one-row data frame: `n`, `k`, `reliability`, `n_pos`/`n_neg`
#'   (how many people were truly positive/negative in this sample),
#'   `tpr`/`fnr` (true positive / false negative rate, among the truly
#'   positive), `fpr`/`tnr` (false positive / true negative rate, among
#'   the truly negative).
evaluate_art3_rule <- function(n, k, reliability = NULL, seed = NULL) {
  stopifnot(
    "`n` must be a positive integer" = n >= 1,
    "`k` must be between 0 and 60 (the number of author items in ART3)" = k >= 0 && k <= 60
  )
  sp  <- default_test_specs()$ART3
  rel <- if (is.null(reliability)) unname(default_reliability_vec()[["ART3"]]) else reliability
  cal <- get_calibration(sp$mu, rel, sp$n_items, sp$guessing)

  d <- simulate_battery(n, reliability = list(ART3 = rel), seed = seed, keep_latent = TRUE)
  hits               <- round((1 - d$ART3_fn_rate) * sp$split[["authors"]])
  truth_positive     <- d$ART3_theta > cal$theta_mu
  predicted_positive <- hits >= k

  tp <- sum(truth_positive & predicted_positive)
  fn <- sum(truth_positive & !predicted_positive)
  fp <- sum(!truth_positive & predicted_positive)
  tn <- sum(!truth_positive & !predicted_positive)

  data.frame(n = n, k = k, reliability = rel, n_pos = tp + fn, n_neg = fp + tn,
             tpr = tp / (tp + fn), fnr = fn / (tp + fn),
             fpr = fp / (fp + tn), tnr = tn / (fp + tn))
}

#' Repeat `evaluate_art3_rule()` `reps` times at each threshold in
#' `k_grid`, drawing a *fresh* sample of `n` people each time -- the
#' device that turns "one classification's outcome" into "the rule's
#' long-run operating characteristics". Deterministic per-cell seeding
#' mirrors psych-battery-dgp.R's run_sweep().
#'
#' @param n            participants per replicate application of the rule.
#' @param k_grid       thresholds to sweep.
#' @param reps         independent applications of the rule per threshold.
#' @param reliability  ART3's reliability; NULL (default) uses its
#'                      Vermeiren-default value.
#' @param seed_base    base seed for deterministic per-cell seeding.
#' @return a tidy data frame: one row per (k, rep), same columns as
#'   `evaluate_art3_rule()` plus `rep`.
sweep_art3_rule <- function(n, k_grid, reps = 200, reliability = NULL, seed_base = 1234) {
  stopifnot("`reps` must be at least 1" = reps >= 1)
  grid <- expand.grid(k = k_grid, rep = seq_len(reps))
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    seed_i  <- seed_base + grid$rep[i] * 1e5 + grid$k[i]
    rows[[i]] <- cbind(rep = grid$rep[i],
                        evaluate_art3_rule(n, grid$k[i], reliability = reliability, seed = seed_i))
  }
  do.call(rbind, rows)
}

#' The ART3 classification rule's ROC-style error-rate trade-off: pale
#' points are individual applications of the rule (one per replicate,
#' noisy -- a single sample's apparent false-positive/true-positive rate
#' can land almost anywhere nearby), the red line connects each
#' threshold's *average* rate across all replicates -- the rule's actual
#' long-run operating curve. Tightening `k` (fewer people classified
#' "frequent reader") trades a lower false-positive rate for a lower true
#' positive rate (more false negatives), and there is no threshold that
#' minimizes both at once.
#'
#' @inheritParams sweep_art3_rule
#' @return a ggplot object.
plot_art3_roc <- function(n, k_grid, reps = 200, reliability = NULL, seed_base = 1234) {
  raw  <- sweep_art3_rule(n, k_grid, reps, reliability, seed_base)
  summ <- aggregate(cbind(fpr, tpr) ~ k, data = raw, FUN = function(x) mean(x, na.rm = TRUE))
  summ <- summ[order(summ$k), ]

  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
    geom_point(data = raw, aes(fpr, tpr), color = "#3B6EA5", alpha = 0.15, size = 1) +
    geom_path(data = summ, aes(fpr, tpr), color = "#B23A48", linewidth = 1) +
    geom_point(data = summ, aes(fpr, tpr), color = "#B23A48", size = 2) +
    geom_text(data = summ, aes(fpr, tpr, label = k), color = "#B23A48", nudge_y = 0.03, size = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "False-positive rate (classified \"frequent reader\", isn't)",
         y = "True-positive rate (classified \"frequent reader\", is)",
         subtitle = sprintf(
           "n = %d per replicate, %d replicates per threshold k -- pale dots: one single application of the rule; red line: the rule's long-run operating curve",
           n, reps)) +
    theme_minimal()
}

# ---- Part B: the same logic, abstracted to a fixed-alpha NHST procedure ----

#' A sweep engine paralleling psych-battery-dgp.R's run_sweep(), but for
#' dial (ii), reliability, swept down to values low enough that holding
#' the *observed*-score correlation fixed (run_sweep()'s usual
#' "observed"-scale interpretation of `true_corr`) becomes impossible --
#' disattenuating an observed r back through a low shared reliability can
#' require a latent correlation above 1, which `disattenuate()` can only
#' handle by capping it (with a warning), silently flattening the figure
#' at low reliability. The fix is the same device
#' `sampling-measurement-error.R`'s `simulate_attenuation_sweep()` uses:
#' hold the *latent* (construct-level) correlation fixed instead, and let
#' the *observed* correlation -- and therefore its significance, sign,
#' and coverage -- respond to reliability, which is the whole point of
#' this figure. The value actually tested for significance/sign/coverage
#' at each cell is the *implied* observed correlation,
#' `true_corr_latent * reliability` (both tests share one reliability
#' value here) -- recomputed per cell and written to the output's
#' `true_corr` column so `summarize_sweep()`/`plot_power_heatmap()` can
#' be reused unmodified.
#'
#' @param var_x, var_y      the focal test pair.
#' @param n_grid             sample sizes to sweep.
#' @param reliability_grid   reliability values to sweep (applied to both
#'                            tests equally).
#' @param true_corr_latent   the *true* (disattenuated) correlation to
#'                            hold fixed; NULL (default) backs it out
#'                            from the pair's Vermeiren-default observed
#'                            correlation and reliabilities, via
#'                            `disattenuate()`.
#' @param reps               replicate studies per grid cell.
#' @param alpha              significance threshold.
#' @param seed_base          base seed for deterministic per-cell seeding.
#' @return a tidy data frame in the same shape as `run_sweep()`'s output
#'   (so `summarize_sweep()` can be reused directly), with the fixed
#'   `true_corr_latent` value attached as an attribute.
run_reliability_sweep <- function(var_x, var_y, n_grid, reliability_grid, true_corr_latent = NULL,
                                   reps = 200, alpha = 0.05, seed_base = 1234) {
  if (is.null(true_corr_latent)) {
    rel_default      <- default_reliability_vec()
    r_obs_default    <- default_r_obs_matrix()[var_x, var_y]
    true_corr_latent <- disattenuate(r_obs_default, rel_default[[var_x]], rel_default[[var_y]])
  }
  grid <- expand.grid(n = n_grid, reliability = reliability_grid, rep = seq_len(reps))
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    n <- grid$n[i]; rel <- grid$reliability[i]; rep_i <- grid$rep[i]
    seed_cell   <- seed_base + rep_i * 1e6 + n * 1e3 + round(rel * 1000)
    implied_obs <- true_corr_latent * rel

    d    <- simulate_pair(n, var_x, var_y, true_corr = true_corr_latent, corr_scale = "latent",
                           reliability_x = rel, reliability_y = rel, seed = seed_cell)
    fit  <- fit_focal_cor(d, var_x, var_y)
    diag <- diagnose_cor_result(fit, implied_obs, alpha)
    rows[[i]] <- cbind(data.frame(n = n, true_corr = implied_obs, reliability = rel, rep = rep_i), diag)
  }
  out <- do.call(rbind, rows)
  attr(out, "true_corr_latent") <- true_corr_latent
  out
}

#' Power/Type I rate, sign consistency, or CI coverage across sample size
#' and reliability, with the pair's *true* (construct-level) correlation
#' held fixed -- see `run_reliability_sweep()` for why that, and not the
#' observed-score correlation, is what has to stay fixed here.
#'
#' @param var_x, var_y      the focal test pair.
#' @param n_grid             sample sizes to sweep.
#' @param reliability_grid   reliability values to sweep (applied to both
#'                            tests equally).
#' @param true_corr_latent   the *true* (disattenuated) correlation to
#'                            hold fixed; NULL (default) backs it out
#'                            from the pair's Vermeiren-default observed
#'                            correlation and reliabilities.
#' @param metric             one of `"power_or_type1"`, `"sign_consistency"`,
#'                            `"coverage"`.
#' @param reps               replicate studies per grid cell.
#' @param alpha              significance threshold.
#' @param seed_base          base seed, passed through to
#'                            `run_reliability_sweep()`.
#' @return a ggplot object.
plot_reliability_heatmap <- function(var_x, var_y, n_grid, reliability_grid, true_corr_latent = NULL,
                                      metric = "power_or_type1", reps = 200, alpha = 0.05,
                                      seed_base = 1234) {
  stopifnot(
    "`metric` must be one of power_or_type1, sign_consistency, coverage" =
      metric %in% c("power_or_type1", "sign_consistency", "coverage")
  )
  sweep <- run_reliability_sweep(var_x, var_y, n_grid, reliability_grid, true_corr_latent,
                                  reps, alpha, seed_base)
  summ  <- summarize_sweep(sweep)
  tc    <- attr(sweep, "true_corr_latent")

  fill_name <- switch(metric, power_or_type1 = "Power / Type I rate",
                       sign_consistency = "Sign consistency", coverage = "CI coverage")
  plot_power_heatmap(summ, x = "n", y = "reliability", fill = metric, fill_name = fill_name) +
    labs(subtitle = sprintf(
      "%s vs. %s, true (construct-level) correlation held fixed at %.2f -- the observed correlation (and its significance) responds to reliability",
      var_x, var_y, tc))
}

#' Power/Type I rate, sign consistency, or CI coverage across sample size
#' and true correlation (including 0, the empirical false-positive
#' benchmark), at one fixed reliability.
#'
#' @inheritParams plot_reliability_heatmap
#' @param true_corr_grid  population *observed*-score correlations to
#'                         sweep; include 0 to see the empirical
#'                         false-positive rate alongside power.
#' @param reliability     reliability applied to both tests; NULL
#'                         (default) uses `var_x`'s Vermeiren-default
#'                         value.
#' @return a ggplot object.
plot_true_corr_heatmap <- function(var_x, var_y, n_grid, true_corr_grid, reliability = NULL,
                                    metric = "power_or_type1", reps = 200, alpha = 0.05,
                                    seed_base = 1234) {
  stopifnot(
    "`metric` must be one of power_or_type1, sign_consistency, coverage" =
      metric %in% c("power_or_type1", "sign_consistency", "coverage")
  )
  rel_grid <- if (is.null(reliability)) NULL else reliability
  sweep <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr_grid, reliability_grid = rel_grid,
                      reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  summ  <- summarize_sweep(sweep)

  fill_name <- switch(metric, power_or_type1 = "Power / Type I rate",
                       sign_consistency = "Sign consistency", coverage = "CI coverage")
  plot_power_heatmap(summ, x = "n", y = "true_corr", fill = metric, fill_name = fill_name)
}

#' All three error-rate diagnostics at once -- power/Type I rate, sign
#' consistency, and CI coverage -- faceted against sample size, the same
#' "quartet" framing used throughout this project's heavier Aims 1-2
#' simulation work, here presented explicitly as long-run properties of
#' one fixed testing procedure (Neyman & Pearson, 1933) rather than as
#' evidence about any one study.
#'
#' @inheritParams plot_true_corr_heatmap
#' @return a ggplot object.
plot_process_control_overview <- function(var_x, var_y, n_grid, true_corr_grid = NULL,
                                           reliability = NULL, reps = 300, alpha = 0.05,
                                           seed_base = 1234) {
  if (is.null(true_corr_grid)) {
    true_corr_grid <- c(0, unname(default_r_obs_matrix()[var_x, var_y]))
  }
  rel_grid <- if (is.null(reliability)) NULL else reliability
  sweep <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr_grid, reliability_grid = rel_grid,
                      reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  summ  <- summarize_sweep(sweep)

  plot_metrics_vs_x(summ, x = "n", metrics = c("power_or_type1", "sign_consistency", "coverage"),
                     color_by = "true_corr") +
    labs(subtitle = sprintf(
      "%s vs. %s -- power/Type I, sign consistency, and CI coverage as long-run rates of one fixed decision procedure",
      var_x, var_y))
}

#' The Type M ("exaggeration ratio") diagnostic: the median
#' `|estimate| / |true_corr|` among *significant* results only, across
#' sample size -- how inflated a reported effect tends to be, conditional
#' on the significance filter having let it through. Deliberately a
#' *separate* figure from `plot_process_control_overview()`'s three
#' panels, not a fourth one folded into it: `type_m_ratio` is not a 0-1
#' proportion like power/sign-consistency/coverage (it is 1 at no
#' inflation and can comfortably exceed 2-3 at low power), so sharing
#' `plot_metrics_vs_x()`'s fixed `c(0, 1)` y-axis would silently clip the
#' part of the story that matters most -- how much *worse* the
#' exaggeration gets as power falls. `true_corr_grid` defaults to just the
#' pair's own Vermeiren-default correlation (not `0`, unlike
#' `plot_power_curves()`/`plot_process_control_overview()`): `type_m_ratio`
#' is `NA` under the null by construction (see `diagnose_cor_result()`),
#' so a `true_corr = 0` line here would just be an empty row with nothing
#' to show.
#'
#' @inheritParams plot_true_corr_heatmap
#' @param true_corr_grid  population *observed*-score correlation(s) to
#'                         compute the exaggeration ratio for; NULL
#'                         (default) uses the pair's Vermeiren-default
#'                         value only (`0` is never useful here -- see
#'                         above).
#' @return a ggplot object.
plot_type_m_curve <- function(var_x, var_y, n_grid, true_corr_grid = NULL,
                               reliability = NULL, reps = 300, alpha = 0.05,
                               seed_base = 1234) {
  if (is.null(true_corr_grid)) {
    true_corr_grid <- unname(default_r_obs_matrix()[var_x, var_y])
  }
  rel_grid <- if (is.null(reliability)) NULL else reliability
  sweep <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr_grid, reliability_grid = rel_grid,
                      reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  summ  <- summarize_sweep(sweep)

  ggplot(summ, aes(n, type_m_ratio, color = factor(true_corr), group = factor(true_corr))) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray40") +
    geom_line() + geom_point() +
    labs(x = "n", y = "Type M ratio (median |estimate| / |true correlation|, among significant results)",
         color = "True correlation",
         subtitle = sprintf(
           "%s vs. %s -- the significance filter's exaggeration ratio, falling toward 1 (dotted line) as power rises",
           var_x, var_y)) +
    theme_minimal()
}
