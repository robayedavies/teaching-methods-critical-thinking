#' ============================================================
#' nhst-power.R
#'
#' Null hypothesis significance testing (NHST) and power, demonstrated on
#' two correlation pairs from the battery simulated by psych-battery-dgp.R
#' (StuVoc1, StuVoc2, GK, ART3; grounded in Vermeiren, Vandendaele &
#' Brysbaert, 2023):
#'
#'   StuVoc1 -- StuVoc2   observed r ≈ .82 (two related vocabulary subtests)
#'   StuVoc1 -- GK        observed r ≈ .67 (vocabulary and general knowledge)
#'
#' Carrying two correlations with genuinely different effect sizes through
#' the same machinery is the point: the StuVoc1-StuVoc2 pair reaches high
#' power at small sample sizes, StuVoc1-GK needs a good deal more -- the
#' same test, the same decision rule, very different N requirements,
#' because power depends on effect size as much as on sample size.
#'
#' This script adds no new simulation machinery of its own -- every
#' function here is a thin orchestration layer over psych-battery-dgp.R's
#' simulate_pair()/run_sweep()/summarize_sweep()/plot_metrics_vs_x(), in
#' the order:
#'
#'   one simulated study  -> summarize_single_study() / plot_single_study()
#'   many studies, one N  -> plot_corridor()   (replicates Vermeiren et al.'s
#'                                               own Fig. 1 "corridor" plot)
#'   many studies, many N -> plot_power_curves()
#'
#' The corridor and power-curve figures both default to sweeping true_corr
#' = c(0, <the pair's Vermeiren-default observed correlation>) -- 0 so the
#' empirical false-positive (Type I) rate is always visible alongside
#' power, which is the explicit point raised in the walkthrough: "the
#' study wasn't significant" and "the null is true" are not the same
#' thing, and a flat curve at 0 pinned near alpha is what a well-behaved
#' *false* positive rate looks like, not evidence either way about any
#' particular pair.
#'
#' Requires: psych-battery-dgp.R sourced first (uses default_r_obs_matrix(),
#' simulate_pair(), fit_focal_cor(), diagnose_cor_result(), run_sweep(),
#' summarize_sweep(), plot_metrics_vs_x()); ggplot2.
#'
#' Quarto usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("nhst-power.R")
#'
#'   # one study, one decision
#'   summarize_single_study("StuVoc1", "GK", n = 200, seed = 1)
#'   plot_single_study("StuVoc1", "GK", n = 200, seed = 1)
#'
#'   # many studies at one N -- the Vermeiren Fig. 1 style "corridor"
#'   plot_corridor("StuVoc1", "GK", n_grid = c(40, 100, 200, 400), reps = 150, seed = 2)
#'
#'   # power curves for both focal pairs
#'   plot_power_curves("StuVoc1", "StuVoc2", n_grid = c(10, 20, 40, 80), reps = 300, seed = 3)
#'   plot_power_curves("StuVoc1", "GK", n_grid = c(20, 50, 100, 200, 400), reps = 300, seed = 4)
#'   ```
#' ============================================================

library(ggplot2)

# ---- shared helper: turning a p-value into a plain-language decision -------

#' A one-line plain-language NHST decision, given a p-value and alpha.
#' Used by both summarize_single_study() and plot_single_study() so the
#' wording is identical everywhere it appears.
format_pvalue <- function(p_value) {
  if (p_value < 0.001) "p < .001" else sprintf("p = %.3f", p_value)
}

interpret_decision <- function(p_value, alpha = 0.05) {
  p_label <- format_pvalue(p_value)
  if (p_value < alpha) {
    sprintf("reject H0 (%s, threshold %.2f) -- call this result \"significant\"", p_label, alpha)
  } else {
    sprintf("fail to reject H0 (%s, threshold %.2f) -- call this result \"not significant\"", p_label, alpha)
  }
}

# ---- one study ---------------------------------------------------------------

#' Simulate one study of `n` participants on a focal test pair, fit a
#' Pearson correlation test, and summarize it as a single-row data frame
#' with a plain-language decision attached.
#'
#' @param var_x, var_y    the focal pair of tests (e.g. "StuVoc1", "GK").
#' @param n                sample size.
#' @param true_corr        population *observed*-score correlation to
#'                          simulate from; NULL (default) uses the
#'                          Vermeiren-default value for this pair.
#' @param reliability      reliability override applied to both tests;
#'                          NULL (default) uses each test's own
#'                          Vermeiren-default value.
#' @param alpha            significance threshold.
#' @param seed             for reproducibility.
#' @return a data frame: `var_x`, `var_y`, `n`, `true_corr`, `estimate`,
#'   `ci_low`, `ci_high`, `p_value`, `significant`, `decision` (the
#'   plain-language string from `interpret_decision()`).
summarize_single_study <- function(var_x, var_y, n, true_corr = NULL, reliability = NULL,
                                    alpha = 0.05, seed = NULL) {
  tc <- if (is.null(true_corr)) default_r_obs_matrix()[var_x, var_y] else true_corr
  d   <- simulate_pair(n, var_x, var_y, true_corr = tc, corr_scale = "observed",
                        reliability_x = reliability, reliability_y = reliability, seed = seed)
  fit <- fit_focal_cor(d, var_x, var_y)
  data.frame(var_x = var_x, var_y = var_y, n = n, true_corr = tc,
             estimate = fit$estimate, ci_low = fit$ci_low, ci_high = fit$ci_high,
             p_value = fit$p_value, significant = fit$p_value < alpha,
             decision = interpret_decision(fit$p_value, alpha))
}

#' Scatterplot of one simulated study's two test scores, with the fitted
#' line and the NHST decision (estimate, 95% CI, p-value, significant?)
#' annotated directly on the figure -- the single clearest picture of
#' "what a p-value is attached to": one scatter of points, one line fit
#' through them, one yes/no call about whether that line's slope could
#' plausibly be flat.
#'
#' @inheritParams summarize_single_study
#' @return a ggplot object.
plot_single_study <- function(var_x, var_y, n, true_corr = NULL, reliability = NULL,
                               alpha = 0.05, seed = NULL) {
  tc <- if (is.null(true_corr)) default_r_obs_matrix()[var_x, var_y] else true_corr
  d   <- simulate_pair(n, var_x, var_y, true_corr = tc, corr_scale = "observed",
                        reliability_x = reliability, reliability_y = reliability, seed = seed)
  fit <- fit_focal_cor(d, var_x, var_y)

  ggplot(d, aes(.data[[var_x]], .data[[var_y]])) +
    geom_point(color = "#3B6EA5", alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "#B23A48") +
    labs(x = var_x, y = var_y,
         subtitle = sprintf(
           "n = %d \u2014 r = %.2f, 95%% CI [%.2f, %.2f], %s\n%s",
           n, fit$estimate, fit$ci_low, fit$ci_high, format_pvalue(fit$p_value),
           interpret_decision(fit$p_value, alpha))) +
    theme_minimal()
}

# ---- many studies, one N: the "corridor" plot --------------------------------

#' A Vermeiren-et-al.-Fig.-1-style "corridor" plot: many independent
#' simulated studies at each of several sample sizes, one point per
#' study, plotted as (estimate, n) and coloured by whether that study's
#' result was significant. Directly reuses psych-battery-dgp.R's
#' run_sweep() rather than re-simulating -- this figure is just a
#' different view of run_sweep()'s raw (per-replicate) output, not a new
#' simulation.
#'
#' @param var_x, var_y   the focal pair of tests.
#' @param n_grid          sample sizes to compare.
#' @param reps            independent studies simulated per sample size.
#' @param true_corr       population *observed*-score correlation to
#'                         simulate from; NULL (default) uses the
#'                         Vermeiren-default value for this pair.
#' @param reliability     reliability override applied to both tests;
#'                         NULL (default) uses Vermeiren-default values.
#' @param alpha           significance threshold.
#' @param seed_base       base seed, passed through to run_sweep() for
#'                         deterministic per-cell seeding.
#' @return a ggplot object.
plot_corridor <- function(var_x, var_y, n_grid, reps = 150, true_corr = NULL,
                           reliability = NULL, alpha = 0.05, seed_base = 1234) {
  tc <- if (is.null(true_corr)) unname(default_r_obs_matrix()[var_x, var_y]) else true_corr
  rel_grid <- if (is.null(reliability)) NULL else reliability

  raw <- run_sweep(n_grid = n_grid, true_corr_grid = tc, reliability_grid = rel_grid,
                    reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  raw$significant_label <- factor(raw$significant, levels = c(FALSE, TRUE),
                                   labels = c("not significant", "significant"))

  ggplot(raw, aes(x = estimate, y = factor(n), color = significant_label, shape = significant_label)) +
    geom_vline(xintercept = tc, color = "gray40", linetype = "dotted") +
    geom_jitter(height = 0.15, alpha = 0.6) +
    scale_color_manual(values = c("not significant" = "#3B6EA5", "significant" = "#B23A48"), name = NULL) +
    scale_shape_manual(values = c("not significant" = 16, "significant" = 17), name = NULL) +
    labs(x = sprintf("Observed r(%s, %s)", var_x, var_y), y = "Sample size (n)",
         subtitle = sprintf(
           "%d simulated studies per n; true (population) correlation = %.2f (dotted line)",
           reps, tc)) +
    theme_minimal()
}

# ---- many studies, many N: power curves --------------------------------------

#' Power curves (and, at true_corr = 0, the empirical Type I/false-
#' positive rate) across a sample-size grid, for one focal pair. A thin
#' wrapper around run_sweep() + summarize_sweep() + psych-battery-dgp.R's
#' plot_metrics_vs_x(), with a dashed reference line at `alpha` added so
#' the true_corr = 0 curve's relationship to the nominal false-positive
#' rate is visible directly on the figure, not just inferable from where
#' the curve happens to sit.
#'
#' @param var_x, var_y     the focal pair of tests.
#' @param n_grid            sample sizes to sweep.
#' @param true_corr_grid    population *observed*-score correlations to
#'                          compare; NULL (default) uses
#'                          `c(0, <this pair's Vermeiren-default r>)`, so
#'                          the Type I curve is always included alongside
#'                          the pair's own realistic power curve.
#' @param reliability_grid  reliability override applied to both tests;
#'                          NULL (default) uses Vermeiren-default values.
#' @param reps              replicate studies per grid cell.
#' @param alpha             significance threshold.
#' @param seed_base         base seed, passed through to run_sweep().
#' @return a ggplot object.
plot_power_curves <- function(var_x, var_y, n_grid, true_corr_grid = NULL,
                               reliability_grid = NULL, reps = 300, alpha = 0.05,
                               seed_base = 1234) {
  if (is.null(true_corr_grid)) {
    true_corr_grid <- c(0, unname(default_r_obs_matrix()[var_x, var_y]))
  }
  sweep <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr_grid,
                      reliability_grid = reliability_grid, reps = reps,
                      var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  summ <- summarize_sweep(sweep)

  plot_metrics_vs_x(summ, x = "n", metrics = "power_or_type1", color_by = "true_corr") +
    geom_hline(yintercept = alpha, linetype = "dotted", color = "gray40") +
    labs(y = "Power (true_corr != 0)  /  false-positive rate (true_corr = 0)",
         subtitle = sprintf(
           "%s vs. %s \u2014 dotted line = nominal alpha (%.2f); %d replicate studies per cell",
           var_x, var_y, alpha, reps),
         color = "True correlation") +
    theme_minimal()
}
