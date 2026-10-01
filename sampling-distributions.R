#' ============================================================
#' sampling-distributions.R
#'
#' Simulates scores on a hypothetical 0-40 item vocabulary test from a
#' trial level up, and draws the figure sequence: one dot -> two dots ->
#' a "pile of dots" showing a sample's distribution -> a histogram -> a
#' grid of histograms (repeated studies) -> a scatterplot against a
#' continuous covariate.
#'
#' Model (trial level up, so every observed score is *structurally*
#' bounded to [0, n_items] -- no post-hoc clipping is needed):
#'
#'   Ability (logit scale):  theta_i     ~ N(theta_mu, theta_sigma)
#'   Per-item probability:   p_i         =  plogis(theta_i)
#'   Trial-level responses:  Y_i1..Y_i40 ~  Bernoulli(p_i)   (independent)
#'   Observed score:         score_i     =  sum_j Y_ij        in {0,...,40}
#'
#' `theta_mu` and `theta_sigma` aren't set directly -- they're solved
#' numerically (see `calibrate_ability()`) so that the *emergent* mean
#' score matches `mu`, and the emergent parallel-forms reliability
#' (corr(score, score') for two independent administrations drawn from
#' the same theta_i, i.e. the same person re-tested) matches
#' `reliability`. This mirrors real test theory: for a fixed test
#' length, how reliable a test is and how spread out scores are, are
#' linked to each other and to the mean, not independently free knobs --
#' fixing `mu` and one of {`reliability`, `sd_true`} determines the
#' other. `calibrate_ability()` (and therefore `simulate_vocab_data()` /
#' `simulate_parallel_forms()`) accepts either: pass `reliability` (the
#' default) when a figure is about repeatability, or pass `sd_true`
#' directly when a figure needs independent control of spread (e.g. two
#' panels with the same mean but visibly different spread). Whichever
#' one wasn't the calibration target is still available afterwards as an
#' emergent value via `calibrate_ability()$moments`.
#'
#' `r_with_x` is specified on the *observed*-score scale, i.e. it's
#' (approximately) the correlation you'll actually see in
#' `plot_scatter()`. Because measurement error attenuates any
#' correlation involving a noisy observed variable, hitting a target
#' observed-score correlation requires a *larger* correlation between
#' the underlying ability and x; classical test theory gives
#' corr(observed, true ability) = sqrt(reliability), so that attenuation
#' is inverted internally. If `r_with_x` exceeds what's achievable at
#' the requested `reliability` (i.e. `abs(r_with_x) > sqrt(reliability)`),
#' it's capped and a warning is issued.
#'
#' Requires: ggplot2, ggdist (for the dot-pile figure).
#'
#' Quarto usage:
#'   ```{r}
#'   source("sampling-distributions.R")
#'   d <- simulate_vocab_data(n = 40, seed = 1)
#'   plot_dots(d$observed_score[1])
#'   plot_dots(d$observed_score[1:2])
#'   plot_dotpile(d$observed_score)
#'   plot_histogram(d$observed_score)
#'   plot_histogram_grid(simulate_studies(6, 40, seed = 2))
#'   plot_scatter(simulate_vocab_data(n = 60, r_with_x = 0.5, seed = 3))
#'   ```
#' ============================================================

library(ggplot2)

# ---- calibration ------------------------------------------------------

#' Expectation of `f(theta)` under theta ~ N(theta_mu, theta_sigma),
#' via numerical integration (base R, no extra dependency).
expect_theta <- function(f, theta_mu, theta_sigma) {
  integrate(function(theta) f(theta) * dnorm(theta, theta_mu, theta_sigma),
            lower = theta_mu - 10 * theta_sigma,
            upper = theta_mu + 10 * theta_sigma)$value
}

#' Moments of the observed score implied by a given (theta_mu,
#' theta_sigma, n_items), under the trial-level model described above.
#' `var_error` uses the binomial variance n*p*(1-p) as the expected
#' within-person sampling variance across the 40 items.
implied_moments <- function(theta_mu, theta_sigma, n_items) {
  p_bar  <- expect_theta(plogis, theta_mu, theta_sigma)
  p2_bar <- expect_theta(function(th) plogis(th)^2, theta_mu, theta_sigma)
  pq_bar <- expect_theta(function(th) plogis(th) * (1 - plogis(th)), theta_mu, theta_sigma)

  mean_score  <- n_items * p_bar
  var_true    <- n_items^2 * (p2_bar - p_bar^2)
  var_error   <- n_items * pq_bar
  reliability <- var_true / (var_true + var_error)

  list(mean_score = mean_score, var_true = var_true, var_error = var_error,
       reliability = reliability, sd_true = sqrt(var_true),
       sd_total = sqrt(var_true + var_error))
}

#' Solve for the (theta_mu, theta_sigma) of the ability distribution
#' that make the trial-level model's emergent mean score match `mu`,
#' and match *either* a target `reliability` *or* a target `sd_true`
#' (the emergent SD of the "true"/expected score across people) --
#' whichever is supplied. These two are not independently settable: for
#' a fixed mu and n_items, specifying one determines the other, so pick
#' whichever is more natural for the figure you're building (reliability
#' for "how repeatable is this test", sd_true for "how spread out is
#' this population").
#'
#' @param mu           target mean observed score (0 - n_items).
#' @param reliability  target parallel-forms / test-retest reliability
#'                      (0-1, exclusive). Used only when `sd_true` is
#'                      NULL. Values above about 0.99 are not achievable
#'                      for a 40-item test (they'd require an implausibly
#'                      large ability spread) and raise an informative
#'                      error rather than silently returning nonsense.
#' @param sd_true      target SD, on the 0-n_items score scale, of
#'                      people's expected/"true" scores (i.e. the
#'                      between-person spread, excluding item-sampling
#'                      noise). When supplied, this is calibrated
#'                      directly and `reliability` is ignored. Must be
#'                      less than the theoretical ceiling
#'                      `sqrt(mu * (n_items - mu))` (approached only as
#'                      the ability spread goes to infinity); requesting
#'                      a value at or above that raises an informative
#'                      error.
#' @param n_items       number of trials/items (default 40).
#' @return list(theta_mu, theta_sigma, moments) -- `moments` is the
#'   full set of implied moments at the solution, including the
#'   emergent `sd_true`/`sd_total`/`reliability` (whichever wasn't the
#'   target), for inspection.
calibrate_ability <- function(mu, reliability = 0.80, sd_true = NULL, n_items = 40) {
  stopifnot("`mu` must be strictly between 0 and `n_items`" = mu > 0 && mu < n_items)

  theta_mu_for_sigma <- function(theta_sigma) {
    uniroot(function(theta_mu) implied_moments(theta_mu, theta_sigma, n_items)$mean_score - mu,
            lower = -15, upper = 15, extendInt = "yes")$root
  }

  # Robustly solve for the theta_sigma that zeroes `gap_fn`, widening the
  # search bracket by hand rather than relying on uniroot's own
  # extendInt = "yes" -- which occasionally overshoots into a theta_sigma
  # so large that the inner theta_mu search and integrate() start
  # returning NaN, or (near mid-scale mu) integrate() itself reports a
  # divergent integral. A target that's only reachable out there implies
  # an implausibly large, numerically unstable ability spread, so both
  # failure modes are treated the same way: fail with an informative
  # message instead of propagating a cryptic numerical error.
  solve_theta_sigma <- function(gap_fn, fail_message) {
    gap_at <- function(theta_sigma) tryCatch(gap_fn(theta_sigma), error = function(e) NA_real_)
    lower <- 0.02
    upper <- 8
    max_upper <- 128
    g <- gap_at(upper)
    while ((is.na(g) || g < 0) && upper < max_upper) {
      upper <- upper * 2
      g <- gap_at(upper)
    }
    if (is.na(g) || g < 0) stop(fail_message)
    tryCatch(uniroot(gap_fn, lower = lower, upper = upper)$root,
             error = function(e) stop(fail_message))
  }

  if (!is.null(sd_true)) {
    stopifnot("`sd_true` must be positive" = sd_true > 0)
    max_sd <- sqrt(mu * (n_items - mu))
    if (sd_true >= max_sd) {
      stop(sprintf(
        paste("sd_true = %.2f is not achievable at mu = %.1f with",
              "n_items = %d: the theoretical ceiling (approached only as",
              "the ability spread -> infinity) is %.2f. Try a smaller sd_true."),
        sd_true, mu, n_items, max_sd
      ))
    }
    sd_true_gap <- function(theta_sigma) {
      theta_mu <- theta_mu_for_sigma(theta_sigma)
      implied_moments(theta_mu, theta_sigma, n_items)$sd_true - sd_true
    }
    theta_sigma <- solve_theta_sigma(sd_true_gap, sprintf(
      paste("sd_true = %.2f is not achievable at mu = %.1f with",
            "n_items = %d (it would require a numerically unstable",
            "ability spread, even though it's below the theoretical",
            "ceiling of %.2f). Try a smaller sd_true."),
      sd_true, mu, n_items, max_sd
    ))
  } else {
    stopifnot(
      "`reliability` must be strictly between 0 and 1" =
        reliability > 0 && reliability < 1
    )
    reliability_gap <- function(theta_sigma) {
      theta_mu <- theta_mu_for_sigma(theta_sigma)
      implied_moments(theta_mu, theta_sigma, n_items)$reliability - reliability
    }
    theta_sigma <- solve_theta_sigma(reliability_gap, sprintf(
      paste("reliability = %.3f is not achievable at mu = %.1f with",
            "n_items = %d (it would require an implausibly large,",
            "numerically unstable ability spread). Try a lower",
            "reliability (<= ~0.99 is usually safe) or a larger n_items."),
      reliability, mu, n_items
    ))
  }

  theta_mu <- theta_mu_for_sigma(theta_sigma)
  list(theta_mu = theta_mu, theta_sigma = theta_sigma,
       moments = implied_moments(theta_mu, theta_sigma, n_items))
}

# ---- simulation -----------------------------------------------------------

#' Simulate one sample of vocabulary-test scores, trial level up.
#'
#' @param n            sample size.
#' @param mu           target mean observed score (0-40 scale).
#' @param reliability  target parallel-forms / test-retest reliability
#'                      (0-1); jointly with `mu` and `n_items`, this
#'                      determines the emergent between-person ability
#'                      spread (see the file header). Ignored if
#'                      `sd_true` is supplied.
#' @param sd_true      alternate calibration target: the SD (0-n_items
#'                      scale) of people's expected/"true" scores,
#'                      calibrated directly instead of `reliability` --
#'                      use this when a figure needs independent control
#'                      of spread rather than of repeatability. See
#'                      `calibrate_ability()` for the achievable range.
#' @param r_with_x     desired correlation between the *observed* score
#'                      and the covariate x, from -1 to 1 (0 = no
#'                      relationship). Capped at sqrt(reliability) if
#'                      the requested value isn't achievable, with a
#'                      warning (`reliability` here is the emergent
#'                      value when calibrating via `sd_true`).
#' @param x            optional pre-specified covariate vector (length n);
#'                      if NULL, one is drawn from N(x_mean, x_sd).
#' @param x_mean, x_sd  mean/SD used to generate x when not supplied.
#' @param n_items      number of trials/items (default 40, matching a
#'                      0-40 item test).
#' @param seed         for reproducibility.
#'
#' @return a data frame: id, x, theta (ability, logit scale), p_correct
#'   (per-item probability implied by theta), observed_score (integer,
#'   always in [0, n_items] by construction).
simulate_vocab_data <- function(n, mu = 30, reliability = 0.80, sd_true = NULL, r_with_x = 0,
                                 x = NULL, x_mean = 0, x_sd = 1,
                                 n_items = 40, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(sd_true) && !missing(reliability)) {
    warning(paste(
      "Both `reliability` and `sd_true` were supplied to",
      "simulate_vocab_data(); they aren't independently settable",
      "(fixing mu and sd_true determines reliability, and vice versa).",
      "Using `sd_true` and ignoring `reliability`."
    ))
  }
  cal <- calibrate_ability(mu, reliability = reliability, sd_true = sd_true, n_items = n_items)
  reliability <- cal$moments$reliability

  if (is.null(x)) x <- rnorm(n, x_mean, x_sd)

  max_r <- sqrt(reliability)
  if (abs(r_with_x) > max_r) {
    warning(sprintf(
      paste("r_with_x = %.2f is not achievable at reliability = %.2f",
            "(measurement error caps the observed-score correlation at",
            "about %.2f); using %.2f instead."),
      r_with_x, reliability, max_r, sign(r_with_x) * max_r
    ))
    r_with_x <- sign(r_with_x) * max_r
  }
  if (r_with_x != 0 && n < 2) {
    warning("r_with_x is not meaningful for n < 2 (correlation is undefined for a single observation); ignoring it.")
    r_with_x <- 0
  }
  r_theta_x <- if (r_with_x == 0) 0 else r_with_x / max_r

  # Only standardise x (and only use it) when it's actually needed: with
  # r_theta_x == 0, scale(x) may be NaN (e.g. n < 2, or a constant x) and
  # 0 * NaN is NaN in R, not 0 -- multiplying it in unconditionally would
  # silently turn every theta (and hence every score) into NA.
  z <- rnorm(n)
  if (r_theta_x == 0) {
    theta_z <- z
  } else {
    xz <- as.numeric(scale(x))
    if (anyNA(xz)) {
      warning("x has zero variance (e.g. all identical values); r_with_x cannot be induced. Ignoring r_with_x.")
      theta_z <- z
    } else {
      theta_z <- r_theta_x * xz + sqrt(1 - r_theta_x^2) * z
    }
  }
  theta <- cal$theta_mu + cal$theta_sigma * theta_z
  p_correct <- plogis(theta)

  observed_score <- rbinom(n, size = n_items, prob = p_correct)

  data.frame(id = seq_len(n), x = x, theta = theta, p_correct = p_correct,
             observed_score = observed_score)
}

#' Simulate `n_studies` independent samples of size `n_per_study`, for
#' the "no two samples look identical" histogram-grid figure.
#' Returns one long data frame with a `study` id column.
simulate_studies <- function(n_studies, n_per_study, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  out <- lapply(seq_len(n_studies), function(s) {
    d <- simulate_vocab_data(n_per_study, seed = NULL, ...)
    d$study <- paste("Study", s)
    d
  })
  do.call(rbind, out)
}

#' Two independent administrations (parallel forms / test-retest) per
#' person, sharing the same ability theta_i, to check/demonstrate the
#' calibrated reliability directly (corr(form1, form2) ~ reliability).
#' @param sd_true  alternate calibration target in place of `reliability`
#'                  -- see `simulate_vocab_data()` / `calibrate_ability()`.
simulate_parallel_forms <- function(n, mu = 30, reliability = 0.80, sd_true = NULL,
                                     n_items = 40, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(sd_true) && !missing(reliability)) {
    warning(paste(
      "Both `reliability` and `sd_true` were supplied to",
      "simulate_parallel_forms(); they aren't independently settable.",
      "Using `sd_true` and ignoring `reliability`."
    ))
  }
  cal <- calibrate_ability(mu, reliability = reliability, sd_true = sd_true, n_items = n_items)
  theta <- rnorm(n, cal$theta_mu, cal$theta_sigma)
  p_correct <- plogis(theta)
  data.frame(id = seq_len(n), theta = theta, p_correct = p_correct,
             form1 = rbinom(n, n_items, p_correct),
             form2 = rbinom(n, n_items, p_correct))
}

# ---- figures ----------------------------------------------------------

#' One or two scores plotted as points on a number line.
plot_dots <- function(scores, xlim = c(0, 40)) {
  # y nudged just above 0 so the dot rests on top of the axis line
  # rather than being bisected by it.
  df <- data.frame(score = scores, y = 0.04)
  ggplot(df, aes(score, y)) +
    geom_point(size = 6, color = "#3B6EA5") +
    scale_x_continuous(name = "Score (items correct)") +
    scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(xlim = xlim) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
}

#' A sample's scores as a "pile of dots" (ggdist::stat_dots) -- shows
#' the distribution as discrete points, not a smoothed density/histogram.
plot_dotpile <- function(scores, xlim = c(0, 40), ...) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plot_dotpile(). Install with install.packages('ggdist').")
  }
  df <- data.frame(score = scores)
  ggplot(df, aes(x = score)) +
    ggdist::stat_dots(color = "#3B6EA5", fill = "#3B6EA5", ...) +
    scale_x_continuous(name = "Score (items correct)") +
    coord_cartesian(xlim = xlim) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
}

#' A grid of "pile of dots" distributions, one per "study" -- the
#' ggdist::stat_dots analogue of plot_histogram_grid(), for showing
#' sampling variability as discrete points rather than binned counts.
#' @param data a data frame with a `study` and `observed_score` column,
#'   e.g. from simulate_studies().
#' @param binwidth dot diameter in data units (score points). ggdist's
#'   default computes a separate binwidth per facet panel, sized to
#'   stack that panel's tallest pile within the available height --
#'   which makes dots different sizes across panels. Fixing binwidth
#'   here keeps every panel's dots the same size.
plot_dotpile_grid <- function(data, value_col = "observed_score",
                               xlim = c(0, 40), ncol = NULL, binwidth = 2.5, ...) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plot_dotpile_grid(). Install with install.packages('ggdist').")
  }
  data$.value <- data[[value_col]]
  ggplot(data, aes(x = .value)) +
    ggdist::stat_dots(color = "#3B6EA5", fill = "#3B6EA5", binwidth = binwidth, ...) +
    scale_x_continuous(name = "Score (items correct)") +
    coord_cartesian(xlim = xlim) +
    facet_wrap(~study, ncol = ncol) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(2, "lines"),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
}

#' A histogram of one sample's scores.
plot_histogram <- function(scores, binwidth = 1, xlim = c(0, 40)) {
  df <- data.frame(score = scores)
  ggplot(df, aes(x = score)) +
    geom_histogram(binwidth = binwidth, fill = "#3B6EA5", color = "white") +
    scale_x_continuous(name = "Score (items correct)") +
    coord_cartesian(xlim = xlim) +
    ylab("Count") +
    theme_minimal()
}

#' A grid of histograms, one per "study", to show sampling variability
#' across repeated samples from the same population.
#' @param data a data frame with a `study` and `observed_score` column,
#'   e.g. from simulate_studies().
plot_histogram_grid <- function(data, value_col = "observed_score",
                                 binwidth = 1, xlim = c(0, 40), ncol = NULL) {
  data$.value <- data[[value_col]]
  ggplot(data, aes(x = .value)) +
    geom_histogram(binwidth = binwidth, fill = "#3B6EA5", color = "white") +
    scale_x_continuous(name = "Score (items correct)") +
    coord_cartesian(xlim = xlim) +
    ylab("Count") +
    facet_wrap(~study, ncol = ncol) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
}

#' Scatterplot of score against the continuous covariate, with a fitted
#' line, to illustrate the (settable) correlation.
plot_scatter <- function(data, x_col = "x", y_col = "observed_score", show_fit = TRUE) {
  df <- data.frame(x = data[[x_col]], y = data[[y_col]])
  p <- ggplot(df, aes(x, y)) +
    geom_point(color = "#3B6EA5", alpha = 0.8) +
    xlab(x_col) + ylab("Score (items correct)") +
    theme_minimal()
  if (show_fit) p <- p + geom_smooth(method = "lm", se = TRUE, color = "#B23A48")
  p
}
