#' ============================================================
#' testing-sequence.R
#'
#' Extends person-grid.R with a "person at a computer" composite icon,
#' for the narrative sequence: one person doing the task -> two people
#' doing the task -> a whole sample doing the task (grid, no scores
#' shown). Same silhouette/age/disability/colour parametrization as
#' person-grid.R carries straight over, including person-grid.R's two
#' fill modes: "diversity" (sampled skin tones) or "value" (continuous
#' diverging scale driven by a numeric score/effect per icon).
#'
#' The screen is sized and positioned to face the head (vertically
#' centred on it, at a fixed horizontal offset) so the icon reads as
#' "a person looking at a screen" rather than two unrelated shapes.
#'
#' Requires person-grid.R to be sourced first (same folder).
#'
#' Quarto usage:
#'   ```{r}
#'   source("person-grid.R")
#'   source("testing-sequence.R")
#'   testing_icons(n = 1, scores = 27)$plot
#'   testing_icons(n = 2, scores = c(31, 24))$plot
#'   testing_icons(n = 30, seed = 7)$plot   # no `scores` -> sample, no labels
#'   # colour-by-score instead of sampled diversity:
#'   testing_icons(n = 6, mode = "value", value = c(-2, -1, 0, 0.5, 1, 2))$plot
#'   ```
#' ============================================================

library(ggplot2)

# ---- monitor/computer geometry -------------------------------------------

#' Screen geometry in the same local coordinate frame as person_parts()
#' (head centred at x = 0, y = 1.62), so it can be positioned relative
#' to the head and scaled/offset identically to the person. Sized and
#' vertically centred on the head so the pair reads as "looking at a
#' screen" rather than a monitor floating near the feet.
monitor_parts <- function(scale = 1, ox = 0, oy = 0,
                           center_y = 1.62, x_offset = 0.80) {
  bezel_hw <- 0.40
  bezel_hh <- 0.34
  inset <- 0.05
  neck_h <- 0.10
  neck_hw <- 0.04
  base_hw <- 0.22
  base_h <- 0.05

  bezel_top <- center_y + bezel_hh
  bezel_bottom <- center_y - bezel_hh
  neck_top <- bezel_bottom
  neck_bottom <- neck_top - neck_h
  base_top <- neck_bottom
  base_bottom <- base_top - base_h

  parts <- list(
    bezel  = data.frame(x = x_offset + c(-bezel_hw, bezel_hw, bezel_hw, -bezel_hw),
                         y = c(bezel_top, bezel_top, bezel_bottom, bezel_bottom)),
    screen = data.frame(x = x_offset + c(-bezel_hw + inset, bezel_hw - inset,
                                          bezel_hw - inset, -bezel_hw + inset),
                         y = c(bezel_top - inset, bezel_top - inset,
                               bezel_bottom + inset, bezel_bottom + inset)),
    # inverted-T stand: narrow neck down to a wider foot
    neck   = data.frame(x = x_offset + c(-neck_hw, neck_hw, neck_hw, -neck_hw),
                         y = c(neck_top, neck_top, neck_bottom, neck_bottom)),
    base   = data.frame(x = x_offset + c(-base_hw, base_hw, base_hw, -base_hw),
                         y = c(base_top, base_top, base_bottom, base_bottom))
  )
  out <- lapply(names(parts), function(nm) {
    p <- parts[[nm]]
    data.frame(part = nm, x = ox + scale * p$x, y = oy + scale * p$y)
  })
  do.call(rbind, out)
}

#' One person + screen, as a single composite polygon set with a `role`
#' column ("person" vs "monitor") so they can be filled differently.
person_at_computer <- function(kind = "A", age = "adult", disabled = FALSE,
                                scale = 1, ox = 0, oy = 0,
                                person_color = "#3B6EA5",
                                screen_color = "#CFE8FF",
                                frame_color = "#555555") {
  p <- person_parts(kind, age, disabled, scale, ox, oy)
  p$role <- "person"
  p$fill <- person_color

  m <- monitor_parts(scale = scale, ox = ox, oy = oy)
  m$role <- "monitor"
  m$fill <- ifelse(m$part == "screen", screen_color, frame_color)

  rbind(p, m)
}

# ---- main sequence function -----------------------------------------------

#' Draw n copies of "person at computer".
#'
#' @param n        how many icons (1 = single participant, 2 = a pair,
#'                 larger = a sample).
#' @param scores   optional numeric vector of length n; if supplied, each
#'                 icon gets an "x = <score>" label underneath (used for
#'                 the n = 1 and n = 2 steps). Leave NULL for the sample
#'                 grid step, where individual scores aren't shown.
#' @param ncol     icons per row; default is n itself for n <= 4, else
#'                 a compact grid.
#' @param comp_kind, comp_age, disabled_prop, seed: same composition
#'                 controls as person_grid().
#' @param mode     "diversity" (default; sampled skin-tone fill, as
#'                 before) or "value" (continuous diverging fill driven
#'                 by `value`, e.g. to colour icons by score/effect).
#' @param skin_tones hex colours to sample from when mode = "diversity".
#' @param value    numeric vector of length n; required when
#'                 mode = "value". 0 = the null/reference point (see
#'                 `diverging_colors()` in person-grid.R).
#' @return list(plot = ggplot object, data = polygon data frame)
testing_icons <- function(n, scores = NULL, ncol = NULL,
                           comp_kind = c(A = 0.45, B = 0.45, N = 0.10),
                           comp_age = c(adult = 1.0),
                           disabled_prop = 0,
                           mode = c("diversity", "value"),
                           skin_tones = c("#4A3728", "#7A5230", "#B07A4E",
                                          "#D9A066", "#F0C892", "#F5E1C8"),
                           value = NULL,
                           seed = 1, dx = 1.9, dy = 2.5) {
  mode <- match.arg(mode)
  if (!is.null(scores)) stopifnot(length(scores) == n)
  if (is.null(ncol)) ncol <- if (n <= 4) n else min(6, ceiling(sqrt(n)))
  nrow_ <- ceiling(n / ncol)

  kinds <- sample_counts(n, comp_kind, seed)
  ages <- sample_counts(n, comp_age, seed + 1)
  disabled <- sample_counts(n, c(yes = disabled_prop, no = 1 - disabled_prop), seed + 2) == "yes"

  if (mode == "diversity") {
    set.seed(seed + 3)
    colors <- sample(skin_tones, n, replace = TRUE)
  } else {
    stopifnot("`value` must be supplied and length n when mode = 'value'" =
                !is.null(value) && length(value) == n)
    colors <- diverging_colors(value)
  }

  polys <- vector("list", n)
  labels <- vector("list", n)
  for (i in seq_len(n)) {
    r <- (i - 1) %/% ncol
    c <- (i - 1) %% ncol
    ox <- c * dx
    oy <- -r * dy
    p <- person_at_computer(kinds[i], ages[i], disabled[i], 1.0, ox, oy, colors[i])
    p$icon_id <- i
    p$poly_id <- paste(i, p$part, p$role)
    polys[[i]] <- p
    if (!is.null(scores)) {
      labels[[i]] <- data.frame(x = ox + 0.2, y = oy - 0.35,
                                 label = paste0("x = ", scores[i]))
    }
  }
  df <- do.call(rbind, polys)

  plot <- ggplot(df, aes(x, y, group = poly_id, fill = fill)) +
    geom_polygon() +
    scale_fill_identity() +
    coord_equal() +
    theme_void()

  if (!is.null(scores)) {
    lab_df <- do.call(rbind, labels)
    plot <- plot + geom_text(data = lab_df, aes(x, y, label = label),
                              inherit.aes = FALSE, size = 5)
  }

  list(plot = plot, data = df)
}

#' A grid of icon-grids, one per "study" (facet_wrap(~study)) -- the icon
#' analogue of sampling-distributions.R's plot_histogram_grid()/
#' plot_dotpile_grid(), for showing "many studies, each with their own
#' independently-sampled group of people".
#'
#' @param n_studies    number of study panels.
#' @param n_per_study  number of icons (people) per study panel.
#' @param ncol         panels per row (passed to facet_wrap).
#' @param seed         base seed; study s uses seed + s so each panel's
#'                     sample is independently drawn but reproducible.
#' @param ...          passed through to testing_icons() for each panel's
#'                     composition/colour controls (comp_kind, comp_age,
#'                     disabled_prop, mode, skin_tones, value, dx, dy).
#'                     `scores` isn't supported here (no labels, as for
#'                     the single-study sample grid).
#' @return list(plot = ggplot object, data = polygon data frame)
testing_icons_studies <- function(n_studies, n_per_study, ncol = NULL,
                                   seed = 1, ...) {
  panels <- vector("list", n_studies)
  for (s in seq_len(n_studies)) {
    g <- testing_icons(n = n_per_study, seed = seed + s, ...)
    g$data$study <- paste("Study", s)
    panels[[s]] <- g$data
  }
  df <- do.call(rbind, panels)
  df$poly_id <- paste(df$study, df$poly_id)

  plot <- ggplot(df, aes(x, y, group = poly_id, fill = fill)) +
    geom_polygon() +
    scale_fill_identity() +
    coord_equal() +
    facet_wrap(~study, ncol = ncol) +
    theme_void() +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing = unit(2, "lines"))

  list(plot = plot, data = df)
}
