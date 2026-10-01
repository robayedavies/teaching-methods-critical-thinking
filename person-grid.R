#' ============================================================
#' person-grid.R
#'
#' A parametrized "restroom-sign style" pictogram generator, built with
#' plain polygons in ggplot2 (no external icon fonts, so it renders
#' identically everywhere Quarto builds the book).
#'
#' Each icon is assembled from a few flat shapes:
#'   - head          (circle)
#'   - torso         (trapezoid)
#'   - legs          ("A" = trouser silhouette, "B" = skirt silhouette,
#'                     "N" = half trouser / half skirt -- an inclusive,
#'                     non-binary-coded silhouette)
#'   - wheelchair    (swapped in for the legs when disabled = TRUE)
#'   - cane          (added when age = "elder")
#'   - overall scale (reduced when age = "child")
#'
#' Colour is decided by ONE of two modes, chosen by whether you pass
#' `value=`:
#'   mode = "diversity" (default): each icon's fill is sampled from a
#'          discrete set of skin tones -- for showing sample diversity.
#'   mode = "value": each icon's fill is placed on a continuous,
#'          symmetric diverging scale (dark negative -> light null ->
#'          dark positive) driven by a numeric `value` per icon -- for
#'          showing an effect size, a score, a probability, etc.
#'
#' Use in a Quarto R chunk:
#'
#'   ```{r}
#'   #| label: fig-diversity-sample
#'   #| fig-cap: "A hypothetical study sample"
#'   #| fig-alt: !expr diversity_grid$alt_text
#'   source("person-grid.R")
#'   diversity_grid <- person_grid(n = 30, seed = 7)
#'   diversity_grid$plot
#'   ```
#' ============================================================

library(ggplot2)

# ---- low-level geometry -------------------------------------------------

circle_pts <- function(cx, cy, r, n = 24) {
  t <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}

#' Build the list of sub-polygons for one icon, in local coordinates,
#' then offset by (ox, oy) and multiplied by `scale`.
person_parts <- function(kind = c("A", "B", "N"),
                          age = c("adult", "child", "elder"),
                          disabled = FALSE,
                          scale = 1, ox = 0, oy = 0) {
  kind <- match.arg(kind)
  age <- match.arg(age)
  parts <- list()

  # head
  parts[["head"]] <- circle_pts(0, 1.62, 0.30)

  if (disabled) {
    # Redesigned for legibility at small sizes: one large, unbroken wheel
    # (no small second circle to be mistaken for a stray dot), seated
    # directly under the torso with no gap, plus a single low footrest
    # peg. Fewer, bigger shapes read better than several small ones once
    # icons are shrunk down in a grid.
    parts[["torso"]] <- data.frame(
      x = c(-0.30, 0.30, 0.30, -0.30),
      y = c(1.28, 1.28, 0.62, 0.62)
    )
    parts[["wheel"]] <- circle_pts(0.00, 0.28, 0.34, n = 32)
    parts[["footrest"]] <- data.frame(
      x = c(0.34, 0.50, 0.50, 0.34),
      y = c(0.02, 0.02, 0.12, 0.12)
    )
  } else {
    parts[["torso"]] <- data.frame(
      x = c(-0.34, 0.34, 0.20, -0.20),
      y = c(1.30, 1.30, 0.72, 0.72)
    )

    if (kind == "A") {           # trouser legs: two bars with a slim gap
      parts[["leg_l"]] <- data.frame(x = c(-0.20, -0.03, -0.06, -0.20),
                                      y = c(0.72, 0.72, 0.00, 0.00))
      parts[["leg_r"]] <- data.frame(x = c(0.03, 0.20, 0.20, 0.06),
                                      y = c(0.72, 0.72, 0.00, 0.00))
    } else if (kind == "B") {    # skirt
      parts[["legs"]] <- data.frame(x = c(-0.20, 0.20, 0.46, -0.46),
                                     y = c(0.72, 0.72, 0.00, 0.00))
    } else if (kind == "N") {    # half trouser (left) / half skirt (right)
      parts[["leg_l"]] <- data.frame(x = c(-0.20, -0.02, -0.05, -0.20),
                                      y = c(0.72, 0.72, 0.00, 0.00))
      parts[["legs_skirt_half"]] <- data.frame(x = c(0.00, 0.20, 0.46, 0.00),
                                                y = c(0.72, 0.72, 0.00, 0.00))
    }

    if (age == "elder") {
      parts[["cane"]] <- data.frame(x = c(0.36, 0.40, 0.58, 0.54),
                                     y = c(0.78, 0.78, 0.02, 0.02))
    }
  }

  # apply scale (child = smaller) then offset into the grid position
  out <- lapply(names(parts), function(nm) {
    p <- parts[[nm]]
    data.frame(part = nm, x = ox + scale * p$x, y = oy + scale * p$y)
  })
  do.call(rbind, out)
}

# ---- composition sampling -----------------------------------------------

#' Largest-remainder sampling: turns a named vector of proportions into
#' exactly `n` labels, hitting the target proportions as closely as
#' integer counts allow, then shuffles their order.
sample_counts <- function(n, props, seed) {
  stopifnot(
    "`props` must have non-negative values, with at least one positive" =
      all(props >= 0) && sum(props) > 0
  )
  props <- props / sum(props)
  set.seed(seed)
  raw <- props * n
  counts <- floor(raw)
  remainder <- n - sum(counts)
  if (remainder > 0) {
    order <- order(-(raw - counts))
    for (i in seq_len(remainder)) {
      idx <- order[((i - 1) %% length(order)) + 1]
      counts[idx] <- counts[idx] + 1
    }
  }
  sample(rep(names(props), counts))
}

# ---- colour -----------------------------------------------------------

#' Map a numeric value to a symmetric diverging colour, vivid at both
#' extremes and light at zero (the "null effect" point). `value` need
#' not be centred on the data -- 0 is always treated as the reference
#' point, so effect sizes / differences map intuitively.
#'
#' Default endpoints are blue/orange from the Okabe-Ito colour-blind-safe
#' palette (distinguishable under deuteranopia, protanopia, and
#' tritanopia) rather than a red/blue or red/green diverging scale.
diverging_colors <- function(value, neg = "#0072B2", mid = "#F2F1EC", pos = "#E69F00") {
  span <- max(abs(value), 1e-9)
  tt <- value / span                     # -1 .. 1, 0 = null
  ramp_lo <- grDevices::colorRamp(c(neg, mid))
  ramp_hi <- grDevices::colorRamp(c(mid, pos))
  to_hex <- function(m) grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  vapply(tt, function(t) {
    if (t <= 0) to_hex(ramp_lo(t + 1))
    else        to_hex(ramp_hi(t))
  }, character(1))
}

# ---- automatic grid shape -----------------------------------------------

#' Pick the `ncol` (icons per row) that makes an n-icon grid's bounding
#' box (`ncol * dx` wide by `ceiling(n / ncol) * dy` tall) match a target
#' width:height `aspect` ratio as closely as possible -- e.g. so the
#' figure fills a 16:9 slide with minimal leftover margin either way.
choose_ncol <- function(n, dx, dy, aspect) {
  ncols <- seq_len(n)
  nrows <- ceiling(n / ncols)
  grid_ratio <- (ncols * dx) / (nrows * dy)
  # log-ratio distance treats "twice too wide" and "twice too tall"
  # as equally bad, so it doesn't systematically favour wide grids.
  ncols[which.min(abs(log(grid_ratio / aspect)))]
}

# ---- main function ------------------------------------------------------

#' Build a grid of person pictograms.
#'
#' @param n           number of icons.
#' @param ncol        icons per row (rows are computed from n). Leave as
#'                    `"auto"` (the default) to pick the layout that best
#'                    fills a window of the given `aspect` ratio, or pass
#'                    an integer to fix it yourself.
#' @param aspect      target width:height ratio used only when
#'                    `ncol = "auto"`. Defaults to 16/9. Set this to match
#'                    your actual slide size (e.g. `revealjs`'s default
#'                    is 960/700; a custom `width`/`height` in the format
#'                    YAML should be reflected here) so the grid fills
#'                    the slide with as little leftover margin as
#'                    possible.
#' @param comp_kind   named proportions for silhouette: "A" (trouser-
#'                    coded), "B" (skirt-coded), "N" (half/half --
#'                    reads as a non-binary / trans-inclusive icon).
#' @param comp_age    named proportions for "adult", "child", "elder".
#' @param disabled_prop proportion of icons shown using the wheelchair
#'                    silhouette.
#' @param mode        "diversity" (discrete skin-tone fill, sampled) or
#'                    "value" (continuous diverging fill driven by
#'                    `value`).
#' @param skin_tones  hex colours to sample from when mode = "diversity".
#' @param value       numeric vector of length n; required when
#'                    mode = "value". 0 = the null/reference point.
#' @param seed        for reproducibility.
#'
#' @return a list with `plot` (a ggplot object), `data` (the full
#'   polygon data frame, if you want to customise further), `meta`
#'   (per-icon kind/age/disabled/fill), and `alt_text` (a plain-language
#'   summary of the composition, generated from the same sampling used
#'   to draw the icons -- pass it straight into `fig-alt`).
person_grid <- function(n = 30, ncol = "auto", aspect = 16 / 9,
                         comp_kind = c(A = 0.45, B = 0.45, N = 0.10),
                         comp_age = c(adult = 0.70, child = 0.15, elder = 0.15),
                         disabled_prop = 0.12,
                         mode = c("diversity", "value"),
                         skin_tones = c("#4A3728", "#7A5230", "#B07A4E",
                                        "#D9A066", "#F0C892", "#F5E1C8"),
                         value = NULL,
                         seed = 1,
                         dx = 1.15, dy = 2.15) {
  mode <- match.arg(mode)
  if (identical(ncol, "auto")) ncol <- choose_ncol(n, dx, dy, aspect)
  nrow_ <- ceiling(n / ncol)

  kinds <- sample_counts(n, comp_kind, seed)
  ages <- sample_counts(n, comp_age, seed + 1)
  disabled <- sample_counts(n, c(yes = disabled_prop, no = 1 - disabled_prop), seed + 2) == "yes"

  if (mode == "diversity") {
    set.seed(seed + 3)
    fill <- sample(skin_tones, n, replace = TRUE)
  } else {
    stopifnot("`value` must be supplied and length n when mode = 'value'" =
                !is.null(value) && length(value) == n)
    fill <- diverging_colors(value)
  }

  polys <- vector("list", n)
  for (i in seq_len(n)) {
    r <- (i - 1) %/% ncol
    c <- (i - 1) %% ncol
    ox <- c * dx
    oy <- -r * dy
    scale <- if (ages[i] == "child") 0.65 else 1.0
    p <- person_parts(kinds[i], ages[i], disabled[i], scale, ox, oy)
    p$icon_id <- i
    p$poly_id <- paste(i, p$part)
    p$fill <- fill[i]
    polys[[i]] <- p
  }
  df <- do.call(rbind, polys)

  # No expansion/margins: the panel should be exactly the icons' bounding
  # box, so that when Reveal stretches the image to fill the slide there
  # is no leftover white border eating into the "fill the window" effect.
  plot <- ggplot(df, aes(x, y, group = poly_id, fill = fill)) +
    geom_polygon() +
    scale_fill_identity() +
    coord_equal(expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  meta <- data.frame(icon_id = seq_len(n), kind = kinds, age = ages,
                      disabled = disabled, fill = fill)

  alt_text <- if (mode == "diversity") {
    sprintf(
      "Grid of %d person icons of varying skin tone, gender presentation, and age, illustrating sample diversity. %d%% shown as children, %d%% as older adults, and %d%% using a wheelchair icon.",
      n,
      round(100 * mean(ages == "child")),
      round(100 * mean(ages == "elder")),
      round(100 * mean(disabled))
    )
  } else {
    sprintf(
      "Grid of %d person icons coloured on a diverging scale from %.2f to %.2f, blue for negative values, near-white at zero, orange for positive values.",
      n, min(value), max(value)
    )
  }

  list(plot = plot, data = df, meta = meta, alt_text = alt_text)
}
