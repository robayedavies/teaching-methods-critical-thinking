#' ============================================================
#' concept-map.R
#'
#' A general "N columns of ellipses, many-to-many curved links between
#' adjacent columns" diagram, in the style of McElreath's Statistical
#' Rethinking Figure 1.2 (theories -> process models -> statistical
#' models). Built from plain polygons + ggplot2::geom_curve, so it
#' needs nothing beyond ggplot2.
#'
#' Quarto usage (the exact figure you described):
#'
#'   ```{r}
#'   #| label: fig-theory-process-stats
#'   #| fig-cap: "Many theories and processes can imply the same statistical result"
#'   source("concept-map.R")
#'
#'   columns <- list(
#'     c("Words we experience more are easier to access",
#'       "Experience does not matter to access",
#'       "Words we experience early in life are easier to access"),
#'     c("Neural network connections change",
#'       "Neural word representations change",
#'       "Word recognition decisions change"),
#'     c("Experience correlates with word recognition speed",
#'       "Experience does not correlate with word recognition speed")
#'   )
#'   titles <- c("Theories", "Process models", "Statistical models")
#'
#'   # edges[[k]]: which item in column k links to which item in column k+1
#'   # (1-based indices, one row per link -- repeat a `from` or `to` for
#'   # many-to-many links)
#'   edges <- list(
#'     data.frame(from = c(1,1,2,3,3), to = c(1,2,3,1,2)),   # theories -> processes
#'     data.frame(from = c(1,2,3,3),   to = c(1,1,2,1))       # processes -> stats
#'   )
#'
#'   concept_map(columns, edges, titles)$plot
#'   ```
#' ============================================================

library(ggplot2)

# ---- geometry ------------------------------------------------------------

ellipse_pts <- function(cx, cy, rx, ry, n = 72) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + rx * cos(t), y = cy + ry * sin(t))
}

wrap_label <- function(label, width = 16) {
  paste(strwrap(label, width = width), collapse = "\n")
}

# ---- main function ---------------------------------------------------

#' @param columns    a list of character vectors, one per column, each
#'                   holding that column's node labels top-to-bottom.
#' @param edges      a list of length length(columns)-1; edges[[k]] is a
#'                   data frame with columns `from` and `to`: 1-based
#'                   indices into columns[[k]] and columns[[k+1]]
#'                   respectively. Repeat an index for many-to-many links.
#' @param titles     optional character vector of column headings.
#' @param rx, ry     ellipse half-width / half-height. Either a single
#'                   value (applied to every column) or a vector of
#'                   length length(columns) giving a separate size per
#'                   column -- useful when one column's labels are much
#'                   longer than another's (e.g. long theory sentences
#'                   vs. short process names).
#' @param dx, dy     horizontal spacing between columns / vertical
#'                   spacing between ellipses within a column. Note `dx`
#'                   is a single fixed gap regardless of per-column `rx`;
#'                   a warning is issued if it's too small to clear two
#'                   adjacent columns' ellipses.
#' @param wrap_width characters per line before wrapping label text.
#'                   Single value or one per column, as for `rx`/`ry`.
#' @param fill_colors  one fill colour per column; recycled/extended
#'                   with a default pastel palette if too few given.
#' @param curvature  passed to geom_curve (0 = straight lines).
#' @param arrows     if TRUE (default), draw a small arrowhead at each
#'                   curve's target end. Set FALSE for plain connecting
#'                   lines, closer to McElreath's original figure (which
#'                   shows correspondence between columns, not
#'                   directional causality).
#'
#' @return list(plot = ggplot object, nodes = node data frame,
#'   edges = edge data frame used for the curves)
concept_map <- function(columns, edges, titles = NULL,
                         rx = 1.15, ry = 0.62, dx = 3.3, dy = 1.55,
                         wrap_width = 16,
                         fill_colors = c("#F6DCC7", "#CFE3F5", "#D9E8CE", "#E9D8F0", "#F5E3B0"),
                         edge_color = "#888888", edge_alpha = 0.7, curvature = 0.2,
                         arrows = TRUE) {
  k <- length(columns)
  stopifnot(length(edges) == k - 1)
  if (length(fill_colors) < k) {
    fill_colors <- grDevices::colorRampPalette(fill_colors)(k)
  }

  # rx/ry/wrap_width can be a single value (applied to every column) or
  # a vector with one entry per column.
  recycle_per_column <- function(val, name) {
    if (length(val) == 1) return(rep(val, k))
    if (length(val) == k) return(val)
    stop(sprintf("`%s` must have length 1 or length(columns) = %d, not %d.",
                 name, k, length(val)))
  }
  rx <- recycle_per_column(rx, "rx")
  ry <- recycle_per_column(ry, "ry")
  wrap_width <- recycle_per_column(wrap_width, "wrap_width")

  # dx is a single fixed column spacing; if per-column rx makes two
  # neighbouring ellipses wider than that gap, warn rather than silently
  # rendering overlapping ellipses.
  for (ci in seq_len(k - 1)) {
    min_dx <- rx[ci] + rx[ci + 1]
    if (dx <= min_dx) {
      warning(sprintf(
        paste("dx = %.2f may be too small for columns %d and %d",
              "(rx = %.2f and %.2f): their ellipses may overlap.",
              "Consider increasing dx."),
        dx, ci, ci + 1, rx[ci], rx[ci + 1]
      ))
    }
  }

  # Validate edge indices up front, so a typo (e.g. an off-by-one, or a
  # `from`/`to` swapped) fails with a message naming the column and the
  # offending index -- rather than silently dropping the edge (as NA
  # coordinates) later and only surfacing a generic ggplot warning.
  for (ci in seq_len(k - 1)) {
    e <- edges[[ci]]
    n_from <- length(columns[[ci]])
    n_to <- length(columns[[ci + 1]])
    bad_from <- e$from[!(e$from %in% seq_len(n_from))]
    bad_to <- e$to[!(e$to %in% seq_len(n_to))]
    if (length(bad_from) > 0) {
      stop(sprintf(
        "edges[[%d]]$from has out-of-range value(s) %s: column %d ('%s') only has %d item(s), so valid indices are 1..%d.",
        ci, paste(unique(bad_from), collapse = ", "), ci,
        if (!is.null(titles)) titles[ci] else paste("column", ci), n_from, n_from
      ))
    }
    if (length(bad_to) > 0) {
      stop(sprintf(
        "edges[[%d]]$to has out-of-range value(s) %s: column %d ('%s') only has %d item(s), so valid indices are 1..%d.",
        ci, paste(unique(bad_to), collapse = ", "), ci + 1,
        if (!is.null(titles)) titles[ci + 1] else paste("column", ci + 1), n_to, n_to
      ))
    }
  }

  # node layout: column i at x = (i-1)*dx, items centred vertically,
  # using that column's own rx/ry/wrap_width.
  nodes <- do.call(rbind, lapply(seq_len(k), function(ci) {
    labels <- columns[[ci]]
    m <- length(labels)
    data.frame(
      col = ci, idx = seq_len(m),
      x = (ci - 1) * dx,
      y = (m - 1) / 2 * dy - (seq_len(m) - 1) * dy,
      label = vapply(labels, wrap_label, character(1), width = wrap_width[ci]),
      fill = fill_colors[ci],
      rx = rx[ci], ry = ry[ci],
      stringsAsFactors = FALSE
    )
  }))
  nodes$node_id <- paste(nodes$col, nodes$idx)

  node_xy <- function(ci, i) {
    row <- nodes[nodes$col == ci & nodes$idx == i, ]
    c(x = row$x, y = row$y, rx = row$rx)
  }

  # edge endpoints: right edge of source ellipse -> left edge of target,
  # each using its own node's rx (so differently-sized ellipses still
  # meet at their actual boundaries).
  edge_df <- do.call(rbind, lapply(seq_len(k - 1), function(ci) {
    e <- edges[[ci]]
    do.call(rbind, lapply(seq_len(nrow(e)), function(j) {
      p1 <- node_xy(ci, e$from[j])
      p2 <- node_xy(ci + 1, e$to[j])
      data.frame(x = p1["x"] + p1["rx"], y = p1["y"],
                 xend = p2["x"] - p2["rx"], yend = p2["y"])
    }))
  }))

  # ellipse polygons, each at its own node's rx/ry
  ellipse_df <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
    n <- nodes[i, ]
    e <- ellipse_pts(n$x, n$y, n$rx, n$ry)
    e$node_id <- n$node_id
    e$fill <- n$fill
    e
  }))

  title_df <- NULL
  if (!is.null(titles)) {
    top_y <- max(nodes$y + nodes$ry) + 0.6
    title_df <- data.frame(x = unique(nodes$x), y = top_y, label = titles)
  }

  plot <- ggplot() +
    geom_curve(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
               curvature = curvature, color = edge_color, alpha = edge_alpha,
               linewidth = 0.6,
               arrow = if (arrows) arrow(length = unit(0.12, "cm"), type = "closed") else NULL) +
    geom_polygon(data = ellipse_df, aes(x, y, group = node_id, fill = fill),
                 color = "#555555", linewidth = 0.4) +
    scale_fill_identity() +
    geom_text(data = nodes, aes(x, y, label = label), size = 3.2, lineheight = 0.9) +
    coord_equal() +
    theme_void()

  if (!is.null(title_df)) {
    plot <- plot + geom_text(data = title_df, aes(x, y, label = label),
                              fontface = "bold", size = 4.5)
  }

  list(plot = plot, nodes = nodes, edges = edge_df)
}
