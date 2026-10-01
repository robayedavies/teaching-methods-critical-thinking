#' ============================================================
#' peer-review-workflow.R
#'
#' Conceptual flow diagram of the peer-reviewed publication cycle,
#' built from the "person at a computer" icon (testing-sequence.R)
#' reused across four roles: author, editor, reviewer x2.
#'
#' The flow snakes left to right:
#'   1. Author submits manuscript to editor.
#'   2. Editor sends manuscript out to two reviewers (fan-out).
#'   3. Reviewers submit reviews back to editor (fan-in).
#'   4. Editor makes a decision and informs the author.
#' A dashed loop-back arrow underneath represents "reject with
#' revise & resubmit -> repeat" until accept/reject.
#'
#' Requires person-grid.R and testing-sequence.R sourced first
#' (same folder).
#'
#' Quarto usage:
#'   ```{r}
#'   source("person-grid.R")
#'   source("testing-sequence.R")
#'   source("peer-review-workflow.R")
#'   peer_review_workflow()$plot
#'   ```
#' ============================================================

library(ggplot2)

#' Build the peer-review workflow diagram.
#' @return list(plot = ggplot object)
peer_review_workflow <- function() {

  researcher_col <- list(person = "#3B6EA5", screen = "#CFE8FF")
  editor_col     <- list(person = "#B0562D", screen = "#FFDCC5")
  reviewer_col   <- list(person = "#3B9E5F", screen = "#D3F5DF")

  pos <- list(
    researcher1 = c(0,    0),
    editor1     = c(5.0,  0),
    reviewerA   = c(10.5, 3.4),
    reviewerB   = c(10.5, -3.4),
    editor2     = c(16.0, 0),
    researcher2 = c(21.0, 0)
  )

  icon_df <- function(name, kind, cols) {
    p <- pos[[name]]
    d <- person_at_computer(kind = kind, ox = p[1], oy = p[2],
                             person_color = cols$person, screen_color = cols$screen)
    d$icon <- name
    d
  }

  icons <- rbind(
    icon_df("researcher1", "A", researcher_col),
    icon_df("editor1",     "B", editor_col),
    icon_df("reviewerA",   "A", reviewer_col),
    icon_df("reviewerB",   "B", reviewer_col),
    icon_df("editor2",     "B", editor_col),
    icon_df("researcher2", "A", researcher_col)
  )
  icons$poly_id <- paste(icons$icon, icons$part, icons$role)

  right_of <- function(name, dy = 0.9) c(pos[[name]][1] + 1.25, pos[[name]][2] + dy)
  left_of  <- function(name, dy = 0.9) c(pos[[name]][1] - 0.55, pos[[name]][2] + dy)

  straight_arrows <- rbind(
    data.frame(rbind(c(right_of("researcher1"), left_of("editor1")))),
    data.frame(rbind(c(right_of("editor2"),     left_of("researcher2"))))
  )
  names(straight_arrows) <- c("x", "y", "xend", "yend")

  curve_arrows <- data.frame(rbind(
    c(right_of("editor1")[1],   right_of("editor1")[2],   left_of("reviewerA")[1], left_of("reviewerA")[2],  0.30),
    c(right_of("editor1")[1],   right_of("editor1")[2],   left_of("reviewerB")[1], left_of("reviewerB")[2], -0.30),
    c(right_of("reviewerA")[1], right_of("reviewerA")[2], left_of("editor2")[1],   left_of("editor2")[2],   -0.30),
    c(right_of("reviewerB")[1], right_of("reviewerB")[2], left_of("editor2")[1],   left_of("editor2")[2],    0.30)
  ))
  names(curve_arrows) <- c("x", "y", "xend", "yend", "curvature")

  loop_arrow <- data.frame(
    x = pos$researcher2[1] + 0.35, y = -7.0,
    xend = pos$researcher1[1] + 0.35, yend = -7.0
  )

  step_captions <- data.frame(
    x = c(pos$researcher1[1] + 0.35, pos$editor1[1] + 0.35,
          (pos$reviewerA[1] + pos$reviewerB[1]) / 2, pos$editor2[1] + 0.35),
    y = -6.0,
    label = c("1. Submit\nmanuscript",
              "2. Editor sends\nto review",
              "3. Reviewers\nsubmit reviews",
              "4. Decision: accept /\nR&R / reject")
  )

  role_labels <- data.frame(
    x = c(pos$researcher1[1] + 0.35, pos$editor1[1] + 0.35,
          pos$reviewerA[1] + 0.35, pos$reviewerB[1] + 0.35,
          pos$editor2[1] + 0.35, pos$researcher2[1] + 0.35),
    # Reviewer A's label needs the same +2.4 clearance from icon centre
    # used for the on-centerline roles; the previous +1.4 offset landed
    # on the icon itself. Reviewer B keeps its original offset since its
    # label already clears the icon without crowding the step captions
    # below it.
    y = c(2.4, 2.4, pos$reviewerA[2] + 2.4, -4.8, 2.4, 2.4),
    label = c("Author", "Editor", "Reviewer 1", "Reviewer 2", "Editor", "Author")
  )

  plot <- ggplot() +
    geom_polygon(data = icons, aes(x, y, group = poly_id, fill = fill)) +
    geom_segment(data = straight_arrows,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                 linewidth = 0.6, color = "grey30")

  for (i in seq_len(nrow(curve_arrows))) {
    plot <- plot + geom_curve(data = curve_arrows[i, ],
                               aes(x = x, y = y, xend = xend, yend = yend),
                               arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                               linewidth = 0.6, color = "grey30",
                               curvature = curve_arrows$curvature[i])
  }

  plot <- plot +
    geom_curve(data = loop_arrow,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
               linewidth = 0.5, color = "grey55", linetype = "dashed",
               curvature = -0.15) +
    annotate("text", x = (pos$researcher1[1] + pos$researcher2[1]) / 2 + 0.35,
             y = -7.9, label = "reject with revise & resubmit \u2192 repeat the cycle",
             size = 3.3, color = "grey40", fontface = "italic") +
    geom_text(data = role_labels, aes(x, y, label = label), size = 3.1, color = "grey20") +
    geom_text(data = step_captions, aes(x, y, label = label), size = 3.4,
              lineheight = 0.95, fontface = "bold") +
    scale_fill_identity() +
    coord_equal(clip = "off", ylim = c(-8.8, 5.8)) +
    theme_void() +
    theme(plot.margin = margin(20, 25, 20, 25))

  list(plot = plot)
}
