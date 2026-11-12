# =============================================================================
# Stage 10 | DAG.R
# Directed acyclic graph of the planned PFAS -> imprinted isoform -> SRS
# mediation analysis (dagitty).
# Out: FIG_DIR/PFAS_POE_DAG.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
library(dagitty)
library(ggplot2)
library(ggrepel)
library(grid)

# --- DAG ---
g <- dagitty("
dag {
bb=\"0,0,1,1\"
X   [exposure,  pos=\"0.35,0.55\"]
Y   [pos=\"0.55,0.55\"]
SRS [outcome,   pos=\"0.85,0.55\"]
Z   [pos=\"0.15,0.60\"]
W   [pos=\"0.25,0.25\"]
U1  [unobserved,pos=\"0.48,0.40\"]
U2  [unobserved,pos=\"0.08,0.82\"]
Z  -> X
X  -> Y
Y  -> SRS
W  -> SRS
Z  -> W
U1 -> X
U1 -> Y
U2 -> Z
U2 -> Y
U2 -> SRS
U2 -> W
}
")

# Human-readable labels
labs <- c(
  X   = "X: Fetal PFAS (cord blood)",
  Y   = "Y: Placental transcription",
  SRS = "Infant SRS-2",
  Z   = "Z: Maternal PFAS (blood)",
  W   = "W: Maternal-predicted transcription",
  U1  = "U1: Unmeasured confounder (X–Y)",
  U2  = "U2: Unmeasured family/SES"
)

# --- Coordinates (robust across dagitty versions) ---
coords <- dagitty::coordinates(g)
if (is.list(coords) && all(c("x","y") %in% names(coords))) {
  node_names <- names(coords$x)
  nodes <- data.frame(name = node_names,
                      x = as.numeric(coords$x),
                      y = as.numeric(coords$y),
                      stringsAsFactors = FALSE)
} else if (is.matrix(coords) || is.data.frame(coords)) {
  node_names <- rownames(coords)
  xcol <- if ("x" %in% colnames(coords)) "x" else colnames(coords)[1]
  ycol <- if ("y" %in% colnames(coords)) "y" else colnames(coords)[2]
  nodes <- data.frame(name = node_names,
                      x = as.numeric(coords[, xcol]),
                      y = as.numeric(coords[, ycol]),
                      stringsAsFactors = FALSE)
} else stop("Unrecognized format from dagitty::coordinates(g).")

nodes$label <- ifelse(is.na(labs[nodes$name]), nodes$name, labs[nodes$name])

# --- Edge list (explicit) ---
edges_pairs <- data.frame(
  from = c("Z","X","Y","W","Z","U1","U1","U2","U2","U2","U2"),
  to   = c("X","Y","SRS","SRS","W","X","Y","Z","Y","SRS","W"),
  stringsAsFactors = FALSE
)

# Helper: recompute edges from nodes + edges_pairs
compute_edges <- function(nodes_df, pairs_df) {
  data.frame(
    x    = nodes_df$x[ match(pairs_df$from, nodes_df$name) ],
    y    = nodes_df$y[ match(pairs_df$from, nodes_df$name) ],
    xend = nodes_df$x[ match(pairs_df$to,   nodes_df$name) ],
    yend = nodes_df$y[ match(pairs_df$to,   nodes_df$name) ]
  )
}

# --- 1) SHIFT THE MATERNAL PFAS NODE (point + arrows) ---
# amount to move (10% to the right; slight down to clear incoming edges)
shift_spec <- list(
  Z = c(dx = 0.10, dy = -0.01)   # change these if you want
)

nodes_shifted <- nodes
for (nm in names(shift_spec)) {
  nodes_shifted$x[nodes_shifted$name == nm] <- nodes_shifted$x[nodes_shifted$name == nm] + shift_spec[[nm]]["dx"]
  nodes_shifted$y[nodes_shifted$name == nm] <- nodes_shifted$y[nodes_shifted$name == nm] + shift_spec[[nm]]["dy"]
}
edges_shifted <- compute_edges(nodes_shifted, edges_pairs)

# --- 2) OPTIONAL: per-node label nudges (label only; points/arrows stay) ---
label_nudge <- list(
  Z   = c(nx = +0.02, ny = 0.00),  # small extra nudge after moving the point
  X   = c(nx =  -0.01, ny = -0.02),
  Y   = c(nx =  0.00, ny = +0.02),
  U1  = c(nx =  0.00, ny = +0.02),
  U2  = c(nx =  0.00, ny = +0.02),
  W   = c(nx =  0.00, ny = 0.00),
  SRS = c(nx =  0.00, ny =  0.00)
)

labels_df <- within(nodes_shifted, {
  label_x <- x; label_y <- y
})
for (nm in names(label_nudge)) {
  labels_df$label_x[labels_df$name == nm] <- labels_df$label_x[labels_df$name == nm] + label_nudge[[nm]]["nx"]
  labels_df$label_y[labels_df$name == nm] <- labels_df$label_y[labels_df$name == nm] + label_nudge[[nm]]["ny"]
}

# --- 3) Plot: cleaner repel + padding, no axes/grid ---
g <- ggplot() +
  geom_point(
    data = nodes_shifted,
    aes(x = x, y = y),
    size = 3
  ) +
  geom_segment(
    data = edges_shifted,
    aes(x = x, y = y, xend = xend, yend = yend),
    lineend = "round",
    color="red",
    arrow = grid::arrow(length = unit(6, "pt"), type = "closed")
  ) +
  ggrepel::geom_label_repel(
    data = labels_df,
    aes(x = label_x, y = label_y, label = label),
    size = 3.6,
    label.size = 0.25,
    max.overlaps = Inf,
    segment.colour = NA,
    box.padding = 0.6,
    point.padding = 0.25,
    force = 2,           # stronger push
    force_pull = 0.2,
    min.segment.length = 0,
    seed = 42
  ) +
  coord_equal(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02)) +  # tiny margin so labels aren’t clipped
  theme_void()
g

ggsave(fig("PFAS_POE_DAG.png"), g, width = 7.5, height = 6.5, units = "in", dpi = 300, bg = "white")
