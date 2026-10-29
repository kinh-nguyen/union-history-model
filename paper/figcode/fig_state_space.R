# Figure: the eight states and the permitted transitions of the union history.
#
#   Rscript fig_state_space.R          (run from paper1/figcode/)
#
# Writes ../figures/state_space.png. No input data: the diagram is the state
# space itself. Boxes are drawn as rectangles in data units so the arrows can be
# clipped exactly to the box boundary; geom_label sizes its box in points, which
# leaves the arrows floating.
#
# Full journal width (17.5 cm). The chain is eight boxes wide, each carrying a
# two or three line label, so one column width is not available to it.
suppressPackageStartupMessages({ library(ggplot2); library(tibble) })
source("theme_kinh.R")   # default_font, okcol, k.theme, W_FULL, DPI

LAB_PT <- 9                     # label text on the page
HW <- 1.24; HH <- 0.55          # box half-width and half-height, data units,
                                # set to hold three lines of LAB_PT text

st <- tibble::tribble(
  ~id, ~x,  ~y,   ~lab,
  "V", 0.0,  0.0,  "V\nnever had sex",
  "X", 2.6,  1.45, "X\nsexually active,\nnever married",
  "M", 5.2,  0.0,  "M\nfirst marriage",
  "D1",7.8,  1.45, "D1\ndivorced",
  "W1",7.8, -1.45, "W1\nwidowed",
  "R", 10.4, 0.0,  "R\nlater union",
  "D2",13.0, 1.45, "D2\ndivorced",
  "W2",13.0,-1.45, "W2\nwidowed")

ed <- tibble::tribble(
  ~from,~to, ~off,
  "V","X", 0, "X","M", 0, "V","M", 0,
  "M","D1",0, "M","W1",0, "D1","R",0, "W1","R",0,
  "R","D2", 0.18, "D2","R",-0.18, "R","W2",-0.18, "W2","R", 0.18)

g <- function(k, v) st[[v]][match(k, st$id)]
edge <- function(i){
  x1 <- g(ed$from[i],"x"); y1 <- g(ed$from[i],"y")
  x2 <- g(ed$to[i],  "x"); y2 <- g(ed$to[i],  "y")
  dx <- x2-x1; dy <- y2-y1; L <- sqrt(dx^2+dy^2); ux <- dx/L; uy <- dy/L
  t <- min(if (abs(ux)>1e-9) HW/abs(ux) else Inf,
           if (abs(uy)>1e-9) HH/abs(uy) else Inf) + 0.08
  o <- ed$off[i]
  data.frame(x = x1+ux*t-uy*o, y = y1+uy*t+ux*o,
             xend = x2-ux*t-uy*o, yend = y2-uy*t+ux*o)
}
seg <- do.call(rbind, lapply(seq_len(nrow(ed)), edge))

p <- ggplot() +
  geom_segment(data = seg, aes(x,y,xend=xend,yend=yend),
               arrow = arrow(length = unit(0.20,"cm"), type="closed"),
               linewidth = 0.6, colour = "grey15") +
  geom_rect(data = st, aes(xmin=x-HW, xmax=x+HW, ymin=y-HH, ymax=y+HH),
            fill = "white", colour = "grey15", linewidth = 0.6) +
  geom_text(data = st, aes(x, y, label = lab), family = default_font,
            size = LAB_PT / .pt, lineheight = 0.95, colour = "black") +
  coord_fixed(ratio = 1, xlim = c(-1.30, 14.30), ylim = c(-2.1, 2.1),
              clip = "off") +
  theme_void(base_family = default_font)

dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
ggsave("../figures/state_space.png", plot = p, device = ragg::agg_png,
       width = W_FULL, height = W_FULL * 4.2 / 15.6 + 0.12, dpi = DPI,
       bg = "white")
cat("wrote ../figures/state_space.png\n")
