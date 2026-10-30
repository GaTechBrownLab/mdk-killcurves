# =====================================================================
# MDK kill curves - ggplot2 version
#
# Reads mdk_points.csv (written by mdk_killcurves.R) and draws the
# figure. Kept deliberately layer-by-layer so it reads as a worked
# example: each geom_/scale_/theme_ block does one visible thing.
#
#   y = log10( N(t) / N0 )   so 0 = the inoculum, -2 = 99% killed,
#                            -4 = 99.99% killed. Those two dotted
#                            lines are the MDK thresholds: MDK99 is
#                            where a curve first crosses -2, MDK99.99
#                            where it first crosses -4.
# =====================================================================

library(ggplot2)

pts <- read.csv("mdk_points.csv", stringsAsFactors = FALSE)

# ---- 1. dose is ORDERED, so it gets an ordered colour scale ---------
# Dose is not a set of unrelated categories - it runs low to high. So
# it gets one hue going light -> dark, not six different hues. The
# no-drug control is not a dose at all, so it sits outside the ramp in
# neutral grey and reads as the reference it is.
dose_levels <- c(0, 0.4, 2, 4, 24, 100)
dose_cols <- c(
  "0"   = "#7a7974",   # control - neutral, deliberately not on the ramp
  "0.4" = "#86b6ef",
  "2"   = "#5598e7",
  "4"   = "#2a78d6",
  "24"  = "#184f95",
  "100" = "#0d366b"
)
dose_labs <- c("0  (no drug)", "0.4", "2", "4", "24", "100")

pts$Dose <- factor(pts$Treatment, levels = dose_levels)
pts$Detected <- factor(ifelse(pts$censored, "below detection", "counted"),
                       levels = c("counted", "below detection"))
pts <- pts[order(pts$Lineage, pts$Treatment, pts$t), ]

# ---- 2. the plot ----------------------------------------------------
p <- ggplot(pts, aes(x = t, y = log_surv,
                     colour = Dose, group = Dose)) +

  # reference lines first, so data draws on top of them
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_hline(yintercept = c(-2, -4), linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +

  geom_line(linewidth = 0.7) +

  # open triangle = the plate was blank, so the point is the detection
  # limit and the true value is somewhere below it
  geom_point(aes(shape = Detected), size = 2.1, stroke = 0.7,
             fill = "white") +

  facet_wrap(~ Lineage) +

  scale_colour_manual(values = dose_cols, labels = dose_labs,
                      drop = FALSE) +
  scale_shape_manual(values = c("counted" = 16, "below detection" = 6)) +
  scale_x_continuous(breaks = c(0, 8, 24, 48)) +
  # the right-hand axis names the two threshold lines, so the student
  # does not have to hold "-2 means 99%" in their head while reading
  scale_y_continuous(
    breaks = seq(-8, 2, 2),
    sec.axis = dup_axis(breaks = c(-2, -4),
                        labels = c("99% killed", "99.99% killed"),
                        name = NULL)) +

  labs(
    title    = "Tobramycin kill curves, evolved N. subflava lines",
    subtitle = "dotted lines: 99% and 99.99% killing - the MDK thresholds",
    x = "time (h)",
    y = expression(log[10] ~ N(t) / N[0]),
    colour = "Tobramycin\n(ug/mL)",
    shape  = NULL,
    caption = "Open triangles sit at the detection limit; true density is lower."
  ) +

  # ---- 3. recessive furniture, so the curves carry the chart --------
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
    strip.text         = element_text(face = "bold", hjust = 0),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "grey35", size = 9),
    plot.caption       = element_text(colour = "grey45", size = 8, hjust = 0),
    legend.key.height  = unit(11, "pt"),
    axis.title         = element_text(colour = "grey25")
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape  = guide_legend(order = 2, override.aes = list(colour = "grey30"))
  )

ggsave("mdk_killcurves_gg.png", p, width = 9.5, height = 3.8, dpi = 200, bg = "white")
print(p)
