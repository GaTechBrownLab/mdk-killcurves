# =====================================================================
# MDK kill curves - Neisseria subflava evolved T-lines, tobramycin
# Metrics follow Balaban et al. 2019 Nat Rev Microbiol 17:441-448
#   MDK99    = minimum duration of killing for 99%   of the population
#   MDK99.99 = minimum duration of killing for 99.99% of the population
#   tolerance   -> MDK99 elevated, MIC unchanged
#   persistence -> MDK99 normal but MDK99.99 disproportionately long
#                  (biphasic curve)
#
# figures need ggplot2
# =====================================================================
# run from the repository folder; the plate counts go in data/ (see data/README.md)

# ---- 0. input -------------------------------------------------------
d <- read.csv("data/(complete)_20261019_1058_evolvedT1-3_MDK_TOB_completed.csv", stringsAsFactors = FALSE, check.names = FALSE)

names(d) <- trimws(names(d))

# ---- transcription fixes --------------------------------------------
# Applied in code, never to the source workbook, so the raw sheet stays
# the record and every change here is printed when the script runs.
# Add a row per fix, with the date the plate was re-read and what it said.
# None outstanding. 2026-10-20: plates re-read against the sheet for
# 24 h S13/S16 and all 48 h high-dose wells - the sheet was correct
# throughout, including 24 h S16 Rep4 = 6.
# To add one, follow this shape:
#   CORRECTIONS <- data.frame(Time = "24h", SampleID = "S16", Dilution = 0,
#                             Rep = "Rep4", value = 16,
#                             reason = "plate re-read <date>: ...",
#                             stringsAsFactors = FALSE)
CORRECTIONS <- data.frame(
  Time = character(), SampleID = character(), Dilution = numeric(),
  Rep = character(), value = numeric(), reason = character(),
  stringsAsFactors = FALSE
)

# ---- 1. settings ----------------------------------------------------
SPOT_REPS      <- c("Rep1", "Rep2", "Rep3", "Rep4")
COUNT_MIN      <- 3     # countable window per spot
COUNT_MAX      <- 100
DROP_4H_S01_12 <- TRUE  # 4 h dilutions were wrong for S01-S12 (see notes)
DROP_4H_ALL    <- FALSE  # S13-S18 also misbehave at 4 h - drop the timepoint
APPLY_CORR     <- TRUE # use Correction_factor instead of dropping 4 h
N0_FROM        <- c(0, 0.4, 2, 4)  # doses used to define the t=0 inoculum
# --- what counts as a real zero -------------------------------------
# A missing replicate is NOT a zero. Blanks in Rep1-4 happen because a
# spot was uncountable, smeared, or mispipetted, and those rows are
# averaged over whatever replicates ARE present.
# A real zero is only a plate that was inspected and found empty.
BLANK_IS_ZERO  <- TRUE  # honour the verified blanks below
VERIFIED_ZERO_TIMES <- c(48)  # timepoints inspected and confirmed empty.
                              # Rows outside these times are never zeroed,
                              # whatever the Notes column says.
RECORDED_ZERO_IS_ZERO <- TRUE # a replicate written down as "0" is an
                              # observation; FALSE drops those rows instead
BLANK_SPOTS    <- 4     # spots plated on an all-empty row. Applies both to
                        # inspected-blank plates and to rows where a single
                        # "0" was written to mean "no colonies in any
                        # replicate" - in both cases every spot was examined,
                        # so the pooled detection limit uses all of them.
POOL_SPOTS_LOD <- TRUE  # LOD from total volume examined, not one spot
# Extra neat volume plated at 48 h on the separate confirmation plate, in uL.
# Counts toward the volume examined, so a blank there lowers the 48 h LOD.
CONFIRM_48H_UL <- 10

# ---- 2. tidy --------------------------------------------------------
# merged cells in the sheet leave blanks - carry the last value down
fill_down <- function(x) {
  x[x == ""] <- NA
  idx <- cumsum(!is.na(x)); idx[idx == 0] <- NA
  c(NA, x[!is.na(x)])[idx + 1]
}
for (v in c("Time", "Sample ID", "Lineage", "Treatment"))
  d[[v]] <- fill_down(as.character(d[[v]]))

d$Treatment <- as.numeric(d$Treatment)
d$Dilution  <- as.numeric(d$Dilution)
d$Volume_ul <- as.numeric(d$Volume_ul)
d$t         <- as.numeric(sub("h$", "", d$Time))
for (v in SPOT_REPS) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))

# apply the transcription fixes now that the columns are filled and numeric
if (nrow(CORRECTIONS)) {
  cat("--- transcription corrections applied ---\n")
  for (i in seq_len(nrow(CORRECTIONS))) {
    cx <- CORRECTIONS[i, ]
    k <- which(d$Time == cx$Time & d$`Sample ID` == cx$SampleID &
               !is.na(d$Dilution) & d$Dilution == cx$Dilution)
    if (length(k) != 1) {
      warning(sprintf("correction %d matched %d rows - not applied", i, length(k)))
      next
    }
    old <- d[[cx$Rep]][k]
    d[[cx$Rep]][k] <- cx$value
    cat(sprintf("  %s %s dil %g %s: %s -> %g   (%s)\n", cx$Time, cx$SampleID,
                cx$Dilution, cx$Rep, ifelse(is.na(old), "NA", old),
                cx$value, cx$reason))
  }
}

# Dilution is recorded as a FRACTION (1e-3 = 1000-fold).
# A recorded 0 means neat/undiluted -> factor 1, not zero.
d$dil <- ifelse(is.na(d$Dilution) | d$Dilution == 0, 1, d$Dilution)

d$cf <- suppressWarnings(as.numeric(d$Correction_factor))
d$cf[is.na(d$cf)] <- 1
if (!APPLY_CORR) d$cf <- 1

# An inspected-empty plate is an observation (density below detection).
# An absent replicate is a gap and is simply not averaged.
if (!"Notes" %in% names(d)) d$Notes <- NA_character_
n_present <- rowSums(!is.na(d[SPOT_REPS]))

d$blank <- BLANK_IS_ZERO &
  grepl("No colony counts recorded", d$Notes, ignore.case = TRUE) &
  n_present == 0 &
  d$t %in% VERIFIED_ZERO_TIMES

# rows with no counts that were NOT verified empty (uncountable, TNTC,
# pipetting slips) carry no information and are dropped, never zeroed
d$unusable <- n_present == 0 & !d$blank

keep <- !is.na(d$Volume_ul) & (n_present > 0 | d$blank)
d <- d[keep, ]

d$n_spots <- rowSums(!is.na(d[SPOT_REPS]))
# na.rm = TRUE: the mean is over the replicates actually read, so a row
# with 2 of 4 spots is a 2-spot mean, not a 4-spot mean with two zeros
d$mean_ct <- rowMeans(d[SPOT_REPS], na.rm = TRUE)
d$mean_ct[d$blank] <- 0

# A single "0" on the sheet means no colonies in ANY replicate - the
# other cells were left blank rather than written out. So such a row
# had every spot examined, like an inspected-blank plate, and both get
# the full spot count for the detection limit.
d$recorded_zero <- !d$blank & d$mean_ct == 0
if (!RECORDED_ZERO_IS_ZERO) {
  d <- d[!d$recorded_zero, ]
} else {
  d$n_spots[d$recorded_zero] <- BLANK_SPOTS
}
d$n_spots[d$blank] <- BLANK_SPOTS

d$cfu_ml  <- d$mean_ct / (d$Volume_ul / 1000) / d$dil * d$cf

# limit of detection. Per spot it is 1 colony / spot volume; pooling the
# spots actually examined lowers it by n_spots. (A Poisson 95% upper bound
# on a zero count is ~3x this - quoted in the console output.)
vol_examined <- (d$Volume_ul / 1000) * (if (POOL_SPOTS_LOD) pmax(d$n_spots, 1) else 1)
# the 48 h confirmation plate was neat, so it only helps undiluted samples
extra <- ifelse(d$blank & d$t == 48 & d$dil == 1, CONFIRM_48H_UL / 1000, 0)
vol_examined <- vol_examined + extra
d$lod_ml <- 1 / vol_examined / d$dil

if (DROP_4H_ALL) {
  d <- d[d$t != 4, ]
} else if (DROP_4H_S01_12) {
  bad <- d$t == 4 & d$`Sample ID` %in% sprintf("S%02d", 1:12)
  d   <- d[!bad, ]
}

# ---- 3. one density per lineage x dose x time -----------------------
# pick the dilution whose spots fall in the countable window;
# if none qualifies take the one closest to the middle of it.
pick_dilution <- function(g) {
  ok <- g[g$mean_ct >= COUNT_MIN & g$mean_ct <= COUNT_MAX, ]
  if (nrow(ok)) return(ok[which.min(abs(ok$mean_ct - (COUNT_MIN + COUNT_MAX) / 2)), ])
  nz <- g[g$mean_ct > 0, ]
  if (nrow(nz)) return(nz[which.min(abs(nz$mean_ct - (COUNT_MIN + COUNT_MAX) / 2)), ])
  g[which.min(g$lod_ml), ]            # all blank: keep the most sensitive plate
}

key <- interaction(d$Lineage, d$Treatment, d$t, drop = TRUE)
pts <- do.call(rbind, lapply(split(d, key), pick_dilution))
pts <- pts[order(pts$Lineage, pts$Treatment, pts$t), ]
rownames(pts) <- NULL

# censor zero counts at the detection limit and flag them
pts$censored <- pts$mean_ct == 0
pts$cfu_ml[pts$censored] <- pts$lod_ml[pts$censored]
pts$log10_cfu <- log10(pts$cfu_ml)

# ---- audit: every row treated as a zero, so it can be checked by eye
cat("\n--- treated as real zeros (censored at the detection limit) ---\n")
aud <- pts[pts$censored, c("Time", "Sample ID", "Lineage", "Treatment",
                           "Dilution", "Volume_ul", "n_spots", "blank",
                           "recorded_zero", "lod_ml")]
if (nrow(aud)) print(aud, row.names = FALSE) else cat("none\n")
cat("blank = plate inspected and empty;",
    "recorded_zero = a replicate written down as 0\n")
cat("\nreplicate coverage of the rows that WERE counted:\n")
print(table(spots_read = pts$n_spots[!pts$censored]))
cat("(rows with fewer than 4 spots are averaged over the spots read,",
    "not padded with zeros)\n")

# ---- 4. survival relative to the inoculum ---------------------------
# every dose of a lineage was split from one culture, so N0 is a property
# of the lineage, not of the well. Using the well's own t=0 would hide
# any killing that happened before the first sample was taken.
n0_tab <- sapply(split(pts, pts$Lineage), function(g) {
  z <- g[g$t == 0 & g$Treatment %in% N0_FROM, ]
  if (!nrow(z)) NA_real_ else median(z$log10_cfu, na.rm = TRUE)
})

pts$log10_N0  <- n0_tab[as.character(pts$Lineage)]
pts$log_surv  <- pts$log10_cfu - pts$log10_N0   # log10( N(t) / N0 )

# ---- 5. MDK by interpolation ----------------------------------------
# First time the survival curve crosses -2 (MDK99) or -4 (MDK99.99).
# Linear interpolation in log-survival between the bracketing samples.
# Returns the crossing time (t = NA if the drop is never reached) and
# whether it was reached via a censored (below-detection) point. A
# censored crossing is an UPPER bound: the true density is below the
# plotted LOD, so the real crossing is at or before it.
mdk <- function(time, logS, cens, target) {
  o <- order(time); time <- time[o]; logS <- logS[o]; cens <- cens[o]
  ok <- is.finite(time) & is.finite(logS)
  time <- time[ok]; logS <- logS[ok]; cens <- cens[ok]
  if (length(time) < 2) return(list(t = NA_real_, cens = FALSE))
  for (i in 2:length(time)) {
    if (logS[i] <= target) {
      if (logS[i - 1] <= target)
        return(list(t = time[i - 1], cens = cens[i - 1]))
      f <- (logS[i - 1] - target) / (logS[i - 1] - logS[i])
      return(list(t = time[i - 1] + f * (time[i] - time[i - 1]), cens = cens[i]))
    }
  }
  list(t = NA_real_, cens = FALSE)
}

grp <- split(pts, interaction(pts$Lineage, pts$Treatment, drop = TRUE))
res <- do.call(rbind, lapply(grp, function(g) {
  m2 <- mdk(g$t, g$log_surv, g$censored, -2)
  m4 <- mdk(g$t, g$log_surv, g$censored, -4)
  data.frame(
    Lineage   = g$Lineage[1],
    Dose      = g$Treatment[1],
    n_times   = nrow(g),
    t_max     = max(g$t),
    min_logS  = round(min(g$log_surv, na.rm = TRUE), 2),
    MDK99     = round(m2$t, 2), MDK99_cens   = m2$cens,
    MDK9999   = round(m4$t, 2), MDK9999_cens = m4$cens,
    n_blank   = sum(g$censored),
    stringsAsFactors = FALSE)
}))
res <- res[order(res$Lineage, res$Dose), ]
rownames(res) <- NULL

# ">" = threshold never reached inside the window
# "<=" = reached only via a blank plate, so the true crossing is earlier
fmt <- function(x, cens, tmax)
  ifelse(is.na(x), paste0("> ", tmax, " h"),
         ifelse(cens, paste0("<= ", x), format(x)))
out <- data.frame(res[c("Lineage", "Dose", "n_times", "t_max", "min_logS")],
                  MDK99   = fmt(res$MDK99,   res$MDK99_cens,   res$t_max),
                  MDK9999 = fmt(res$MDK9999, res$MDK9999_cens, res$t_max),
                  n_blank = res$n_blank, stringsAsFactors = FALSE)

cat("\n--- inoculum (log10 CFU/mL, median of doses",
    paste(N0_FROM, collapse = "/"), "at t=0) ---\n")
print(round(n0_tab, 2))
cat("\n--- MDK summary (h) ---\n")
print(out, row.names = FALSE)

#write data
write.csv(pts, "mdk_points.csv",  row.names = FALSE)
write.csv(out, "mdk_summary.csv", row.names = FALSE)

# ---- 6. figure ------------------------------------------------------
# mdk_plot_ggplot.R draws the kill curves from mdk_points.csv.
