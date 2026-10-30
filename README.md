# mdk-killcurves

Tobramycin kill curves and minimum duration for killing (MDK) for evolved lines of *Neisseria subflava*.

## Question

Do the evolved lines survive tobramycin longer, and if so, is it tolerance (the whole population dies more slowly) or persistence (a small subpopulation survives)? Following Balaban et al. (2019, *Nat Rev Microbiol* 17:441-448), the scripts estimate MDK99, the time to kill 99% of the population, and MDK99.99, the time to kill 99.99%, for three evolved lineages (T1-T3) across six tobramycin doses, including a no-drug control, over 48 h. A raised MDK99 points to tolerance; a normal MDK99 with a disproportionately long MDK99.99 (a biphasic curve) points to persistence.

## My role

I wrote and ran the kill-curve analysis.

## Status

Ongoing work; not yet published.

## What the scripts do

`mdk_killcurves.R` turns the spot-plate colony counts into one density per lineage, dose and time point. It handles countable windows, missing versus truly empty spots, and pooled detection limits, and prints every correction it applies. It then calculates survival relative to the inoculum and the MDK99 and MDK99.99 for each lineage and dose, and writes `mdk_points.csv` and `mdk_summary.csv`.

`mdk_plot_ggplot.R` draws the kill-curve figure from `mdk_points.csv`, layer by layer (`mdk_killcurves_gg.png`).

The analysis settings (countable range, how blanks and recorded zeros are treated, which doses define the inoculum) are listed and explained at the top of `mdk_killcurves.R`.

## How to run

Last run in October 2026 with R 4.6.0 and ggplot2 4.0.3. `mdk_killcurves.R` needs only base R; `mdk_plot_ggplot.R` needs ggplot2.

1. Put the plate-count sheet in `data/` (see Data availability). `data/README.md` describes its columns.
2. Set the working directory to the repository folder, run `mdk_killcurves.R`, then `mdk_plot_ggplot.R`.

## Data availability

The colony counts are not included while the work is unpublished. They are available from me on request.

## License

The code is released under the MIT License (see `LICENSE`).
