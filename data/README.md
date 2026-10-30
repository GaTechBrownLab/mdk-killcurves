# Data file

The data are not included in this repository (see the main README). `mdk_killcurves.R` expects this file here:

**`(complete)_20261019_1058_evolvedT1-3_MDK_TOB_completed.csv`**: spot-plate colony counts, one row per time point, sample, dilution and plate. Comma-separated.
Columns: `Time` (for example 0h, 1h, 48h), `Sample ID`, `Lineage` (T1-T3), `Treatment` (tobramycin dose; 0 is the no-drug control), `Dilution`, `Volume_ul` (spot volume), `Rep1`-`Rep4` (colonies per spot; blank = not countable), `Correction_factor`, `Notes`.
