"""Subset the GYB cooler file to include only chromosomes 4 and 11 for the purpose of writing an R package."""

# %%
from pathlib import Path

import cooler

from utils import subset_chr_cooler


# %%
file_mcool = Path("~/projects/ONT/analysis/3c/cooler/GYB/GYB.chrs.mcool").expanduser()
resolution = 100000
clr = cooler.Cooler(f"{file_mcool}::/resolutions/{resolution}")


# %%
clr_chr4_11 = subset_chr_cooler(clr, ["chr4", "chr11"], "../inst/extdata/cooler/chr4_11-1mb.cool")


# %%
