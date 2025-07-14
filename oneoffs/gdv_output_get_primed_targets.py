import glob
import polars as pl
import sys

"""
Get the targets that were primed in the specified runs, based on GoDig Viewer visualization export.
"""

dfs = []
for file in sys.argv[1:]:
    dfs.append(
        pl.read_csv(file)
        .select(
            pl.col.GeneSymbol,
            pl.col.Sequence.alias("Peptide"),
            pl.col.Z.alias("z"),
            pl.col.NumIDSuccess,
        )
    )

df = (
    pl.concat(dfs, how="vertical")
    .group_by(pl.col.GeneSymbol, pl.col.Peptide, pl.col.z)
    .agg(pl.col.NumIDSuccess.max())
    .filter(pl.col.NumIDSuccess > 0)
    .drop("NumIDSuccess")
    .unique()
    .sort("GeneSymbol", "Peptide", "z")
)

df.write_csv("primed.csv", quote_style="non_numeric")
