import glob
import polars as pl

"""
Get the targets that were missed in the specified runs, based on GoDig Viewer visualization export.
"""

dfs = []
for file in glob.glob("GoDigExperiments/Shin18plex/*_TargetTable.csv"):
    dfs.append(
        pl.read_csv(file)
        .select(
            pl.col.GeneSymbol,
            pl.col.Sequence.alias("Peptide"),
            pl.col.Z.alias("z"),
            pl.col.NumPassSumSN,
            pl.lit(file.split("Caleb-Shin")[1]).alias("name"),
        )
    )

df = (
    pl.concat(dfs, how="vertical")
    .group_by(pl.col.name, pl.col.GeneSymbol, pl.col.Peptide, pl.col.z)
    .agg(pl.col.NumPassSumSN.max())
    .filter(pl.col.NumPassSumSN == 0)
    .drop("NumPassSumSN", "name")
    .unique()
    .sort("GeneSymbol", "Peptide", "z")
)

df.write_csv(f"250402_Caleb-Shin_HCC44_4runLib-missed_{df.shape[0]}-targets.csv", quote_style="non_numeric")
