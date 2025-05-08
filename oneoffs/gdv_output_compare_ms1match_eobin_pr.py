import polars as pl
import sys

dfs = []
for i, pr in enumerate(sys.argv[1:]):
    dfs.append(
        pl.read_csv(pr)
        .select(
            pl.lit(i).alias("pr"),
            pl.col.Sequence,
            pl.col.GeneSymbol,
            pl.col.MZ,
            pl.col.Z,
            pl.col.NumIDSuccess,
            pl.col.NumMonitorScans,
        )
    )

df = (
    pl.concat(dfs, how="vertical")
    .pivot(
        on="pr",
        index=["Sequence", "GeneSymbol", "MZ", "Z"],
        values=["NumIDSuccess", "NumMonitorScans"],
    )
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=10000, fmt_str_lengths=1000):
    print(
        df.filter(pl.col.NumIDSuccess_0 == 0)
        .sort(pl.col.NumMonitorScans_1, pl.col.NumIDSuccess_1)
    )

    print(
        df.filter(pl.col.NumIDSuccess_0 == 1)
        .sort(pl.col.NumMonitorScans_1, pl.col.NumIDSuccess_1)
    )
