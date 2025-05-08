import polars as pl
import sys

files = {
	"MS1match": "3200from2000super1285primed/ed06279_CML_MS1_match-1min_pep-30sec_3PR-10bins-noNarrow-savePRLib_60min_3200_GoDig2p0_SN180_600-1k_A1_TargetTable.csv",
	"EObin": "3200from2000super1285primed/ed06296_CML_EOBin_AR-1bin_3PR-10bins-noNarrow-savePRLib-onlyTargetPrimed_60min_3200-with-super_GoDig2p0_SN180_600-1k_A1_TargetTable.csv",
}

dfs = []
for name, run in files.items():
    dfs.append(
        pl.read_csv(run)
        .select(
            pl.lit(name).alias("run"),
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
        on="run",
        index=["Sequence", "GeneSymbol", "MZ", "Z"],
        values=["NumIDSuccess", "NumMonitorScans"],
    )
    .filter((pl.col.NumMonitorScans_MS1match > 0) & (pl.col.NumMonitorScans_EObin > 0))
    .drop("NumMonitorScans_MS1match", "NumMonitorScans_EObin")
    .with_columns(
    	pl.when((pl.col.NumIDSuccess_MS1match > 0) & (pl.col.NumIDSuccess_EObin > 0)).then(pl.lit("both"))
    	.when((pl.col.NumIDSuccess_MS1match == 0) & (pl.col.NumIDSuccess_EObin == 0)).then(pl.lit("neither"))
    	.when((pl.col.NumIDSuccess_MS1match > 0) & (pl.col.NumIDSuccess_EObin == 0)).then(pl.lit("ms1match_only"))
    	.when((pl.col.NumIDSuccess_MS1match == 0) & (pl.col.NumIDSuccess_EObin > 0)).then(pl.lit("eobin_only"))
    	.alias("found_in")
    )
    .sort("found_in", "Sequence")
)

print(df.get_column("found_in").value_counts())

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=10000, fmt_str_lengths=1000):
	print(df)
