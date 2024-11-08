import os
import polars as pl
import sys

dfs = []

for path in sys.argv[1:]:

	dfs.append(
		pl.scan_csv(path, separator="\t")
		.with_columns(pl.lit(os.path.basename(path)).alias("table"))
		.select("ProteinID", "peptide_sequence", "table")
		.unique()
		.collect()
	)

if len(dfs) != 2:
	raise ValueError("我們目前只能處理兩個dataframes")

df = (
	dfs[0]
	.join(
		dfs[1],
		on=["ProteinID", "peptide_sequence"],
		how="full",
		coalesce=True,
	)
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=1000, fmt_str_lengths=100):
	print(
		df.filter(
			(
				pl.col("table").is_null()
				#pl.col("table_right").is_null()
			) & ~pl.col("peptide_sequence").is_duplicated()
		)
		.sort("peptide_sequence")
	)

#for col in ["table", "table_right"]:
#	print(df.get_column(col).is_null().sum())
