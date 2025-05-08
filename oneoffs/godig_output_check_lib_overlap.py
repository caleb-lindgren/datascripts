import glob
import polars as pl

"""
Check the overlap between GoDig output files from experiments done with the same target list but two different
spectral libraries
"""

dfs = []

schema_overrides = {
	"Precursor Int (Log10)": pl.Float64,
	"MS2 IT": pl.Float64,
	"MS3 IT": pl.Float64,
	"Fraction of Fragments": pl.Float64,
	"Quant 3": pl.Float64,
}

for file in glob.glob("names_peptides_appended/*HCC44*Result.csv"):
	dfs.append(
		pl.read_csv(file, schema_overrides=schema_overrides)
		.select(
			pl.col.GeneSymbol,
			pl.lit("HCC44").alias("lib"),
		)
	)

for file in glob.glob("names_peptides_appended/*DTB*Result.csv"):
	dfs.append(
		pl.read_csv(file, schema_overrides=schema_overrides)
		.select(
			pl.col.GeneSymbol,
			pl.lit("20k").alias("lib"),
		)
	)

df = (
	pl.concat(dfs, how="vertical")
	.unique()
	.with_columns(pl.lit(1).alias("exists"))
	.pivot(
		on="lib",
		index="GeneSymbol",
		values="exists",
	)
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=10000, fmt_str_lengths=1000, fmt_table_cell_list_len=100):
	print(df)
