import polars as pl

dfs = {
	#"base": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_20240927.tsv", separator="\t"),
	#"optClipM": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_optClipM_20240927.tsv", separator="\t"),
	#"optClipM_varNAcet": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_optClipM_varNAcet_20240927.tsv", separator="\t"),
	#"varClipM": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_varClipM_20240927.tsv", separator="\t"),
	#"varClipM_varNAcet": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_varClipM_varNAcet_20240927.tsv", separator="\t"),
	#"varNAcet": pl.scan_csv("clipNMet_test/proteins_processed/az03057_2024FASTA_varNAcet_20240927.tsv", separator="\t"),
}

all_df = None
for k in dfs.keys():
	df = (
		dfs[k].select(pl.col("protein_id").alias(f"id_{k}"))
		.with_columns(pl.col(f"id_{k}").alias("id"))
	)

	if all_df is None:
		all_df = df
	else:
		all_df = all_df.join(
			df,
			on="id",
			how="full",
			coalesce=True,
		)

all_df = (
	all_df.filter(pl.any_horizontal(pl.all().is_null()))
	.sort("id")
	.collect()
)

all_df = all_df.select(sorted(all_df.columns))

#with pl.Config(tbl_rows=50, tbl_cols=-1, tbl_width_chars=1000, fmt_str_lengths=100):
#	print(all_df)

for col in all_df.columns:
	print(f"{col}: {all_df.get_column(col).is_null().sum()}")

# Output of lines 40-41 for a variety of combinations of files
# The labels may be misleading. The number next to each label is the number of proteins from the *other* file that
# weren't in this file. So a lower number means this file did better.
#
# (pyenv) data; python tmp.py
# id: 0
# id_base: 23
# id_optClipM: 19
# (pyenv) data; python tmp.py
# id: 0
# id_base: 26
# id_varClipM: 15
# (pyenv) data; python tmp.py
# id: 0
# id_base: 25
# id_varNAcet: 20
# (pyenv) data; python tmp.py
# id: 0
# id_base: 74
# id_optClipM_varNAcet: 51
# (pyenv) data; python tmp.py
# id: 0
# id_base: 81
# id_varClipM_varNAcet: 58
# (pyenv) data; python tmp.py
# id: 0
# id_optClipM: 18
# id_varClipM: 11
# (pyenv) data; python tmp.py
# id: 0
# id_optClipM_varNAcet: 27
# id_varClipM_varNAcet: 27
# (pyenv) data; python tmp.py
# id: 0
# id_optClipM: 56
# id_optClipM_varNAcet: 37
# (pyenv) data; python tmp.py
# id: 0
# id_varClipM: 56
# id_varClipM_varNAcet: 44

