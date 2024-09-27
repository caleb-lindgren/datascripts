import polars as pl
import sys

f = sys.argv[1]
t = sys.argv[2]

if t == "0":
	lf = (
		pl.scan_csv(f, separator="\t")
	)

elif t == "1":
	lf = (
		pl.scan_csv(f, separator="\t")
		.select(pl.all().n_unique())
	)

elif t == "2": # Check what the difference is between duplicate peptides--apparently charge state
	lf = (
		pl.scan_csv(f, separator="\t")
		.sort(by=["PeptideSequence", "ProteinId"])
		.filter(pl.col("PeptideSequence").is_duplicated())
	)

elif t == "3": # Avg peptides per protein
	lf = (
		pl.scan_csv(f, separator="\t")
		.select("ProteinId", "PeptideSequence")
		.unique()
		.group_by("ProteinId")
		.len()
		.select(pl.col("len").mean())
	)

elif t == "4": # Avg proteins is each peptide sequence assigned to
	lf = (
		pl.scan_csv(f, separator="\t")
		.select("ProteinId", "PeptideSequence")
		.unique()
		.group_by("PeptideSequence")
		.len()
		.select(pl.col("len").mean())
	)

with pl.Config(tbl_cols=20, tbl_rows=10):
	print(lf.collect())
