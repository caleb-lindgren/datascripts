import polars as pl
import sys

"""
Takes a peptide output file from:
Core>Browse Data>View Peptide-Protein Mapping>View>Protein Map Information>Download w/Peptides
"""

file = sys.argv[1]
count_type = sys.argv[2]

if count_type == "0":
	lf = (
		pl.scan_csv(file, separator="\t")
	)

elif count_type == "1":
	lf = (
		pl.scan_csv(file, separator="\t")
		.select(pl.all().n_unique())
	)

elif count_type == "2": # Check what the difference is between duplicate peptides--apparently charge state
	lf = (
		pl.scan_csv(file, separator="\t")
		.sort(by=["PeptideSequence", "ProteinId"])
		.filter(pl.col("PeptideSequence").is_duplicated())
	)

elif count_type == "3": # Avg peptides per protein
	lf = (
		pl.scan_csv(file, separator="\t")
		.select("ProteinId", "PeptideSequence")
		.unique()
		.group_by("ProteinId")
		.len()
		.select(pl.col("len").mean())
	)

elif count_type == "4": # Avg proteins is each peptide sequence assigned to
	lf = (
		pl.scan_csv(file, separator="\t")
		.select("ProteinId", "PeptideSequence")
		.unique()
		.group_by("PeptideSequence")
		.len()
		.select(pl.col("len").mean())
	)

with pl.Config(tbl_cols=20, tbl_rows=10):
	print(lf.collect())
