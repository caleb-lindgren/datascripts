import altair as alt
import glob
import os
import polars as pl
import sys

"""
Remove contaminant and decoy sequences from output, and make some charts or print some statistics

Usage: python process.py BASE_DIR

BASE_DIR should contain two directories: peptides_raw and proteins_raw
"""

basedir = sys.argv[1]
for type in ["peptides", "proteins"]:

	in_dir = os.path.join(basedir, type + "_raw")
	out_dir = os.path.join(basedir, type + "_processed")

	os.makedirs(out_dir, exist_ok=True)

	if type == "peptides":
		prot_col = "ProteinID"
	elif type == "proteins":
		prot_col = "protein_id"

	for filename in sorted(glob.glob(os.path.join(in_dir, "*.tsv"))):

		basename = filename.split(os.sep)[-1].split(".")[0]

		lf = (
			pl.scan_csv(filename, separator="\t")
			.filter(
				~pl.col(prot_col).str.ends_with("_contaminant") & 
				~pl.col(prot_col).str.starts_with("#")
			)
		)

		if type == "pep":
			uniq = lf.collect().select(pl.col("peptide_sequence").unique())
			print(f"{basename: >25} - {uniq.shape[0]: >8}")

		lf.collect().write_csv(os.path.join(out_dir, f"{basename}.tsv"), separator="\t")
