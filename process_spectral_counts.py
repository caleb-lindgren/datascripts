import altair as alt
import glob
import os
import polars as pl
import sys

"""
Takes peptide and protein output files from:
Core>Browse Data>View Peptide-Protein Mapping>View>Protein Map Information>Download w/Peptides | Download Protein List

Remove contaminant and decoy sequences from output, and print some statistics

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
		df = (
			pl.scan_csv(filename, separator="\t")
			.with_columns(
				contaminant=pl.col(prot_col).str.ends_with("_contaminant"),
				decoy=pl.col(prot_col).str.starts_with("#")
			)
		)

		#print(df.collect())

		# Print amounts of contaminants and decoys

		# Filter out contaminants and decoys
		df = (
			df.filter(
				~pl.col("contaminant"),
				~pl.col("decoy"),
			)
			.drop("contaminant", "decoy")
			.collect()
		)

		# Print number of unique and total peptides or proteins
		if type == "peptides":
			uniq = df.select(pl.col("peptide_sequence").unique())
			print(f"{basename: >25} peptides: {uniq.shape[0]: >8} unique, {df.shape[0]: >8} total")
		elif type == "proteins":
			uniq = df.select(pl.col("protein_id").unique())

			if uniq.shape[0] != df.shape[0]:
				print(f"WARNING: different number of unique ({uniq.shape[0]}) and total ({df.shape[0]}) proteins for {filename}")

			print(f"{basename: >25} proteins: {df.shape[0]: >8}")

		df.write_csv(os.path.join(out_dir, f"{basename}.tsv"), separator="\t")
