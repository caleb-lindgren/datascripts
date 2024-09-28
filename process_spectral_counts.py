import glob
import os
import polars as pl
import sys

"""
Takes peptide and protein output files from:
Core>Browse Data>View Peptide-Protein Mapping>View>Protein Map Information>Download w/Peptides | View Proteins>Download Table

Remove contaminant and decoy sequences from output, and print some statistics

Usage: python process.py BASE_DIR

BASE_DIR should contain two directories: peptides_raw and proteins_raw
"""

basedir = sys.argv[1]

msgs = {
	"contaminant_protein": [],
	"contaminant_peptide": [],
	"decoy_protein": [],
	"decoy_peptide": [],
	"protein_count": [],
	"peptide_count": [],
}

for type in ["peptides", "proteins"]:

	in_dir = os.path.join(basedir, type + "_raw")
	out_dir = os.path.join(basedir, type + "_processed")

	os.makedirs(out_dir, exist_ok=True)

	if type == "peptides":
		prot_col = "ProteinID"
	elif type == "proteins":
		prot_col = "protein_id"

	filenames = sorted(glob.glob(os.path.join(in_dir, "*.tsv")))

	# Get the maximum amount we'll need to pad the basenames by
	max_basename_len = max([len(os.path.basename(filename).split(".")[0]) for filename in filenames])

	for filename in filenames:

		basename = os.path.basename(filename).split(".")[0]

		df = (
			pl.scan_csv(filename, separator="\t")
			.with_columns(
				contaminant=pl.col(prot_col).str.ends_with("_contaminant"),
				decoy=pl.col(prot_col).str.starts_with("#")
			)
			.collect()
		)

		# Print amounts of contaminants and decoys
		contam_ct = df.filter(pl.col("contaminant")).shape[0]
		decoy_ct = df.filter(pl.col("decoy")).shape[0]

		contam_msg = f"{basename: >{max_basename_len}} contaminant {type}: {contam_ct: >4}/{df.shape[0]: <6} ({contam_ct / df.shape[0]:%})"
		decoy_msg = f"{basename: >{max_basename_len}} decoy {type}: {decoy_ct: >4}/{df.shape[0]: <6} ({decoy_ct / df.shape[0]:%})"

		if type == "peptides":
			msgs["contaminant_peptide"].append(contam_msg)
			msgs["decoy_peptide"].append(decoy_msg)
		elif type == "proteins":
			msgs["contaminant_protein"].append(contam_msg)
			msgs["decoy_protein"].append(decoy_msg)

		# Filter out contaminants and decoys
		df = (
			df.filter(
				~pl.col("contaminant"),
				~pl.col("decoy"),
			)
			.drop("contaminant", "decoy")
		)

		# Print number of unique and total peptides or proteins
		if type == "peptides":
			uniq = df.select(pl.col("peptide_sequence").unique())
			msgs["peptide_count"].append(f"{basename: >{max_basename_len}} peptides: {uniq.shape[0]: >6} unique, {df.shape[0]: >6} total")

		elif type == "proteins":
			uniq = df.select(pl.col("protein_id").unique())

			if uniq.shape[0] != df.shape[0]:
				print(f"WARNING: different number of unique ({uniq.shape[0]}) and total ({df.shape[0]}) proteins for {filename}")

			msgs["protein_count"].append(f"{basename: >{max_basename_len}} proteins: {df.shape[0]: >6}")

		df.write_csv(os.path.join(out_dir, f"{basename}.tsv"), separator="\t")

for msg_type in sorted(msgs.keys()):
	[print(msg) for msg in msgs[msg_type]]
	print()
