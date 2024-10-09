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

df_sort_order = []
df_basenames = []
df_types = []
df_stat_names = []
df_counts = []
df_decoy_counts = []
df_contaminant_counts = []
df_real_counts = []

def add_entry(sort_order, basename, type, stat_name, filter):

	df_sort_order.append(sort_order)
	df_basenames.append(basename)
	df_types.append(type)
	df_stat_names.append(stat_name)

	df_counts.append(
		filter(df)
	)
	df_decoy_counts.append(
		filter(df.filter(pl.col("decoy")))
	)
	df_contaminant_counts.append(
		filter(df.filter(pl.col("contaminant")))
	)
	df_real_counts.append(
		filter(df.filter(~pl.col("decoy") & ~pl.col("contaminant")))
	)

for type in ["peptides", "proteins"]:

	# Get input and output dirs
	in_dir = os.path.join(basedir, type + "_raw")
	out_dir = os.path.join(basedir, type + "_processed")

	os.makedirs(out_dir, exist_ok=True)

	# Set column names
	if type == "peptides":
		prot_col = "ProteinID"
	elif type == "proteins":
		prot_col = "protein_id"

	# Get all the input files
	filenames = sorted(glob.glob(os.path.join(in_dir, "*.tsv")))

	for filename in filenames:

		basename = os.path.basename(filename).split(".")[0]

		# Mark the rows that are decoys and contaminants
		df = (
			pl.scan_csv(filename, separator="\t")
			.with_columns(
				contaminant=pl.col(prot_col).str.ends_with("_contaminant"),
				decoy=pl.col(prot_col).str.starts_with("##"),
			)
			.collect()
		)

		if type == "peptides":

			# Number of successfully identified peptides (number of unique combinations of search ID and peptide ID)
			add_entry(
				sort_order=0,
				basename=basename,
				type=type,
				stat_name="successful_peptide_ids",
				filter=lambda df: df.n_unique(subset=["SearchID", "PeptideID"]),
			)

			# Number of unique peptides (number of unique values in sequence column)
			add_entry(
				sort_order=1,
				basename=basename,
				type=type,
				stat_name="unique_peptides",
				filter=lambda df: df.n_unique(subset="peptide_sequence"),
			)

			# Total PSMs
			add_entry(
				sort_order=2,
				basename=basename,
				type=type,
				stat_name="total_psms",
				filter=lambda df: df.shape[0],
			)

			# Picked PSMs
			add_entry(
				sort_order=3,
				basename=basename,
				type=type,
				stat_name="picked_psms",
				filter=lambda df: df.filter(pl.col("Peptide_parsimony").is_not_null()).shape[0],
			)

			# Unique PSMs
			add_entry(
				sort_order=4,
				basename=basename,
				type=type,
				stat_name="unique_psms",
				filter=lambda df: df.filter(pl.col("Peptide_parsimony") == "U").shape[0],
			)

			# Razor PSMs
			add_entry(
				sort_order=5,
				basename=basename,
				type=type,
				stat_name="razor_psms",
				filter=lambda df: df.filter(pl.col("Peptide_parsimony") == "R").shape[0],
			)

		elif type == "proteins":
			with pl.Config(tbl_rows=20, tbl_cols=-1, tbl_width_chars=1000, fmt_str_lengths=100):
				print(df)

			# Protein count
			add_entry(
				sort_order=6,
				basename=basename,
				type=type,
				stat_name="proteins",
				filter=lambda df: df.shape[0],
			)

			if df.n_unique("protein_id") != df.shape[0]:
				raise ValueError(f"Protein table {filename} has non-unique rows!!")

		# Filter out contaminants and decoys
		df = (
			df.filter(
				~pl.col("contaminant"),
				~pl.col("decoy"),
			)
			.drop("contaminant", "decoy")
		)

		df.write_csv(os.path.join(out_dir, f"{basename}.tsv"), separator="\t")

df_all = (
	pl.DataFrame({
		"basename": df_basenames,
		"type": df_types,
		"stat_name": df_stat_names,
		"count": df_counts,
		"decoy_count": df_decoy_counts,
		"contaminant_count": df_contaminant_counts,
		"real_count": df_real_counts,
		"sort_order": df_sort_order,
	})
	.sort("basename", "sort_order")
	.drop("sort_order")
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=1000, fmt_str_lengths=50):
	print(df_all)
