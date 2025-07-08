import glob
import os
import polars as pl
import sys

"""
Takes peptide output files from:
Core>Browse Data>View Peptide-Protein Mapping>View>Protein Map Information>Download w/Peptides

Remove contaminant and decoy sequences from output, and print some statistics

Usage: python process.py BASE_DIR

BASE_DIR should contain the peptide files
"""

basedir = sys.argv[1]

df_sort_order = []
df_basenames = []
df_stat_names = []
df_counts = []
df_decoy_counts = []
df_contaminant_counts = []
df_real_counts = []

def add_entry(sort_order, basename, stat_name, filter):

	df_sort_order.append(sort_order)
	df_basenames.append(basename)
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

# Get input and output dirs
in_dir = basedir
out_dir = os.path.join(basedir, "..", "peps_processed")

os.makedirs(out_dir, exist_ok=True)

# Set column names
prot_col = "ProteinID"

# Get all the input files
filenames = sorted(glob.glob(os.path.join(in_dir, "*.tsv")))

for filename in filenames:

	basename = os.path.basename(filename).rsplit(".", maxsplit=1)[0]

	# Mark the rows that are decoys and contaminants
	df = (
		pl.scan_csv(filename, separator="\t")
		.with_columns(
			contaminant=pl.col(prot_col).str.ends_with("_contaminant"),
			decoy=pl.col(prot_col).str.starts_with("##"),
		)
		.collect()
	)

	# Number of successfully identified peptides (number of unique combinations of search ID and peptide ID)
	add_entry(
		sort_order=0,
		basename=basename,
		stat_name="successful_peptide_ids",
		filter=lambda df: df.n_unique(subset=["SearchID", "PeptideID"]),
	)

	# Total PSMs
	add_entry(
		sort_order=1,
		basename=basename,
		stat_name="total_psms",
		filter=lambda df: df.shape[0],
	)

	# Picked PSMs
	add_entry(
		sort_order=2,
		basename=basename,
		stat_name="picked_psms",
		filter=lambda df: df.filter(pl.col("Peptide_parsimony").is_not_null()).shape[0],
	)

	# Unique PSMs
	add_entry(
		sort_order=3,
		basename=basename,
		stat_name="unique_psms",
		filter=lambda df: df.filter(pl.col("Peptide_parsimony") == "U").shape[0],
	)

	# Razor PSMs
	add_entry(
		sort_order=4,
		basename=basename,
		stat_name="razor_psms",
		filter=lambda df: df.filter(pl.col("Peptide_parsimony") == "R").shape[0],
	)

	# Number of unique peptides (number of unique values in sequence column)
	add_entry(
		sort_order=5,
		basename=basename,
		stat_name="unique_peptides",
		filter=lambda df: df.n_unique(subset="peptide_sequence"),
	)

	# Filter out contaminants and decoys
	df = (
		df.filter(
			~pl.col("contaminant"),
			~pl.col("decoy"),
		)
		.drop("contaminant", "decoy")
	)

	# Only keep unique and Razor PSMs
	df = df.filter(pl.col("Peptide_parsimony").is_not_null())

	df.write_csv(os.path.join(out_dir, f"{basename}.tsv"), separator="\t")

def format_col(coltype, max_col_width):
	return pl.concat_str(
		[
			pl.col("stat_name"),
			pl.lit(f" {coltype}: "),
			pl.col(f"{coltype}_count").map_elements(lambda ct: f"{ct: >{max_col_width},}", return_dtype=pl.String),
			pl.lit(" / "),
			pl.col("count").map_elements(lambda ct: f"{ct: <{max_col_width},}", return_dtype=pl.String),
			pl.lit(" ("),
			pl.col(f"{coltype}_prop").map_elements(lambda prop: f"{prop:.2%}", return_dtype=pl.String),
			pl.lit(")")
		],
		separator="",
	)

df_all = (
	pl.DataFrame({
		"basename": df_basenames,
		"stat_name": df_stat_names,
		"count": df_counts,
		"decoy_count": df_decoy_counts,
		"contaminant_count": df_contaminant_counts,
		"real_count": df_real_counts,
		"sort_order": df_sort_order,
	})
	.sort("basename", "sort_order")
	.drop("sort_order")
	.with_columns(
		decoy_prop=pl.col("decoy_count") / pl.col("count"),
		contaminant_prop=pl.col("contaminant_count") / pl.col("count"),
		real_prop=pl.col("real_count") / pl.col("count"),
	)
)

max_col_width = 9

df_all = (
	df_all.with_columns(
		decoy=format_col("decoy", max_col_width=max_col_width),
		contaminant=format_col("contaminant", max_col_width=max_col_width),
		real=format_col("real", max_col_width=max_col_width),
	)
	.select("basename", "real", "decoy", "contaminant")
)

for basename in sorted(df_all.get_column("basename").unique()):

	sel = (
		df_all.filter(pl.col("basename") == basename)
	)

	with pl.Config(
		tbl_rows=-1,
		tbl_cols=-1,
		tbl_width_chars=1000,
		fmt_str_lengths=1000,
		tbl_cell_alignment="RIGHT",
		tbl_hide_column_data_types=True,
		tbl_hide_dataframe_shape=True,
		tbl_hide_column_names=True,
	):
		print(sel)
