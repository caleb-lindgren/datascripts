import os
import polars as pl
import pyopenms as oms
import sys
import time

# Input is protein table(s) downloaded from protein map in Core, all will be combined
dfs = []
for res in sys.argv[1:-1]:
	dfs.append(pl.read_csv(res, separator="\t"))

df = (
	pl.concat(dfs, how="vertical")
	.select(
		pl.col.protein_id,#.replace({"sp|P02768|ALBU_HUMAN_contaminant": "sp|P02768|ALBU_HUMAN"}),
		pl.col.assigned_peptides.alias("spectral_cts"),
	)
	.sort(by=["spectral_cts"], descending=True)
	.unique(subset="protein_id", keep="first")
	.sort(by=["spectral_cts", "protein_id"], descending=True)
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=10000, fmt_str_lengths=1000, fmt_table_cell_list_len=100):
	print(df)

full_fasta_path = sys.argv[-1]
full_fasta = []
oms.FASTAFile().load(full_fasta_path, full_fasta)

xlink_fasta_path = f"from_comet_{time.strftime('%Y%m%d-%H%M%S')}.fasta"
xlink_fasta = []

for prot in df.get_column("protein_id"):
	added = 0
	for entry in full_fasta:
		if entry.identifier == prot:
			xlink_fasta.append(entry)
			added += 1

	if added == 0:
		raise ValueError(f"{prot} not found")

oms.FASTAFile().store(xlink_fasta_path, xlink_fasta)
