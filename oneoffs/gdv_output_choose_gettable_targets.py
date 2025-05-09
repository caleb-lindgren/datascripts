import altair as alt
import glob
import math
import os
import polars as pl
import re
import sys

from lxml import etree

# Parse the GDV target table CSVs
gdv_csvs = []
for file in [
#	"/mnt/cellbio/Gygi Lab/caleb/GoDigExperiments/MS1Matching/1200primed/ed06092_CML_MS1_match-1min_pep-30sec_2PR-10bins-noNarrow-savePRLib_60min_1200primed_GoDig2p0_SN180_600-1k_PR1_TargetTable.csv",
#	"/mnt/cellbio/Gygi Lab/caleb/GoDigExperiments/MS1Matching/1200primed/ed06094_CML_MS1_match-1min_pep-30sec_2PR-10bins-noNarrow-savePRLib_60min_1200primed_GoDig2p0_SN180_600-1k_A1_TargetTable.csv",
	"/mnt/cellbio/Gygi Lab/caleb/GoDigExperiments/MS1Matching/1200primed/ed06386_CML_MS1_match-1min_pep-30sec_1PR-narrow1st-savePRLib_60min_1285primed_SN180_600-1k_PR1_TargetTable.csv",
	"/mnt/cellbio/Gygi Lab/caleb/GoDigExperiments/MS1Matching/1200primed/ed06388_CML_MS1_match-5min_pep-30sec_1PR20250509LibOnly_60min_1285primed_SN180_600-1k_A1_TargetTable.csv",
]:
	run = re.search(r"_([PA]R?\d)_TargetTable.csv", file).group(1)

	gdv_csvs.append(
		pl.read_csv(file, schema_overrides={"TotalSumSN": pl.Float64})
		.select(
			pl.lit(run).alias("run"),
			pl.col.GeneSymbol,
			pl.col.Sequence.alias("Peptide"),
			pl.col.Z.alias("z"),
			pl.col.TotalSumSN,
			pl.col.NumPassSumSN,
			pl.col.NumIDSuccess,
		)
	)

df = (
	pl.concat(gdv_csvs, how="vertical")
	.pivot(
		on="run",
		index=["GeneSymbol", "Peptide", "z"],
		values=["TotalSumSN", "NumPassSumSN", "NumIDSuccess"],
	)
	.drop("TotalSumSN_PR1", "NumPassSumSN_PR1")
	.filter((pl.col.NumIDSuccess_PR1 > 0) & (pl.col.NumPassSumSN_A1 > 0))
	.sort("NumIDSuccess_PR1", "TotalSumSN_A1")
)

# Check overlap
v1_overlap = (
	pl.read_csv("/mnt/cellbio/Gygi Lab/caleb/GoDigMeta/TargetLists/20250509_CML_300gettable-from-1285primed.csv")
	.join(
		df.with_row_index(name="rank").select("GeneSymbol", "Peptide", "z", "rank"),
		on=["GeneSymbol", "Peptide", "z"],
		how="left"
	)
	.sort("rank")
	.with_row_index(name="rank2")
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=10000, fmt_str_lengths=1000, fmt_table_cell_list_len=100):
	print(v1_overlap)
	print(df)

df = (
	df
	.sort("NumIDSuccess_PR1", "TotalSumSN_A1")
	.tail(300)
	.select("GeneSymbol", "Peptide", "z")
)

df.write_csv("/mnt/cellbio/Gygi Lab/caleb/GoDigMeta/TargetLists/20250509_CML_300gettableV2-from-1285primed.csv", quote_style="non_numeric")
