import glob
import os
import polars as pl
import sys

input_dir = sys.argv[-1]

for path in glob.glob(os.path.join(input_dir, "*_filtered.tsv")):

	df = (
		pl.read_csv(path, separator="\t")
		.select(
			pl.col("Peptide").str.split_exact(by="-", n=1).struct.rename_fields(["seq1", "seq2"]),
			pl.col("# XIons").str.extract("\ (\d+-\d+)$").str.split_exact(by="-", n=1).struct.rename_fields(["rel1", "rel2"]).alias("rel"),
			pl.col("Reference").str.split_exact(by="-", n=1).struct.rename_fields(["Protein1", "Protein2"]),
			pl.col("G.Pos.1").alias("AbsPos1"),
			pl.col("G.Pos.2").alias("AbsPos2"),
			pl.col("BinoScore").alias("score"),
		)
		.unnest("Peptide", "rel", "Reference")
		.with_columns(
			pl.col("^seq\d$").str.extract(".?\.(.*)\..?").str.replace_all("[^A-Z]", ""),
			pl.col("^Protein\d$").str.strip_chars(),
		)
		.select(
			(pl.col.seq1 + "-" + pl.col.seq2 + "-a" + pl.col.rel1 + "-b" + pl.col.rel2).alias("Id"),
			pl.col.Protein1,
			pl.col.Protein2,
			pl.col.AbsPos1,
			pl.col.AbsPos2,
			pl.col.score,
		)
		#.filter(pl.all_horizontal(pl.col("^Protein\d$") == "sp|P02768|ALBU_HUMAN"))
		.filter(pl.all_horizontal(pl.col("^Protein\d$") == "sp|P00563|KCRM_RABIT"))
	)

	df.write_csv(path.split(".")[0] + "_XLA.csv")
