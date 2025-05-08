import glob
import polars as pl

"""
Append peptides to gene names in GoDig output to prevent aggregation at gene level by Shuken's GoDig analysis R
scripts
"""

for result in sorted(glob.glob("*_Result.csv")):

	res = (
		pl.read_csv(result)
		.with_columns(pl.concat_str([pl.col.GeneSymbol, pl.col.Peptide], separator="_").alias("GeneSymbol"))
	)

	res.write_csv(f"names_peptides_appended/{result}")
