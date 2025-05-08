import polars as pl

"""
Take the given file, downloaded from Peptide View on Core, containing "Gene Symbol" and "PA Gene Symbol" columns, and
where those columns are null, fill them with the value of the "Protein ID" column. Useful for processing the table
as needed for use in GoDig library generation, when the search includes proteins that didn't have their gene symbols
properly filled in in the headers of the original FASTA the search was done with.
"""

file = "ec10224_CML_ProKAS_RC_peptide_view.csv"

df = (
	pl.read_csv(file)
	.with_columns(
		pl.col("Gene Symbol").fill_null(pl.col("Protein ID")),
		pl.col("PA Gene Symbol").fill_null(pl.col("Protein ID")),
	)
)

df.write_csv(file)
