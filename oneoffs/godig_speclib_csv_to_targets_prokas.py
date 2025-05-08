import polars as pl

"""
Take spectral library CSV export from GoDig and create a target list from the ProKAS peptides that were found
"""

df = (
	pl.read_csv("prokas_CKS6_SpecLib_export.csv")
	.select(
		pl.col.gene.alias("GeneSymbol"),
		pl.col.peptide.alias("Peptide"),
		pl.col.prec_z.alias("z"),
	)
	.filter(
		pl.col.GeneSymbol.str.contains("CKS6_AGA") |
		(pl.col.GeneSymbol == "ProKAS_CKS6_Linker") |
		(pl.col.GeneSymbol == "EGFP") |
		(pl.col.GeneSymbol == "NbALFA") |
		(pl.col.GeneSymbol == "NLS")
	)
	.unique()
	.sort("GeneSymbol", "Peptide", "z")
)

df.write_csv("20250425_prokas_CKS6_targets.csv", quote_style="non_numeric")
