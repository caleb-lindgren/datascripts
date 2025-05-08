import polars as pl

"""
Combine the two listed target lists and take 3200 randomly from them
"""

dfs = []

for i, file in enumerate([
	"231107_SRSII081C_2k-MaxOverlap-RT15-110_5_Overlap15-top2000under31_GoDigTargets.csv",
	"250131_SRSIV001G_GoDig2p0_Primed_GoDigTargets_1285.csv",
]):
	dfs.append(
		pl.read_csv(file)
	)

df = (
	pl.concat(dfs, how="vertical")
	.unique(subset=pl.col.Peptide, keep="none")
	.sample(3200, with_replacement=False, shuffle=True)
)

df.write_csv("250403_CML_3200-from-2000superfliers-1285primed.csv", quote_style="non_numeric")
