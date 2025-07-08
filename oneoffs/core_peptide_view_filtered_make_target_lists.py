import glob
import polars as pl
import re

# Get the list of peptides in our spectral library
sl_path = "../../Libraries/240126_SRSIII001_qy-4cell-24frac-noFAIMS-withContams/240126_SRSIII001_qy-4cell-24frac-withContam_RENAMED_BestPSM_SpecLib.xml"

with open(sl_path, "r") as sl_handle:
	speclib_str = sl_handle.read()

matches = re.findall(
	re.compile(
		r'<LibPeptide Sequence="(.+?)" GeneSymbol="(.+?)">\n\s*<AveragedSpectrum Charge="(\d)" CV="0">',
		re.MULTILINE,
	),
	speclib_str)

peps, genes, charges = zip(*matches)

sdf = (
	pl.DataFrame({
		"Peptide": peps,
		"GeneSymbol": genes,
		"z": charges,
	})
	.unique()
)

# Load all the peptides from the standards
# Unfiltered peptide tables downloaded from peptide view in Core were filtered using the following script:
# ../general/core_pepsOnly_processSpectralCounts.py
dfs = []
for filtered in sorted(glob.glob("filtered_peps/*.tsv")):
	search = re.search(r"ed[0-9]{5}", filtered).group()
	dfs.append(
		pl.read_csv(filtered, separator="\t")
		.select(
			pl.lit(search).alias("search"),
			pl.col.peptide_sequence.str.split(".").list.get(1).alias("Peptide"),
		)
		.unique()
	)

df = (
	pl.concat(dfs, how="vertical")
	.group_by("Peptide")
	.len()
	.filter(
		pl.col.Peptide.is_in(sdf.get_column("Peptide").implode()) &
		(pl.col.len >= 6) &
		(pl.col.len <= 10)
	)
	.join(
		sdf,
		on="Peptide",
		how="inner",
	)
	.sort(by=["len", "z"], descending=[True, False])
	.unique(subset="Peptide", keep="first") # Keep only one charge state per peptide. Take the most often found one.
	.select("GeneSymbol", "Peptide", "z")
)

df5k = (
	df.sample(n=5000, with_replacement=False, seed=0)
	.sort("GeneSymbol", "Peptide")
)
df5k.write_csv("250703_CML_5k-from-5cell15runs-in-qy4cellLib.csv")

df3k = (
	df5k.sample(n=3000, with_replacement=False, seed=0)
	.sort("GeneSymbol", "Peptide")
)
df3k.write_csv("250703_CML_3k-from-5cell15runs-in-qy4cellLib.csv")

df1500 = (
	df3k.sample(n=1500, with_replacement=False, seed=0)
	.sort("GeneSymbol", "Peptide")
)
df1500.write_csv("250703_CML_1500-from-5cell15runs-in-qy4cellLib.csv")

df300 = (
	df1500.sample(n=300, with_replacement=False, seed=0)
	.sort("GeneSymbol", "Peptide")
)
df300.write_csv("250703_CML_300-from-5cell15runs-in-qy4cellLib.csv")
