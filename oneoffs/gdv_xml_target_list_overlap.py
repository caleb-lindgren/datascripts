import altair as alt
import glob
import math
import os
import polars as pl
import re
import sys

from lxml import etree

# Parse the target list
targets = (
	pl.read_csv(sys.argv[1])
	.select(
		pl.col.GeneSymbol.alias("gene"),
		pl.col.Peptide.alias("seq"),
		pl.col.z,
	)
)

# Parse the GDVXMLs
gdvxmls = []
for file in sys.argv[2:]:

	priming_regex = re.search(r"PR(\d)\.gdvxml$", file)
	analytical_regex = re.search(r"A(\d)\.gdvxml$", file)

	if priming_regex:
		run_type = "priming"
		run_idx = int(priming_regex.group(1))
	elif analytical_regex:
		run_type = "analytical"
		run_idx = int(analytical_regex.group(1))
	else:
		raise ValueError(f"Bad filename: {file}")
	
	gdvxmls.append({
		"file": file,
		"type": run_type,
		"idx": run_idx,
	})

ns = {
	"xsd": "http://www.w3.org/2001/XMLSchema",
	"xsi": "http://www.w3.org/2001/XMLSchema-instance",
}

# From each run, get the retention time and target for each successful MS2
run_types = []
run_idcs = []
rts = []
ms_levels = []
genes = []
seqs = []
mzs = []
zs = []

for gdvxml in gdvxmls:

	gdvxml["tree"] = etree.parse(os.path.join(gdvxml["file"]))

	for scan in gdvxml["tree"].xpath(
		"""
		//SerialScanNumToScanInfo
		/SerializableKVP[
		    Value[@xsi:type='SerialIDMS2']
		    /CosineMatchList/IDMS2MatchAttempt/SPSIonsPassIsoFilter[
		        not(@xsi:nil='true') and number(text()) != 0
		    ]
		]/Value
		""",
		namespaces=ns
	):
		rt = float(scan.findtext("RT"))

		for target in scan.xpath("TargetsTriggeredByMasterMonitor/SimpleTarget"):
			run_types.append(gdvxml["type"])
			run_idcs.append(gdvxml["idx"])
			rts.append(rt)
			ms_levels.append(2)
			genes.append(target.findtext("GeneSymbol"))
			seqs.append(target.findtext("Sequence"))
			mzs.append(float(target.findtext("MZ")))
			zs.append(int(target.findtext("Z")))

# From the analytical run GDVXMLs, get the analytical RT for each successful MS3
for gdvxml in gdvxmls:

	if gdvxml["type"] != "analytical":
		continue

	# Get MS3 info
	for scan in gdvxml["tree"].xpath(
		"""
		//SerialScanNumToScanInfo
		/SerializableKVP[Value[@xsi:type='SerialOTMS3']]
		/Value
		""",
		namespaces=ns
	):
		rt = float(scan.findtext("RT"))

		for target in scan.xpath("TargetAtIndex"):
			run_types.append(gdvxml["type"])
			run_idcs.append(gdvxml["idx"])
			ms_levels.append(3)
			rts.append(rt)
			genes.append(target.findtext("GeneSymbol"))
			seqs.append(target.findtext("Sequence"))
			mzs.append(float(target.findtext("MZ")))
			zs.append(int(target.findtext("Z")))

scans = (
	pl.DataFrame({
		"run_type": run_types,
		"run_idx": run_idcs,
		"ms_level": ms_levels,
		"rt": rts,
		"gene": genes,
		"seq": seqs,
		"mz": mzs,
		"z": zs,
	})
	.group_by("run_type", "run_idx", "ms_level", "gene", "seq", "z")
	.len(name="num_scans")
	.with_columns(pl.format("{}_run{}_ms{}", "run_type", "run_idx", "ms_level").alias("on"))
	.pivot(
		on="on",
		index=["gene", "seq", "z"],
		values="num_scans",
	)
)

df = (
	targets.join(
		scans,
		on=["gene", "seq", "z"],
		how="full",
		coalesce=True,
	)
	.filter(pl.all_horizontal(pl.all().is_not_null()))
	.sort(
		by=["priming_run1_ms2", "analytical_run1_ms2", "analytical_run1_ms3"],
		descending=[False, True, True],
	)
	.tail(300)
	.select(
		pl.col.gene.alias("GeneSymbol"),
		pl.col.seq.alias("Peptide"),
		pl.col.z,
	)
)

df.write_csv("/mnt/clindgren/cellbio/Gygi Lab/caleb/GoDigMeta/TargetLists/20250509_CML_300gettable-from-1285primed.csv", quote_style="non_numeric")
