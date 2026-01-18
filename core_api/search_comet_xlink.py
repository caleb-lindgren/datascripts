import coreapipy
import datetime
import json
import os
import polars as pl
import pyopenms as oms
import re
import sys
import time

"""
Usage:
python search_comet.py COMET_WORKFLOW_PATH XLINK_WORKFLOW_PATH RAW_FILE_GROUP RAW FILE PATHS
"""

comet_workflow_path = sys.argv[1]
xlink_workflow_path = sys.argv[2]
raw_group = sys.argv[3]
raw_paths = sys.argv[4:]

timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
safe_filename_re = r"^[A-Za-z0-9_\-]+$"

if not re.match(safe_filename_re, raw_group):
	raise ValueError(f"RAW_FILE_GROUP must be only alphanumeric plus underscore and hyphen, no spaces")

uploaded = []
for raw_path in raw_paths:

	raw_filename = os.path.basename(raw_path)
	raw_filename_parts = raw_filename.split(".")

	if raw_filename_parts[1] != "raw":
		raise ValueError(f"raw file '{raw_filename}': extension must be '.raw'")

	if not re.match(safe_filename_re, raw_filename_parts[0]):
		raise ValueError(f"raw file '{raw_filename}': filename must be only alphanumeric plus underscore and hyphen, no spaces")

	raw_savename = f"{timestamp}_{raw_group}_{os.path.basename(raw_path)}"

	# Upload raw file
	resp = coreapipy.post_raw(
		path=raw_path,
		name=raw_savename,
	)

	resp.raise_for_status()
	uploaded.append(raw_savename)
	print(f"Uploaded '{raw_path}' as '{raw_savename}'")

# Get the path to the user's directory on the server
user_dir = coreapipy.get_raws_paths()[-1].split(coreapipy.username)[0] + coreapipy.username

# Queue search
last_jobs = []
protein_map_ids = []

for raw in uploaded:

	full_path = f"{user_dir}/{raw}"

	resp = coreapipy.post_search(
		workflow_path=comet_workflow_path,
		raws=[full_path],
	)

	resp_data = resp.json()
	last_jobs.append(resp_data["items"][-1]["targets"][-1]["data"]["assembler_job_id"])
	protein_map_ids.append(resp_data["items"][-1]["targets"][-1]["data"]["protein_map_id"])

# Wait for searches to finish--last job is protein_assembler.load_protein
print("Searching with Comet", end="", flush=True)
while last_jobs:
	running = []
	for job_id in last_jobs:
		status = coreapipy.get_job(job_id).json()["jobs_status"]
		if status != "done":
			running.append(job_id)
			print(".", end="", flush=True)
	last_jobs = running
	time.sleep(10)

print("\nComet search(es) finished")

# Download protein tables and create FASTA for xlink search
dfs = []
for protein_map_id in protein_map_ids:
	dfs.append(
		coreapipy.get_protein_map(protein_map_id)
		.select(pl.col.protein_id)
		.group_by("protein_id")
		.agg(pl.col.protein_id.len().alias("spectral_cts"))
	)

df = (
	pl.concat(dfs, how="vertical")
	.sort(by=["spectral_cts"], descending=True)
	.unique(subset="protein_id", keep="first")

	# Filter out any proteins that don't have at least 1% as many spectral
	# counts as the protein with the most spectral counts
	.filter(pl.col.spectral_cts >= pl.col.spectral_cts.max() * 0.01)
	.sort(by=["spectral_cts", "protein_id"], descending=True)
)

with open(comet_workflow_path, "r") as comet_workflow_file:
	comet_workflow = json.load(comet_workflow_file)

comet_params_path_server = comet_workflow["items"][3]["parameters"]["search_params"]
comet_params = coreapipy.get_search_params(
	type="comet",
	path=comet_params_path_server,
)
fasta_path_server = re.search(r"^database_name = (.*)$", comet_params, re.MULTILINE).group(1)

# TODO: download FASTA, figure out how oms wants to read it

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

# Upload, reverse, and register FASTA

# Make new xlink params with customized FASTA

# Queue xlink search and LDAs
