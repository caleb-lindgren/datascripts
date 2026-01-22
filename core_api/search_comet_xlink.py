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
python search_comet.py COMET_WORKFLOW_PATH XLINK_WORKFLOW_PATH ORGANISM RAW_FILE_GROUP XLINKER RAW FILE PATHS
"""

comet_workflow_path = sys.argv[1]
xlink_workflow_path = sys.argv[2]
organism = sys.argv[3]
raw_group = sys.argv[4]
xlinker = sys.argv[5]
raw_paths = sys.argv[6:]

local_fastas_dir = "fastas"

datestamp = datetime.datetime.now().strftime("%Y%m%d")
timestamp = datetime.datetime.now().strftime("%H%M%S")
safe_filename_re = r"^[A-Za-z0-9_\-]+$"

if not re.match(safe_filename_re, raw_group):
	raise ValueError(f"RAW_FILE_GROUP must be only alphanumeric plus underscore and hyphen, no spaces")

xlinkers = [
	"FHD",
	"BS3",
]

if xlinker not in xlinkers:
	raise ValueError(f"Invalid xlinker '{xlinker}'")

uploaded = []
for raw_path in raw_paths:

	raw_filename = os.path.basename(raw_path)
	raw_filename_parts = raw_filename.split(".")

	if raw_filename_parts[1] != "raw":
		raise ValueError(f"raw file '{raw_filename}': extension must be '.raw'")

	if not re.match(safe_filename_re, raw_filename_parts[0]):
		raise ValueError(f"raw file '{raw_filename}': filename must be only alphanumeric plus underscore and hyphen, no spaces")

	raw_savename = f"{datestamp}_{raw_group}_{os.path.basename(raw_path)}"

	# Upload raw file
	upload_resp = coreapipy.post_raw(
		path=raw_path,
		name=raw_savename,
	)

	uploaded.append(raw_savename)
	print(f"Uploaded '{raw_path}' as '{raw_savename}'")

# Get the path to the user's directory on the server
user_dir = coreapipy.get_raws_paths()[-1].split(coreapipy.username)[0] + coreapipy.username

# Queue search
def queue_track_search(
	type,
	raw_filenames,
	user_dir,
	get_track_job_id,
	get_output_id,
	workflow=None
	workflow_path=None,
):

	if workflow is None and workflow_path is None:
		raise ValueError("Must provide either workflow or workflow path")

	last_jobs = []
	output_ids = []

	for raw in raw_file_names:

		full_path = f"{user_dir}/{raw}"

		if workflow is None:
			post_search_resp = coreapipy.post_search(
				raws=[full_path],
				workflow_path=workflow_path,
			)
		else:
			post_search_resp = coreapipy.post_search(
				raws=[full_path],
				workflow=workflow,
			)

		last_jobs.append(get_track_job_id(post_search_resp))
		output_ids.append(get_output_id(post_search_resp))

	# Wait for searches to finish--last job is protein_assembler.load_protein
	print(f"Searching with {type}", end="", flush=True)
	while last_jobs:
		running = []
		for job_id in last_jobs:
			status = coreapipy.get_job(job_id).json()["jobs_status"]
			if status != "done":
				running.append(job_id)
				print(".", end="", flush=True)
		last_jobs = running
		time.sleep(10)

	print("\nSearches finished")

	return output_ids

protein_map_ids = queue_track_search(
	workflow_path=comet_workflow_path,
	raw_filenames=uploaded,
	user_dir=user_dir,
	get_track_job_id=lambda resp: resp["items"][-1]["targets"][-1]["data"]["assembler_job_id"],
	get_output_id=lambda resp: resp["items"][-1]["targets"][-1]["data"]["protein_map_id"],
)

# Download protein tables and create FASTA for xlink search
comet_protein_map_dfs = []
for protein_map_id in protein_map_ids:
	comet_protein_map_dfs.append(
		coreapipy.get_protein_map(protein_map_id)
		.select(pl.col.protein_id)
		.group_by("protein_id")
		.agg(pl.col.protein_id.len().alias("spectral_cts"))
	)

comet_protein_map_df = (
	pl.concat(comet_protein_map_dfs, how="vertical")
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

comet_fasta_path_server = re.search(r"^database_name = (.*)$", comet_params, re.MULTILINE).group(1)
comet_fasta_path_local = os.path.join(
	local_fastas_dir,
	os.path.basename(comet_fasta_path_server),
)

def load_fasta_maybe_download(local_path, server_path):

	if not os.path.isfile(local_path):
		with open(local_path, "w") as fasta_file:
			fasta_file.write(coreapipy.get_fasta(path=server_path))

	fasta = []
	oms.FASTAFile().load(local_path, fasta)

	return fasta

comet_fasta = load_fasta_maybe_download(
	local_path=comet_fasta_path_local,
	server_path=comet_fasta_path_server,
)

xlink_fasta = []

for prot in comet_protein_map_df.get_column("protein_id"):
	added = 0
	for entry in comet_fasta:
		if entry.identifier == prot:
			xlink_fasta.append(entry)
			added += 1

	if added == 0:
		raise ValueError(f"{prot} not found")

# Save the xlink FASTA
raws_names = sorted([
	os.path.basename(raw_path).strip(".raw")
	for raw_path in raw_paths
])

xlink_fasta_name = f"{datestamp}_{timestamp}_xlink_{raw_group}_{raws_names[0]}-{raws_names[-1]}.fasta"
xlink_fasta_path_local = os.path.join(local_fastas_dir, xlink_fasta_name)

oms.FASTAFile().store(xlink_fasta_path_local, xlink_fasta)

# Upload, reverse, and register the xlink FASTA
coreapipy.post_fasta(
	name=xlink_fasta_name,
	path=xlink_fasta_path_local,
	reverse=True,
	register=True,
	organism=organism,
	type="Uniprot",
	notes=None,
	search_algorithm="Comet",
)

# Make new xlink params with the xlink FASTA
xlink_fasta_path_server = "/".join([
	os.path.dirname(comet_fasta_path_server),
	xlink_fasta_name,
])

with open(xlink_workflow_path, "r") as xlink_workflow_file:
	xlink_workflow = json.load(xlink_workflow_file)

xlink_params_template_path_server = xlink_workflow["items"][3]["parameters"]["search_params"]
xlink_params = coreapipy.get_search_params(
	type="xlink",
	path=xlink_params_template_path_server,
)

xlink_params_fasta_label = "fastaFile ="

xlink_params = re.sub(
	pattern=rf"^{xlink_params_fasta_label} *$",
	repl=f"{xlink_params_fasta_label} {xlink_fasta_path_server}",
	string=xlink_params,
	count=1,
	flags=re.MULTILINE
)

xlink_params_name = f"{datestamp}_xlink-{xlinker}_{raw_group}_{raws_names[0]}-{raws_names[-1]}.params"

xlink_params_resp = coreapipy.post_search_params(
	type="xlink",
	name=xlink_params_name,
	params_str=xlink_params,
)

xlink_params_path_server = xlink_params_resp["path"]

# Queue xlink search and LDAs
xlink_workflow["items"][3]["parameters"]["search_params"] = xlink_params_path_server

xlink_searched_ids = queue_track_search(
	workflow=xlink_workflow,
	raw_filenames=uploaded,
	user_dir=user_dir,
	get_track_job_id=lambda resp: resp["items"][-1]["targets"][-1]["data"]["job_id"],
	get_output_id=lambda resp: resp["items"][-1]["targets"][-1]["data"]["search_id"],
)

xlink_searched_ids = pl.DataFrame({
	"search_id": xlink_searched_ids,
})

xlink_searched_ids.write_csv(
	os.path.join(
		"core_download_params",
		f"{datestamp}_xlink-{xlinker}_{raw_group}_{raws_names[0]}-{raws_names[-1]}.tsv",
	),
	separator="\t",
)
