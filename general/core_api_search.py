import coreapipy
import sys

workflow_path = sys.argv[1]

path = sys.argv[1]
name = sys.argv[2]

resp = coreapipy.post_raw(
	path=path,
	name=name,
)

print(resp)
print(resp.text)

#raws = coreapipy.get_raws_paths()[-4:]
raws = coreapipy.get_raws_paths()[-1:]

for raw in raws:
	resp = coreapipy.post_search(
		workflow_path=workflow_path,
		raws=[raw],
	)

	print(f"File: {raw}\nresp: {resp}\n")
	import pprint; pprint.pprint(resp.json())
