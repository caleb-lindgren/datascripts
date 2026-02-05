import altair as alt
import numpy as np
import os
import polars as pl
import sys

"""
Usage:

python make_bca_curve.py LABELS_FILE_PATH VALUES_FILE_PATH
"""

def fit(labels_file, values_file):

	# Read in labels and values
	rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
	labels = (
		pl.read_csv(labels_file, has_header=False, new_columns=[str(i) for i in range(1, 13)], separator="\t")
		.with_columns(pl.Series(name="row", values=rows))
	)
	values = (
		pl.read_csv(values_file, has_header=False, new_columns=[str(i) for i in range(1, 13)], separator="\t")
		.with_columns(pl.Series(name="row", values=rows))
	)

	# Melt both identically and join
	labels = labels.unpivot(index="row", variable_name="column", value_name="label")
	values = values.unpivot(index="row", variable_name="column")

	df = (
		labels.join(
			values,
			on=["column", "row"],
			how="full",
			coalesce=True,
		)
		.filter(pl.col("label") != "")
		.select("label", pl.col("value").cast(pl.Float64), pl.concat_str("row", "column").alias("well"))
	)

	# Mean blank replicates and subtract blank
	df_means = df.group_by("label").mean()
	if df_means.filter(pl.col("label") == "standard_0").get_column("value").shape[0] > 0:
		blank_mean_absorbance = df_means.filter(pl.col("label") == "standard_0").get_column("value")[0]
		df = df.with_columns(pl.col("value") - blank_mean_absorbance)

	# Separate out standards and samples
	standards_regex = pl.col("label").str.contains("standard_\d{1,4}", literal=False)
	stands = df.filter(standards_regex)
	samples = df.filter(~standards_regex)

	# Expand sample labels
	samples = (
		samples.with_columns(
			pl.col("label").str.split("_1/").list.to_struct(fields=["sample", "dilution_factor"], upper_bound=2)
		)
		.unnest("label")
		.with_columns(pl.col("dilution_factor").cast(pl.Int64))
	)

	# Expand standards labels
	stands = (
		stands.with_columns(pl.col("label").str.split("standard_").list.to_struct(fields=["_", "conc"], upper_bound=2))
		.unnest("label")
		.drop("_")
		.with_columns(pl.col("conc").cast(pl.Float64))
		.sort(by="conc")
	)

	# Fit equation to standards
	coef = np.polynomial.polynomial.polyfit(
		x=stands.get_column("value"),
		y=stands.get_column("conc"),
		deg=2,
	)

	# Predict for samples
	samples = samples.with_columns(
		conc_pred=np.polynomial.polynomial.polyval(
			x=samples.get_column("value"),
			c=coef,
		)
	)

	# Plot
	pldf = (
		pl.concat([
			stands.with_columns(pl.concat_str(
				pl.lit("standard_"),
				pl.col("conc").cast(pl.Int32).cast(pl.String).str.zfill(4),
			).alias("group")),
			samples.select(pl.col("conc_pred").alias("conc"), "value", "well", pl.col("sample").alias("group"))
		])
		.sort("group")
	)

	x = np.linspace(pldf.get_column("value").min(), pldf.get_column("value").max(), 1000)
	linedf = pl.DataFrame({
		"x": x,
		"y": np.polynomial.polynomial.polyval(x=x, c=coef),
	})

	samples_chart = alt.Chart(pldf.filter(~pl.col("group").str.starts_with("standard_"))).mark_point().encode(
		x="value",
		y="conc",
		color=alt.Color(
			"group",
			sort=pldf.get_column("group").to_list(),
		),
		tooltip=["group", "well"],
	)

	standards_base = alt.Chart(pldf.filter(pl.col("group").str.starts_with("standard_"))).encode(
		x="value",
		y="conc",
		color=alt.Color(
			"group",
			sort=pldf.get_column("group").to_list(),
		),
		tooltip=["group", "well"],
	)

	standards_points = standards_base.mark_point(shape="square", filled=True, opacity=1)
	standards_lines = standards_base.mark_line()

	line = alt.Chart(linedf).mark_line().encode(
		x=alt.X(
			"x",
			axis=alt.Axis(title="Absorbance"),
		),
		y=alt.Y(
			"y",
			axis=alt.Axis(title="Concentration (ug/mL)"),
		),
	)

	chart = (
		alt.layer(line, standards_lines, standards_points, samples_chart)
		.properties(
			width=800,
			height=800,
		)
		.interactive()
	)

	labels_basename = os.path.basename(labels_file).split(".")[0]
	values_basename = os.path.basename(values_file).split(".")[0]
	chart.save(f"charts/{labels_basename}_{values_basename}.html")

	# Scale and combine dilutions
	samples = (
		samples.with_columns((pl.col("dilution_factor") * pl.col("conc_pred")).alias("full_conc_pred"))
		.select(
			pl.col("sample"),#.str.split("_").list.slice(0, 2).list.join("_"),
			"well",
			pl.col("full_conc_pred").alias("conc"),
			"dilution_factor",
		)
		.with_columns(pl.col("conc").mean().over("sample").alias("group_mean"))
		.sort("sample")
	)

	# Save results
	samples.write_csv(f"concs/{labels_basename}_{values_basename}.tsv", separator="\t")

fit(labels_file=sys.argv[1], values_file=sys.argv[2])
