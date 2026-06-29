#!/usr/bin/env python3
"""Interactive UMAP+HDBSCAN scatter of miRNA binding profiles.

Reads the ``results`` table from ``rimap_results.duckdb``, aggregates per
miRNA (mean binding_dG, Kd, bridge, etc.), projects to 2-D via UMAP, clusters
with HDBSCAN, and writes a self-contained Plotly HTML at ``viz/index.html``.

Note: this script runs its pipeline at import-time rather than under a ``main()``
guard. Run it directly (`python scripts/mirna_interactive_plot.py`).
"""
# venv: .venv/bin/python3
import duckdb
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from umap import UMAP
from hdbscan import HDBSCAN
import plotly.express as px
import plotly.graph_objects as go

con = duckdb.connect("rimap_results.duckdb", read_only=True)
df = con.execute("SELECT * FROM results").fetchdf()
con.close()

df["log_Kd"] = np.log1p(df["Kd"])

num_agg = (
    df.groupby("mirna_id")
    .agg(
        mean_binding_dG=("binding_dG", "mean"),
        mean_log_Kd=("log_Kd", "mean"),
        mean_bridge=("bridge", "mean"),
        std_binding_dG=("binding_dG", "std"),
        std_log_Kd=("log_Kd", "std"),
        median_binding_dG=("binding_dG", "median"),
        n_interactions=("binding_dG", "count"),
    )
    .fillna(0)
)

top_seed_types = df["seed_type"].value_counts().head(12).index.tolist()
df_seed = df[df["seed_type"].isin(top_seed_types)]
seed_props = (
    df_seed.groupby(["mirna_id", "seed_type"])
    .size()
    .unstack(fill_value=0)
)
seed_props = seed_props.div(seed_props.sum(axis=1), axis=0)
seed_props.columns = [f"seed_{c}" for c in seed_props.columns]

trans_props = (
    df.groupby(["mirna_id", "transcript_id"])
    .size()
    .unstack(fill_value=0)
)
trans_props = trans_props.div(trans_props.sum(axis=1), axis=0)
trans_props.columns = [f"trans_{c}" for c in trans_props.columns]

features = num_agg.join(seed_props, how="left").join(trans_props, how="left").fillna(0)
mirna_ids = features.index.tolist()

X_scaled = StandardScaler().fit_transform(features.values)
X_umap = UMAP(n_components=2, n_neighbors=30, min_dist=0.1, metric="euclidean", random_state=42).fit_transform(X_scaled)
labels = HDBSCAN(min_cluster_size=15, min_samples=5, metric="euclidean", cluster_selection_method="eom").fit_predict(X_umap)

n_clusters = len(set(labels) - {-1})
n_noise = (labels == -1).sum()

mirna_cluster = pd.DataFrame({"mirna_id": mirna_ids, "cluster": labels})
df_clustered = df.merge(mirna_cluster, on="mirna_id")

top_transcripts_per_cluster = {}
for c in sorted(set(labels) - {-1}):
    cluster_df = df_clustered[df_clustered["cluster"] == c]
    top5 = cluster_df["transcript_id"].value_counts().head(5)
    top_transcripts_per_cluster[c] = "<br>".join(
        f"  {tid} ({cnt:,} int.)" for tid, cnt in top5.items()
    )

seed_profile_per_cluster = {}
for c in sorted(set(labels) - {-1}):
    mask = labels == c
    seed_cols = [col for col in features.columns if col.startswith("seed_")]
    top_seeds = features.loc[np.array(mirna_ids)[mask], seed_cols].mean().sort_values(ascending=False).head(4)
    seed_profile_per_cluster[c] = "<br>".join(
        f"  {col.replace('seed_', '')}: {val:.1%}" for col, val in top_seeds.items()
    )

plot_df = pd.DataFrame({
    "mirna_id": mirna_ids,
    "UMAP_1": X_umap[:, 0],
    "UMAP_2": X_umap[:, 1],
    "cluster": labels,
    "cluster_label": [f"Cluster {l}" if l >= 0 else "Noise" for l in labels],
    "mean_binding_dG": features["mean_binding_dG"].values,
    "mean_log_Kd": features["mean_log_Kd"].values,
    "mean_bridge": features["mean_bridge"].values,
    "n_interactions": features["n_interactions"].values,
})

hover_parts = [
    "<b>%{customdata[0]}</b><br>",
    "Cluster: %{customdata[1]}<br>",
    "<br><b>Features:</b><br>",
    "  binding_dG (mean): %{customdata[2]:.2f}<br>",
    "  log Kd (mean): %{customdata[3]:.2f}<br>",
    "  bridge (mean): %{customdata[4]:.2f}<br>",
    "  interactions: %{customdata[5]:,}<br>",
    "<br><b>Top targets:</b><br>%{customdata[6]}<br>",
    "<br><b>Seed profile:</b><br>%{customdata[7]}",
]
hovertemplate = "".join(hover_parts) + "<extra></extra>"

customdata = np.column_stack([
    plot_df["mirna_id"].values,
    plot_df["cluster_label"].values,
    plot_df["mean_binding_dG"].values,
    plot_df["mean_log_Kd"].values,
    plot_df["mean_bridge"].values,
    plot_df["n_interactions"].values,
    [top_transcripts_per_cluster.get(c, "N/A") for c in labels],
    [seed_profile_per_cluster.get(c, "N/A") for c in labels],
])

fig = go.Figure()

noise_mask = labels == -1
if noise_mask.any():
    fig.add_trace(go.Scattergl(
        x=plot_df.loc[noise_mask, "UMAP_1"],
        y=plot_df.loc[noise_mask, "UMAP_2"],
        mode="markers",
        marker=dict(size=5, color="lightgrey", opacity=0.4),
        name=f"Noise ({noise_mask.sum()})",
        customdata=customdata[noise_mask],
        hovertemplate=hovertemplate,
    ))

palette = px.colors.qualitative.Alphabet + px.colors.qualitative.Dark24
for i, c in enumerate(sorted(set(labels) - {-1})):
    mask = labels == c
    color = palette[i % len(palette)]
    fig.add_trace(go.Scattergl(
        x=plot_df.loc[mask, "UMAP_1"],
        y=plot_df.loc[mask, "UMAP_2"],
        mode="markers",
        marker=dict(size=7, color=color, opacity=0.75, line=dict(width=0.3, color="white")),
        name=f"Cluster {c} ({mask.sum()})",
        customdata=customdata[mask],
        hovertemplate=hovertemplate,
    ))

fig.update_layout(
    title=dict(
        text=f"miRNA Clustering — UMAP + HDBSCAN<br>"
             f"<sup>{n_clusters} clusters, {n_noise} noise points, {len(mirna_ids)} miRNAs</sup>",
        font=dict(size=18),
    ),
    xaxis_title="UMAP 1",
    yaxis_title="UMAP 2",
    template="plotly_white",
    width=1400,
    height=900,
    legend_title_text="Cluster",
    legend=dict(itemsizing="constant", font=dict(size=10)),
    hoverlabel=dict(font_size=12, bgcolor="white"),
    dragmode="pan",
)

out_path = "mirna_clusters_interactive.html"
fig.write_html(out_path, include_plotlyjs="cdn")
print(f"Saved: {out_path}")
print(f"Clusters: {n_clusters}, Noise: {n_noise}/{len(labels)}")
