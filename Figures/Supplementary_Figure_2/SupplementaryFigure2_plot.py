"""
Plot oncoplot based on prepared files
"""

# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


# %%
def read_text_list(filepath: str) -> list[str]:
    """Read text file into list of strings."""
    with open(filepath, "r", encoding="utf8") as handle:
        lines = [line.strip() for line in handle.readlines()]
    return lines


# %%
event_types = [
    "Missense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Splice_Site",
    "Translation_Start_Site",
    "Promoter",
]
event_colors = {
    "Missense_Mutation": "#2ca02c",
    "Frame_Shift_Del": "#1f77b4",
    "Frame_Shift_Ins": "#9467bd",
    "In_Frame_Del": "#bcbd22",
    "In_Frame_Ins": "#d62728",
    "Nonsense_Mutation": "#17becf",
    "Nonstop_Mutation": "#e377c2",
    "Splice_Site": "#ff7f0e",
    "Translation_Start_Site": "#8c546b",
    "Promoter": "#00008b",
    "Multi_Hit": "#000000",
    "None": "#f2f2f2",
}


# %% events in matrix form
def event_agg(events):
    events = list(events)
    if len(events) == 1:
        return events[0]
    else:
        return "Multi_Hit"


# %%
list_genes = read_text_list("genes.txt")
list_samples = read_text_list("samples.txt")

metadata = pd.read_csv("metadata.tsv", sep="\t", index_col="STUDY-LABEL")
metadata["TREATMENT"] = metadata["TREATMENT"].replace(
    {"Pall": "palliative", "Adj": "adjuvant"}
)
metadata.loc["PAT-054", "tmb"] = 9.3

# load MAF data
events_raw = pd.read_csv("dat_variants.maf", sep="\t")
events = events_raw[
    ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Protein_Change"]
]
events = events_raw

# filter by gene list
events = events[events["Hugo_Symbol"].isin(list_genes)]
# keep discarded events
events_discarded = events[~events["Variant_Classification"].isin(event_types)]
# filter by event types
events = events[events["Variant_Classification"].isin(event_types)]


events_tert = events_discarded[
    (events_discarded["Hugo_Symbol"] == "TERT")
    & (events_discarded["Variant_Classification"].isin(["5'Flank", "5'UTR"]))
]
events_tert["Variant_Classification"] = "Promoter"

# add some discarded events
events = pd.concat([events, events_tert])

# make matrix genes x samples
events_matrix = events[
    ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Protein_Change"]
].pivot_table(
    index="Hugo_Symbol",
    columns="Tumor_Sample_Barcode",
    values="Variant_Classification",
    aggfunc=event_agg,
)
# add missing samples (samples without events)
events_matrix = events_matrix.reindex(columns=metadata.index)

# %% determine genes and samples order

genes_by_count = (
    events_matrix.apply(lambda r: sum([not pd.isnull(x) for x in r]), axis=1)
    .sort_values(ascending=False)
    .index
)
events_matrix = events_matrix.reindex(index=genes_by_count)

# first, sort by events1
samples_by_events1 = (
    events_matrix.applymap(lambda x: not pd.isnull(x))
    .sort_values(events_matrix.index.tolist(), axis=1, ascending=False)
    .columns
)

# then, by metadata
metadata_sortkeys = ["TREATMENT"]
samples_by_events1_by_metadata = (
    metadata.loc[samples_by_events1].index
)

metadata = metadata.reindex(index=samples_by_events1_by_metadata)

samples1 = ~events_matrix.loc["BRAF"].isna()
samples2 = (events_matrix.loc["BRAF"].isna()) & (~events_matrix.loc["NRAS"].isna())
samples3 = (
    (events_matrix.loc["BRAF"].isna())
    & (events_matrix.loc["NRAS"].isna())
    & (~events_matrix.loc["NF1"].isna())
)
samples4 = (
    (events_matrix.loc["BRAF"].isna())
    & (events_matrix.loc["NRAS"].isna())
    & (events_matrix.loc["NF1"].isna())
)


metadata["sort"] = ""
metadata.loc[samples1, "sort"] = "1_BRAF"
metadata.loc[samples2, "sort"] = "2_NRAS"
metadata.loc[samples3, "sort"] = "3_NF1"
metadata.loc[samples4, "sort"] = "4_TripleWT"

metadata = metadata.sort_values(["TREATMENT", "sort"])

events_matrix = events_matrix.reindex(columns=metadata.index, index=list_genes)


events_matrix = events_matrix.fillna("None")

events_matrix.to_csv("dat_oncoplot.tsv", sep="\t")
events.sort_values(by=["Tumor_Sample_Barcode", "Chromosome", "Start_Position"]).to_csv(
    "dat_oncoplot_long.tsv", sep="\t", index=False
)

# %%
fig, axes = plt.subplots(
    4,
    2,
    figsize=(14, 7),
    gridspec_kw={
        "height_ratios": [4, 13, 0.8, 6],
        "width_ratios": [10, 1],
    },
    sharey="row",
)
fig.suptitle("Selected driver mutations", fontsize=20)


def plot_event_heatmap(events, ax, event_colors, labels=None):
    events_int_mapping = {v: k for k, v in enumerate(event_colors.keys())}
    events_int_repr = events.map(lambda x: events_int_mapping[x])
    events_int_repr.index.name = None
    events_int_repr.columns.name = None

    if not labels:
        labels = []

    sns.heatmap(
        events_int_repr,
        ax=ax,
        linewidths=0.5,
        cmap=list(event_colors.values()),
        cbar=False,
        xticklabels=labels,
    )
    ax.set_yticks(
        [y + 0.5 for y in range(len(events_int_repr.index))],
        events_int_repr.index.tolist(),
        fontsize=8,
        rotation=0,
    )
    ax.set_yticks([], minor=True)


def make_legend(event_colors):
    handles = []

    for level, color in event_colors.items():
        if level != "None":
            handles.append(mpatches.Patch(color=color, label=level))

    return handles


plot_event_heatmap(
    events_matrix,
    axes[1][0],
    event_colors,
)

# row counts
gene_counts = (~events_matrix.isin(["None", "no data"])).sum(axis=1).values
for y, n in enumerate(gene_counts):
    axes[1][1].text(0.1, y + 0.5, f"{n}", fontsize=8, ha="center", va="center")

axes[1][1].text(0.35, 3.5, "count", fontsize=10, rotation=270, ha="center", va="center")

axes[1][0].set_ylabel("gene")
axes[1][0].tick_params(axis="y", length=0)


def plot_annotation(metadata, factor, ax, colors, label=None):
    data = metadata[[factor]].transpose()
    levels_int_mapping = {v: k for k, v in enumerate(colors.keys())}
    data_intrep = data.applymap(lambda x: levels_int_mapping[x])

    sns.heatmap(
        data_intrep,
        ax=ax,
        linewidths=0.5,
        cmap=list(colors.values()),
        cbar=False,
        xticklabels=False,
    )
    ax.yaxis.set_tick_params(length=0)
    ax.set_yticks(
        [0.5], [data.index[0] if label is None else label], rotation=0, fontsize=8
    )
    ax.set_yticks([], minor=True)
    ax.set_xlabel("")


colors_yes_no = {"adjuvant": "red", "palliative": "blue"}
plot_annotation(
    metadata,
    "TREATMENT",
    axes[2][0],
    colors_yes_no,
    "cohort",
)
axes[2][0].set_xticks(
    [y + 0.5 for y in range(len(events_matrix.columns))],
    events_matrix.columns.tolist(),
    rotation=90,
    fontsize=8,
    ha="center",
    va="center_baseline",
)
axes[2][0].tick_params(axis="x", length=0, pad=3)


legends = [
    axes[3][0].legend(
        handles=make_legend(event_colors),
        title="Variant Classification",
        loc="upper left",
        bbox_to_anchor=(0, 0.6),
        ncol=3,
        fontsize=8,
        title_fontsize=8,
    ),
    axes[3][0].legend(
        handles=make_legend(colors_yes_no),
        title="Cohort",
        loc="upper left",
        bbox_to_anchor=(0.45, 0.6),
        ncol=1,
        fontsize=8,
        title_fontsize=8,
    ),
]

for legend in legends:
    axes[3][0].add_artist(legend)


metadata[["tmb"]].plot.bar(
    ax=axes[0][0],
    stacked=True,
    legend=False,
    color="#1f77b4",
    ylabel="TMB\n(mut/Mb)",
    xlabel="",
    xticks=[],
)

for row, col in [
    (0, 0),
]:
    for side in ["top", "bottom", "left", "right"]:
        axes[row][col].spines[side].set_visible(False)

axes[1][1].axis("off")
axes[3][0].axis("off")

for row, col in [
    (0, 1),
    (2, 1),
    (3, 1),
]:
    axes[row][col].remove()


plt.tight_layout()
plt.subplots_adjust(hspace=0.01, wspace=0.01)
plt.savefig("oncoplot.png", dpi=600, bbox_inches="tight")

