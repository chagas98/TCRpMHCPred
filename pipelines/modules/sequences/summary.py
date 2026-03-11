import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple, List
from collections import defaultdict
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# number of samples
# unique peptides
# unique cdr3a/cdr3b
# unique mhcs
# 10x

# - plot
# - pep position diversity
# - pep length distribution
# - score distribution

# in-depth analysis
# tcrdist3 -> TCR
# levenshtein distance -> peptide
# -> compate with test set, between datasets, sobreposition

# test
# all train summaries
# sobreposition with training data (%)

# peptide length distribution

class SummarizeData:
    def __init__(self):
        self.summary = {}
        self.datasets = {}
        self.labels = {}

    def peplen(self, df: pd.Series) -> dict:
        counts = df.astype(str).apply(len).value_counts()
        proportions = (counts / counts.sum()).to_dict()
        return proportions, counts.to_dict()

    def mhc_distribution(self, df: pd.Series, top_alleles: list) -> tuple:
        mhc_data = df.astype(str).apply(lambda x: ':'.join(x.split(':')[:2]) if x.count(':') >= 2 else x)
        mhc_data = mhc_data.apply(lambda x: x if x in top_alleles else 'Others')
        counts = mhc_data.value_counts()
        proportions = (counts / counts.sum()).to_dict()
        return proportions, counts.to_dict()

    def add_dataset(self, df: pd.DataFrame, name: str, label: str):
        if label not in {"train", "test"}:
            raise ValueError("Label must be 'train' or 'test'")
        self.datasets[name] = df
        self.labels[name] = label

    def run(self):

        all_data = pd.concat(self.datasets.values(), ignore_index=True)
        mhc_data = all_data['MHCa'].astype(str).apply(lambda x: ':'.join(x.split(':')[:2]) if x.count(':') >= 2 else x)
        top10_alleles = mhc_data.value_counts().nlargest(10).index.tolist()
        
        for name, df in self.datasets.items():
            mhc_prop, mhc_counts = self.mhc_distribution(df['MHCa'], top10_alleles)
            peplen_prop, peplen_counts = self.peplen(df['epitope'])

            stats = {
                "set_type": self.labels[name],
                "dataset_name": name,
                "n_samples": len(df),
                "n_epitopes": df["epitope"].nunique(),
                "n_cdr3a": df["CDR3A"].nunique() if "CDR3A" in df.columns else None,
                "n_cdr3b": df["CDR3B"].nunique() if "CDR3B" in df.columns else None,
                "n_mhc_alleles": df["MHCa"].nunique(),
                "peptide_length_counts": peplen_counts,
                "peptide_length_prop": peplen_prop,
                "mhc_distribution_prop": mhc_prop,
                "mhc_distribution_counts": mhc_counts,
                "score_distribution": df['Score'].value_counts().to_dict() if self.labels[name] == 'train' else None
            }

            self.summary[name] = stats

    def to_df(self) -> pd.DataFrame:
        """
        Convert the summary to a DataFrame.
        """
        return pd.DataFrame([stats for stats in self.summary.values()])

    def save_table(self, path: str = 'summary.csv'):
        """
        Save the summary DataFrame to a CSV file.
        """
        self.to_df().to_csv(path, index=False)
        print(f"Summary saved to {path}")

    def plot_summary_stats(self, save_dir: str = '.'):
        """
        Plot key statistics for all datasets in the summary.
        """
        df = self.to_df()
        stats_to_plot = [
            ("n_samples", "Number of Samples", 'Number of Samples'),
            ("n_epitopes", "Unique Epitopes", "Epitopes"),
            ("n_cdr3a", "Unique CDR3A", "CDR3A"),
            ("n_cdr3b", "Unique CDR3B", "CDR3B"),
            ("n_mhc_alleles", "Unique MHC Alleles", "MHC Alleles")
        ]
        plt.figure(figsize=(15, 8))
        for i, (col, label, title) in enumerate(stats_to_plot, 1):
            plt.subplot(2, 3, i)
            sns.barplot(x="dataset_name", y=col, data=df, hue="dataset_name", legend=False)
            plt.title(title)
            plt.xlabel("Dataset")
            plt.ylabel(label)
            plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/summary_statistics.png')
    
    def plot_all_distributions(self, save_dir='.', figsize=(10, 3), map_colors=None,
                            width_ratios=(0.4, 1), wspace=0.07,
                            legend_title="Datasets", legend_ncol=1,
                            legend_bbox=(0.9, 0.5)):
        """
        Shared y-axis, controllable subplot widths, single legend on right.
        - figsize: overall figure size
        - width_ratios: relative widths of the two plots, e.g., (1.5, 1)
        - wspace: horizontal space between subplots
        - legend_bbox: (x, y) anchor for legend on the right side
        """

        map_colors = map_colors or 'tab10'
        distribution_keys = [
            ("peptide_length_prop", "Peptide Length", "Proportion"),
            ("mhc_distribution_prop", "MHC Allele", "Proportion")
        ]

        df = self.to_df()
        fig, axes = plt.subplots(
            1, 2, figsize=figsize, sharey=True,
            gridspec_kw={'width_ratios': width_ratios, 'wspace': wspace}
        )

        all_handles, all_labels = [], []

        for i, (key, xlabel, ylabel) in enumerate(distribution_keys):
            ax = axes[i]
            dist_list = []
            for _, row in df.iterrows():
                dist = row.get(key, {}) if isinstance(row.get(key, {}), dict) else {}
                if not dist:
                    continue
                for k, v in dist.items():
                    dist_list.append({"dataset_name": row["dataset_name"], xlabel: str(k), ylabel: v})

            dist_df = pd.DataFrame(dist_list)
            if dist_df.empty:
                ax.set_visible(False)
                continue

            # Category ordering
            if key == "mhc_distribution_prop":
                dist_df[xlabel] = dist_df[xlabel].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)
                dist_df[xlabel] = dist_df[xlabel].str.replace('HLA-A*011', 'HLA-A*11', regex=False)
                categories = sorted([c for c in dist_df[xlabel].unique() if c != "Others"])
                if "Others" in dist_df[xlabel].unique():
                    categories.append("Others")
            else:  # peptide_length_prop / score_distribution
                categories = sorted(dist_df[xlabel].unique(), key=lambda z: int(z))

            sns.barplot(
                data=dist_df, x=xlabel, y=ylabel, hue="dataset_name",
                palette=map_colors, width=0.8, errorbar=None, order=categories, ax=ax
            )

            # Style
            ax.set_xlabel(xlabel, fontsize=12)
            ax.set_ylabel(ylabel if i == 0 else "", fontsize=12)
            ax.set_ylim(0, 1)
            ax.tick_params(axis="x", rotation=45, labelsize=10)
            ax.tick_params(axis="y", labelsize=10)
            sns.despine(ax=ax)
            ax.grid(axis="y", linestyle="--", linewidth=0.5, color="gray", alpha=0.7)

            # collect legend entries once
            h, l = ax.get_legend_handles_labels()
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            all_handles.extend(h)
            all_labels.extend(l)

        # Deduplicate while preserving order
        uniq = {}
        for lbl, h in zip(all_labels, all_handles):
            if lbl not in uniq:
                uniq[lbl] = h

        # Legend on the RIGHT
        fig.legend(
            handles=list(uniq.values()),
            labels=list(uniq.keys()),
            title=legend_title,
            frameon=False,
            loc='center left',
            bbox_to_anchor=legend_bbox,   # (x>1 puts it outside on the right)
            ncol=legend_ncol
        )

        # Save ensuring legend isn’t clipped
        fig.savefig(f"{save_dir}/plot_distribution.png", dpi=300, bbox_inches='tight')
        plt.close(fig)


def plot_all_epitopes_counts(tcr_peptide_counts, output_plot_path, output_data_path, top_n_zoom=10):
    sns.set_theme(style="whitegrid")

    # Sort and rank
    df = tcr_peptide_counts.sort_values("unique_tcr_count", ascending=False).reset_index(drop=True)
    print('Count of peptides with only 1 unique TCR:')
    print(len(tcr_peptide_counts[tcr_peptide_counts['unique_tcr_count'] == 1]))

    df["rank"] = df.index + 1

    # Create figure with better layout handling
    fig, ax = plt.subplots(figsize=(15, 7), constrained_layout=True)

    # Main plot using hue to silence palette warning
    df["top10"] = df["rank"] <= 10
    sns.barplot(
        x="rank", y="unique_tcr_count", hue="top10", data=df,
        width=1, edgecolor="none",
        ax=ax, palette={True: "orange", False: "skyblue"}, legend=False
    )

    ax.set_yscale("log")
    ax.set_xlabel("Peptides", fontsize=20)
    ax.set_ylabel("TCR counts (log scale)", fontsize=20)

    # X-ticks with limited labels
    xticks = [1, 300, 600, 900, 1112]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=0)
    ax.tick_params(axis='y', labelsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray')
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        #spine.set_linewidth(2)
    

    axins = inset_axes(ax, width="40%", height="55%", loc='upper right', borderpad=3.5)
    bars = sns.barplot(
                x="epitope", y="unique_tcr_count", #hue="epitope",
                data=df.head(top_n_zoom),
                ax=axins, color='orange', legend=False
            )

    # Annotate top N
    for i, bar in enumerate(bars.patches[:top_n_zoom]):
        height = bar.get_height()
        axins.annotate(f"{int(height)}", (bar.get_x() + bar.get_width()/2., height),
                        ha='center', va='bottom', fontsize=10, color='black')

    # Inset with correct tick fix
    print(df["epitope"].head(top_n_zoom).tolist())
    axins.set_yscale("log")
    axins.set_xticks(range(top_n_zoom))
    axins.set_xticklabels(df["epitope"].head(top_n_zoom),
                          rotation=45, ha='right', fontsize=12)
    axins.set_title(f"Top {top_n_zoom} Epitopes", fontsize=15)
    axins.set_ylabel("")
    axins.set_xlabel("")
    for spine in axins.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(2)

    # Save output
    df.to_csv(f'{output_data_path}/tcr_peptide_counts.csv', index=False)
    os.makedirs(output_plot_path, exist_ok=True)
    output_file = os.path.join(output_plot_path, "unique_tcr_counts_logscale_fancy_peptide_rank.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
        
    