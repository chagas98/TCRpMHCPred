import numpy as np
from typing import List, Literal, Optional
from tqdm import tqdm
import umap.umap_ as umap
import matplotlib.pyplot as plt
import seaborn as sns
import Levenshtein
import pandas as pd
import json
import argparse
import esm
import torch
from itertools import product
from sklearn.preprocessing import normalize


class PeptideClustering:
    def __init__(self, peptides: List[str], 
                 set_name: Optional[List[str]] = None, 
                 min_dist: float = 0.05,
                 n_neighbors: int = 50,
                 metric: Literal['levenshtein', 'esm'] = 'levenshtein',
                 norm: bool = False):
        """
        Initialize the PeptideClustering class with a list of peptides.
        :param peptides: List of peptide sequences.
        :param set_name: Optional list of labels for each peptide.
        """
        self.device = ("cuda" if torch.cuda.is_available() else "cpu")
        self.peptides = peptides
        self.set_name = set_name
        self.n_neighbors = n_neighbors
        self.min_dist = min_dist
        self.metric = metric.lower()
        self.norm = norm
        self.model = None
        self.dist_matrix = None
        self.coords = None

        self.df_peptides = pd.DataFrame({
            "peptide": self.peptides,
            "Set_name": self.set_name if set_name is not None else [""] * len(peptides)
        })

        if self.metric not in ['levenshtein', 'esm']:
            raise ValueError("Unsupported metric: {}. Use 'levenshtein' or 'esm'.".format(self.metric))
        if self.metric == 'esm':
            print("Using ESM embeddings for peptide clustering.")
            self.embed_peptides_esm(norm=self.norm)
        else:
            print("Using Levenshtein distance for peptide clustering.")
            self.compute_distance_matrix(norm=self.norm)
        
        print("Computing UMAP coordinates...")
        coords = self.get_coords(n_neighbors=self.n_neighbors, min_dist=self.min_dist)

    def save_coords(self):
        df = pd.DataFrame({
            "UMAP1": self.coords[:, 0],
            "UMAP2": self.coords[:, 1]
            })
        df.reset_index(drop=True, inplace=True)
        return df
     
    def compute_distance_matrix(self, norm) -> np.ndarray:
        """
        Compute a distance matrix for a list of peptides using the specified method.
        :return: A square distance matrix.
        """
        n = len(self.peptides)
        dist_matrix = np.zeros((n, n))

        for i in tqdm(range(n), desc="Levenshtein distances"):
            for j in range(i + 1, n):
                dist = Levenshtein.distance(self.peptides[i], self.peptides[j])
                
                if norm:
                    # Normalize the distance by the length of the longer peptide
                    max_len = max(len(self.peptides[i]), len(self.peptides[j]))
                    dist = dist / max_len if max_len > 0 else 0
                
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

        self.dist_matrix = dist_matrix
        return self.dist_matrix
    
    def embed_peptides_esm(self, model_name="esm2_t33_650M_UR50D", norm=False) -> np.ndarray:
        if self.model is None or self.batch_converter is None:
            self.model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
            self.model.eval().to(self.device)
            self.batch_converter = alphabet.get_batch_converter()

        data = [("peptide_{}".format(i), seq) for i, seq in enumerate(self.peptides)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)

        with torch.no_grad():
            results = self.model(batch_tokens, repr_layers=[6])
            token_representations = results["representations"][6]

        embeddings = []
        for i, (_, seq) in enumerate(data):
            emb = token_representations[i, 1:len(seq)+1].mean(0)
            embeddings.append(emb.cpu().numpy())

        if norm:
            embeddings = normalize(embeddings, norm='l2')
        
        self.embeddings_cache = np.vstack(embeddings)
        return self.embeddings_cache
    
    def get_coords(self, n_neighbors: int = 5, min_dist: float = 0.0) -> np.ndarray:
        """
        Compute UMAP coordinates for the peptides based on the distance matrix.
        :return: UMAP coordinates of the peptides.
        """
        if not hasattr(self, 'dist_matrix'):
            print("Computing distance matrix...")
            self.compute_distance_matrix()

        match self.metric:
            case 'esm':
                reducer = umap.UMAP(metric='cosine',
                                    n_neighbors=n_neighbors,        # adjust based on data density
                                    min_dist=min_dist,          # controls compactness of clusters
                                    n_components=2,        # for 2D visualization
                                    random_state=42)
                self.coords = reducer.fit_transform(self.embeddings_cache)
            
            case'levenshtein':
                reducer = umap.UMAP(metric='precomputed',
                                    n_neighbors=n_neighbors,        # adjust based on data density
                                    min_dist=min_dist,          # controls compactness of clusters
                                    n_components=2,        # for 2D visualization
                                    random_state=42)
                self.coords = reducer.fit_transform(self.dist_matrix)
            case _:
                raise ValueError("Unsupported metric: {}".format(self.metric))
        
    def plot(self, test_set: str, labels: Optional[List[str]] = None, 
             label_name: str = "Set_name", save_dir: Optional[str] = None, 
             map_colors: dict = None, peptides_to_annotate: List[str] = None,
             show_legend = True, legend_position: str = 'upper right', 
             figsize: tuple = (3, 3), outfile: str = '', mv_legend: str = 'full') -> None:
        """
        Plot the UMAP coordinates of the peptides, colored by set_name if provided.
        Allowed values for legend_position are 'upper right' and 'lower center'.
        :param output_file: Path to save the plot.
        """
        print("Plotting UMAP coordinates...")

        if legend_position not in ['upper right', 'lower center', 'upper center']:
            raise ValueError("Allowed values for legend_position are 'upper right', 'lower center', and 'upper center'.")

        fig, ax = plt.subplots(figsize=figsize)

        df = pd.DataFrame({
            "UMAP1": self.coords[:, 0],
            "UMAP2": self.coords[:, 1],
            "Set_name": self.set_name,
            "Labels": labels
        })

        if labels is not None:
            markers = {'Train': 'o', 'Test': 'X'}
            df['Split'] = df['Set_name'].apply(lambda x: 'Test' if x == test_set else "Train")
            if map_colors is None:  
                map_colors = 'tab10'
            
            # Plot the scatter without legend
            sns.scatterplot(
                data=df[df["Set_name"] != test_set], x="UMAP1", y="UMAP2",
                hue='Labels', palette=map_colors, hue_order=list(map_colors.keys()),
                s=20, alpha=0.7, legend=show_legend, ax=ax
            )

            handles, labels_ = ax.get_legend_handles_labels()
            
            if show_legend == 'brief':
                ax.get_legend().remove()   # remove legend from main plot
                legend_fig, legend_ax = plt.subplots(figsize=(2, 3))
                legend_ax.axis("off")
                legend_ax.legend(handles, labels_, loc="center", frameon=False, title=label_name)
                legend_fig.savefig(f"{save_dir}/legend.png", dpi=300)
                plt.close(legend_fig)

            #sns.scatterplot(data=df[df["Set_name"] == test_set], x="UMAP1", y="UMAP2", color="black", 
            #                s=40, legend=show_legend, marker='X', ax=ax)
            
            if mv_legend:
                # Position legend based on allowed values
                if legend_position == 'upper right':
                    sns.move_legend(ax, "upper right", bbox_to_anchor=(1.5, 1), title=None, frameon=False)
                    plt.subplots_adjust(left=0.5)
                elif legend_position == 'lower center':
                    sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, -0.4), ncol=len(labels), title=None, frameon=False)
                    plt.subplots_adjust(bottom=0.2)
                elif legend_position == 'upper center':
                    sns.move_legend(ax, "upper center", bbox_to_anchor=(.5, 1.15), ncol=len(labels), title=None, frameon=False)
                    plt.subplots_adjust(top=0.8)

            df['peptide'] = self.peptides
            df.to_csv('data/umap_coords.csv', index=False)
            pep_test = self.df_peptides[self.df_peptides["Set_name"] == test_set].reset_index(drop=True)["peptide"]
            coords_test = df[df["Set_name"] == test_set].reset_index(drop=True)
            for i, row in coords_test.iterrows():
                if peptides_to_annotate is not None:
                    peptide = pep_test.iloc[i]
                    if peptide in peptides_to_annotate:
                        plt.annotate(peptide, (row["UMAP1"], row["UMAP2"]), fontsize=7, color='black', alpha=0.8)
                
            plt.xlabel("UMAP 1")
            plt.ylabel("UMAP 2")
            plt.tight_layout()
            plt.savefig("{}/umap_{}_{}_norm{}_{}_{}_{}.png".format(
                save_dir, outfile, self.metric, self.norm, 
                label_name, self.n_neighbors, self.min_dist), dpi=300)
            plt.close()
