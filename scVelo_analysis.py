import seaborn as sns
from scvelo.plotting.utils import (
    interpret_colorkey,
    is_categorical,
    savefig_or_show,
    set_colors_for_categorical_obs,
    strings_to_categoricals,
    to_list,
)
from scipy.sparse import issparse
import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import unitvelo as utv
import pymannkendall as mk

###############################################
#                                             #
# Whole-testis cell RNA velocity generalized  #
#      through unitVelo modeling              #
#                                             #
###############################################

###############
#   Wildtype  #
###############

WT_loom = anndata.read_loom("./WT.loom")
WT_loom.var_names_make_unique()
sample_obs = pd.read_csv("./cell_ID_obs.csv")
umap = pd.read_csv("./cell_embeddings.csv")
cell_clusters = pd.read_csv("./cell_clusters_obs.csv")
cell_types = pd.read_csv("./cell_type_obs.csv")

adata_WT = WT_loom[np.isin(WT_loom.obs_names, sample_obs["x"]), :]

adata_WT_index = pd.DataFrame(adata_WT.obs.index)
adata_WT_index = adata_WT_index.rename(columns={0: 'Cell ID'})
adata_WT_index = adata_WT_index.rename(columns={"CellID": 'Cell ID'})


def rep(x): return x.split("-")[0]


adata_WT_index["Cell ID"] = adata_WT_index["Cell ID"].apply(rep)

# Change the column name of umap to unify the same column name Cell ID
umap = umap.rename(columns={'Unnamed: 0': 'Cell ID'})
# Filter the cell ID of adata_WT_index in umap
umap = umap[np.isin(umap["Cell ID"], adata_WT_index["Cell ID"])]
umap = umap.drop_duplicates(subset=["Cell ID"])  # Remove duplicate values
# Merge with umap data based on adata_WT_index Cell ID order
umap_ordered = adata_WT_index.merge(umap, on="Cell ID")
umap_ordered = umap_ordered.iloc[:, 1:]
# Incorporate UMAP into the adata_WT object
adata_WT.obsm['X_umap'] = umap_ordered.values


cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_clusters = cell_clusters[np.isin(
    cell_clusters["Cell ID"], adata_WT_index["Cell ID"])]
cell_clusters = cell_clusters.drop_duplicates(
    subset=["Cell ID"])  # Remove duplicate values
cell_clusters_ordered = adata_WT_index.merge(cell_clusters, on="Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:, 1:]
adata_WT.uns['clusters'] = cell_clusters_ordered.values

cell_types = cell_types.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_types = cell_types[np.isin(
    cell_types["Cell ID"], adata_WT_index["Cell ID"])]
cell_types = cell_types.drop_duplicates(
    subset=["Cell ID"])  # Remove duplicate values
cell_types_ordered = adata_WT_index.merge(cell_types, on="Cell ID")
cell_types_ordered = cell_types_ordered.iloc[:, 1:]
adata_WT.obs['celltype'] = cell_types_ordered.values.flatten()

velo = utv.config.Configuration()
velo.R2_ADJUST = True
velo.IROOT = None
velo.FIT_OPTION = '1'
velo.GPU = 0
velo.N_TOP_GENES = 2500
velo.MAX_ITER = 14000

adata_WT = utv.run_model(adata_WT, "celltype", config_file=velo)

# Visualization

scv.tl.paga(adata_WT, groups='celltype')
ident_colours = ["#e69dc5", "#7BAFDE", "#a6d608", "#F2A874", "#456d88",
                 "#52a95f", "#9db7a5", "#9b77a2", "#fbb1a2", "#bbbbff", "#56c4f3"]
scv.pl.scatter(adata_WT,
               save="./figures/1e.Terminal_states_WT.pdf",
               color=["end_points"],
               xlim=[-12, 10],
               ylim=[-12, 14],
               size=20,
               figsize=((7, 7)),
               legend_fontsize=9,
               title='',
               show=False)
scv.pl.velocity_embedding_stream(adata_WT,
                                 save="./figures/WT.Stream.svg",
                                 basis='X_umap',
                                 color="celltype",
                                 palette=ident_colours,
                                 xlim=[-12, 10],
                                 ylim=[-12, 14],
                                 fontsize=14,
                                 add_text_pos=[0.48, 0.52],
                                 size=80,
                                 alpha=0.8,
                                 figsize=((4, 4)),
                                 legend_fontsize=9,
                                 title='',
                                 arrow_size=0.8,
                                 linewidth=0.6,
                                 show=False)
scv.pl.paga(adata_WT,
            basis='umap',
            save="./figures/WT.Paga_WT.pdf",
            size=50,
            alpha=.2,
            min_edge_width=2,
            node_size_scale=1.5,
            palette=ident_colours,
            xlim=[-12, 10],
            ylim=[-12, 14],
            figsize=((7, 7)),
            legend_fontsize=9,
            show=False)
scv.pl.scatter(adata_WT,
               save="./figures/WT.latent_time_WT.pdf",
               color="latent_time",
               color_map="gnuplot",
               xlim=[-12, 10],
               ylim=[-12, 14],
               size=20,
               figsize=((7, 7)),
               legend_fontsize=9,
               title='',
               show=False)

###############
#   Knockout  #
###############

KO_loom = anndata.read_loom("./Hsf5.KO.loom")
KO_loom.var_names_make_unique()
sample_obs = pd.read_csv("./cell_ID_obs.csv")
umap = pd.read_csv("./cell_embeddings.csv")
cell_clusters = pd.read_csv("./cell_clusters_obs.csv")
cell_types = pd.read_csv("./cell_type_obs.csv")

adata_KO = KO_loom[np.isin(KO_loom.obs_names, sample_obs["x"]), :]

adata_KO_index = pd.DataFrame(adata_KO.obs.index)
adata_KO_index = adata_KO_index.rename(columns={0: 'Cell ID'})
adata_KO_index = adata_KO_index.rename(columns={"CellID": 'Cell ID'})


def rep(x): return x.split("-")[0]


adata_KO_index["Cell ID"] = adata_KO_index["Cell ID"].apply(rep)

umap = umap.rename(columns={'Unnamed: 0': 'Cell ID'})
umap = umap[np.isin(umap["Cell ID"], adata_KO_index["Cell ID"])]
umap = umap.drop_duplicates(subset=["Cell ID"])
umap_ordered = adata_KO_index.merge(umap, on="Cell ID")
umap_ordered = umap_ordered.iloc[:, 1:]
adata_KO.obsm['X_umap'] = umap_ordered.values


cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_clusters = cell_clusters[np.isin(
    cell_clusters["Cell ID"], adata_KO_index["Cell ID"])]
cell_clusters = cell_clusters.drop_duplicates(subset=["Cell ID"])
cell_clusters_ordered = adata_KO_index.merge(cell_clusters, on="Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:, 1:]
adata_KO.uns['clusters'] = cell_clusters_ordered.values

cell_types = cell_types.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_types = cell_types[np.isin(
    cell_types["Cell ID"], adata_KO_index["Cell ID"])]
cell_types = cell_types.drop_duplicates(subset=["Cell ID"])
cell_types_ordered = adata_KO_index.merge(cell_types, on="Cell ID")
cell_types_ordered = cell_types_ordered.iloc[:, 1:]
adata_KO.obs['celltype'] = cell_types_ordered.values.flatten()

adata_KO.var['gene'] = adata_KO.var.index

velo = utv.config.Configuration()
velo.R2_ADJUST = True
velo.IROOT = "Spg"
velo.FIT_OPTION = '1'
velo.GPU = 0
velo.N_TOP_GENES = 2000
velo.MAX_ITER = 14000

df = pd.DataFrame(index=adata_KO.obs_names).reset_index()
adata_KO.uns['root_key'] = df.index[df['CellID'].isin(
    adata_KO.obs_names[adata_KO.obs["celltype"] == "Spg"])][0]
adata_KO = utv.run_model(adata_KO, "celltype", config_file=velo)

# Visualization

scv.tl.paga(adata_KO, groups='celltype')
ident_colours = ["L", "P-LIKE", "SPG", "Z", "EP", "MP", "PL"]
ident_colours = ["#7BAFDE", "#ff6e5f", "#52a95f",
                 "#9db7a5", "#9b77a2", "#bbbbff", "#56c4f3"]
scv.pl.scatter(adata_KO,
               save="./figures/KO.Terminal_states_KO.pdf",
               color=["end_points"],
               xlim=[-12, 10],
               ylim=[-12, 14],
               size=20,
               figsize=((7, 7)),
               legend_fontsize=9,
               title='',
               show=False)
scv.pl.velocity_embedding_stream(adata_KO,
                                 save="./figures/KO.Stream.svg",
                                 basis='X_umap',
                                 color="celltype",
                                 palette=ident_colours,
                                 xlim=[-12, 10],
                                 ylim=[-12, 14],
                                 fontsize=14,
                                 add_text_pos=[0.48, 0.52],
                                 size=80,
                                 alpha=0.8,
                                 figsize=((4, 4)),
                                 legend_fontsize=9,
                                 title='',
                                 arrow_size=0.8,
                                 linewidth=0.6,
                                 show=False)
scv.pl.paga(adata_KO,
            basis='umap',
            save="./figures/KO.Paga_KO.pdf",
            size=50,
            alpha=.2,
            min_edge_width=2,
            node_size_scale=1.5,
            edge_width_scale=1.5,
            palette=ident_colours,
            xlim=[-12, 10],
            ylim=[-12, 14],
            figsize=((7, 7)),
            legend_fontsize=9,
            show=False)
scv.pl.scatter(adata_KO,
               save="./figures/KO.latent_time_KO.pdf",
               color="latent_time",
               color_map="gnuplot",
               xlim=[-12, 10],
               ylim=[-12, 14],
               size=20,
               figsize=((7, 7)),
               legend_fontsize=9,
               title='',
               show=False)


###############################################
#                                             #
# Pachytene stage cell RNA velocity generalized#
#      through dynamical modeling             #
#                                             #
###############################################

WT_loom = anndata.read_loom("./data/WT.loom")
WT_loom.var_names_make_unique()
sample_obs = pd.read_csv("./data/pachytene_cell_ID_obs_remove_KO_mP.csv")
umap = pd.read_csv("./data/pachytene_cell_embeddings.csv")
cell_clusters = pd.read_csv("./data/pachytene_cell_clusters_obs.csv")
cell_types = pd.read_csv("./data/pachytene_cell_type_obs.csv")

WT_loom_sub = WT_loom[np.isin(WT_loom.obs_names, sample_obs["x"]), :]

adata_WT = WT_loom_sub

adata_WT_index = pd.DataFrame(adata_WT.obs.index)
adata_WT_index = adata_WT_index.rename(columns={0: 'Cell ID'})
adata_WT_index = adata_WT_index.rename(columns={"CellID": 'Cell ID'})


def rep(x): return x.split("-")[0]


adata_WT_index["Cell ID"] = adata_WT_index["Cell ID"].apply(rep)

umap = umap.rename(columns={'Unnamed: 0': 'Cell ID'})
umap = umap[np.isin(umap["Cell ID"], adata_WT_index["Cell ID"])]
umap = umap.drop_duplicates(subset=["Cell ID"])
umap_ordered = adata_WT_index.merge(umap, on="Cell ID")
umap_ordered = umap_ordered.iloc[:, 1:]
adata_WT.obsm['X_umap'] = umap_ordered.values


cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_clusters = cell_clusters[np.isin(
    cell_clusters["Cell ID"], adata_WT_index["Cell ID"])]
cell_clusters = cell_clusters.drop_duplicates(subset=["Cell ID"])
cell_clusters_ordered = adata_WT_index.merge(cell_clusters, on="Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:, 1:]
adata_WT.uns['clusters'] = cell_clusters_ordered.values

cell_types = cell_types.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_types = cell_types[np.isin(
    cell_types["Cell ID"], adata_WT_index["Cell ID"])]
cell_types = cell_types.drop_duplicates(subset=["Cell ID"])
cell_types_ordered = adata_WT_index.merge(cell_types, on="Cell ID")
cell_types_ordered = cell_types_ordered.iloc[:, 1:]
adata_WT.obs['celltype'] = cell_types_ordered.values.flatten()


scv.pp.filter_and_normalize(adata_WT, min_shared_counts=30)
scv.pp.moments(adata_WT, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata_WT, n_jobs=8)

scv.tl.velocity(adata_WT, mode="dynamical",
                root_key='root_key', end_key='end_key')
scv.tl.velocity_graph(adata_WT, n_jobs=8)
scv.tl.terminal_states(adata_WT)
scv.tl.latent_time(adata_WT, root_key='root_key')

KO_loom = anndata.read_loom("./data/Hsf5.KO.loom")
KO_loom.var_names_make_unique()
sample_obs = pd.read_csv("./data/pachytene_cell_ID_obs_remove_KO_mP.csv")
umap = pd.read_csv("./data/pachytene_cell_embeddings.csv")
cell_clusters = pd.read_csv("./data/pachytene_cell_clusters_obs.csv")
cell_types = pd.read_csv("./data/pachytene_cell_type_obs.csv")
retain_genes = list(adata_WT.var['fit_likelihood'][np.where(
    ~np.isnan(adata_WT.var['fit_likelihood']))[0]].index)

KO_loom_sub = KO_loom[np.isin(KO_loom.obs_names, sample_obs["x"]), :]

adata_KO = KO_loom_sub

adata_KO_index = pd.DataFrame(adata_KO.obs.index)
adata_KO_index = adata_KO_index.rename(columns={0: 'Cell ID'})
adata_KO_index = adata_KO_index.rename(columns={"CellID": 'Cell ID'})
adata_KO_index["Cell ID"] = adata_KO_index["Cell ID"].apply(rep)

umap = umap.rename(columns={'Unnamed: 0': 'Cell ID'})
umap = umap[np.isin(umap["Cell ID"], adata_KO_index["Cell ID"])]
umap = umap.drop_duplicates(subset=["Cell ID"])
umap_ordered = adata_KO_index.merge(umap, on="Cell ID")
umap_ordered = umap_ordered.iloc[:, 1:]
adata_KO.obsm['X_umap'] = umap_ordered.values


cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_clusters = cell_clusters[np.isin(
    cell_clusters["Cell ID"], adata_KO_index["Cell ID"])]
cell_clusters = cell_clusters.drop_duplicates(subset=["Cell ID"])
cell_clusters_ordered = adata_KO_index.merge(cell_clusters, on="Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:, 1:]
adata_KO.uns['clusters'] = cell_clusters_ordered.values

cell_types = cell_types.rename(columns={'Unnamed: 0': 'Cell ID'})
cell_types = cell_types[np.isin(
    cell_types["Cell ID"], adata_KO_index["Cell ID"])]
cell_types = cell_types.drop_duplicates(subset=["Cell ID"])
cell_types_ordered = adata_KO_index.merge(cell_types, on="Cell ID")
cell_types_ordered = cell_types_ordered.iloc[:, 1:]
adata_KO.obs['celltype'] = cell_types_ordered.values.flatten()

df = pd.DataFrame(index=adata_KO.obs_names).reset_index()
adata_KO.uns['root_key'] = df.index[df['CellID'].isin(
    adata_KO.obs_names[adata_KO.obs["celltype"] == "eP"])][0]
adata_KO.uns['end_key'] = df.index[df['CellID'].isin(
    adata_KO.obs_names[adata_KO.obs["celltype"] == "P-like"])][0]

adata_KO.var['gene'] = adata_KO.var.index

adata_KO_uni = adata_KO.copy()

scv.pp.filter_and_normalize(
    adata_KO, min_shared_counts=30, retain_genes=retain_genes)
scv.pp.moments(adata_KO, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata_KO, n_jobs=8)

velo = utv.config.Configuration()
velo.R2_ADJUST = True
velo.IROOT = "4"
velo.FIT_OPTION = '1'
velo.GPU = 0
velo.N_TOP_GENES = 3000
velo.MAX_ITER = 1500
adata_KO_uni = utv.run_model(adata_KO_uni, "celltype", config_file=velo)

scv.tl.velocity(adata_KO, mode="dynamical",
                root_key='root_key', end_key='end_key')
scv.tl.velocity_graph(adata_KO, n_jobs=8)
scv.tl.terminal_states(adata_KO)
scv.tl.latent_time(adata_KO, root_key='root_key')
adata_KO.obs.latent_time = adata_KO_uni.obs.latent_time

ident_colours_WT = ["#9b77a2", "#fbb1a2", "#bbbbff"]
scv.pl.velocity_embedding_stream(
    adata_WT, basis='X_umap', color="celltype", palette=ident_colours_WT, size=80, alpha=0.8,)
scv.pl.scatter(adata_WT, color=["end_points"], show=True)
scv.pl.scatter(adata_WT, color='latent_time', color_map='gnuplot', size=80)

ident_colours_KO = ["#ff6e5f", "#9b77a2", "#bbbbff"]
scv.pl.velocity_embedding_stream(
    adata_KO, basis='X_umap', color="celltype", palette=ident_colours_KO, size=80, alpha=0.8,)
scv.pl.scatter(adata_KO, color=["end_points"], show=True)
scv.pl.scatter(adata_KO, color='latent_time', color_map='gnuplot', size=80)


adata_KO.write('data/adata_KO.h5ad', compression='gzip')
adata_WT.write('data/adata_WT.h5ad', compression='gzip')

adata_KO = scv.read('data/adata_KO.h5ad')
adata_WT = scv.read('data/adata_WT.h5ad')


top_genes = adata_WT.var['fit_likelihood'][np.where(
    ~np.isnan(adata_WT.var['fit_likelihood']))[0]].index
var_names = top_genes
adata = adata_WT
var_names = [
    name for name in var_names if name in adata_WT.var_names and name in adata_KO.var_names]
xkey = "Mu"
tkey = "latent_time"
time = adata.obs[tkey].values
time = time[np.isfinite(time)]
X = (
    adata[:, var_names].layers[xkey]
    if xkey in adata.layers.keys()
    else adata[:, var_names].X
)
df = pd.DataFrame(X[np.argsort(time)], columns=var_names)
n_convolve = 30
weights = np.ones(n_convolve) / n_convolve
for gene in var_names:
    # TODO: Handle exception properly
    try:
        df[gene] = np.convolve(df[gene].values, weights, mode="same")
    except ValueError as e:
        pass  # e.g. all-zero counts or nans cannot be convolved
max_sort = np.argsort(np.argmax(df.values, axis=0))
df = pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort])
sort_top_genes = df.columns.tolist()

palette = "viridis"
col_color = 'celltype'
col_colors = to_list(col_color)
col_color = []
for _, col in enumerate(col_colors):
    if not is_categorical(adata, col):
        obs_col = adata.obs[col]
        cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
        adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
        col += "_categorical"
        set_colors_for_categorical_obs(adata, col, palette)
    col_color.append(interpret_colorkey(adata, col)[np.argsort(time)])
col_color_WT = col_color

adata = adata_KO
time = adata.obs[tkey].values
time = time[np.isfinite(time)]
X = (
    adata[:, var_names].layers[xkey]
    if xkey in adata.layers.keys()
    else adata[:, var_names].X
)
df_KO = pd.DataFrame(X[np.argsort(time)], columns=var_names)
weights = np.ones(n_convolve) / n_convolve
for gene in var_names:
    # TODO: Handle exception properly
    try:
        df_KO[gene] = np.convolve(df_KO[gene].values, weights, mode="same")
    except ValueError as e:
        pass  # e.g. all-zero counts or nans cannot be convolved
max_sort = np.argsort(np.argmax(df_KO.values, axis=0))
df_KO = pd.DataFrame(df_KO.values[:, max_sort],
                     columns=df_KO.columns[max_sort])

palette = "viridis"
col_color = 'celltype'
col_colors = to_list(col_color)
col_color = []
for _, col in enumerate(col_colors):
    if not is_categorical(adata, col):
        obs_col = adata.obs[col]
        cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
        adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
        col += "_categorical"
        set_colors_for_categorical_obs(adata, col, palette)
    col_color.append(interpret_colorkey(adata, col)[np.argsort(time)])

col_color_KO = col_color

cm = sns.clustermap(pd.merge(df.T, df_KO.T.reindex(index=df.T.index), left_index=True, right_index=True),
                    standard_scale=0,
                    cmap="viridis",
                    col_cluster=False,
                    row_cluster=False,
                    yticklabels=False,
                    xticklabels=False,
                    dendrogram_ratio=(0, 0),
                    figsize=(6.8, 8),
                    col_colors=np.concatenate(
                        (col_color_WT, col_color_KO), axis=1),
                    cbar_pos=None
                    )

plt.savefig("Driver_gene_all.png", dpi=600)


#####################################
#                                   #
#          Binding Gene             #
#                                   #
#####################################


Binding_list = pd.read_csv("./data/HSF5_binding_gene.csv")["Gene_name"]

WT_Mu_time_df = pd.DataFrame(adata_WT.layers['Mu'], columns=adata_WT.var_names)
WT_Mu_time_df['latent_time'] = list(adata_WT.obs["latent_time"])
WT_Mu_time_df['cell_type'] = list(adata_WT.obs["celltype"])
WT_Mu_time_df = WT_Mu_time_df.sort_values(by='latent_time')

KO_Mu_time_df = pd.DataFrame(adata_KO.layers['Mu'], columns=adata_KO.var_names)
KO_Mu_time_df['latent_time'] = list(adata_KO.obs["latent_time"])
KO_Mu_time_df['cell_type'] = list(adata_KO.obs["celltype"])
KO_Mu_time_df = KO_Mu_time_df.sort_values(by='latent_time')

var_names = Binding_list
Bind_Mann_Kendall_result = pd.DataFrame(columns=["trend_WT", "h_WT", "p_WT", "z_WT", "Tau_WT", "s_WT", "var_s_WT", "slope_WT",
                                        "intercept_WT", "trend_KO", "h_WT", "p_KO", "z_WT", "Tau_WT", "s_WT", "var_s_WT", "slope_KO", "intercept_KO"])
repression_abnormal_gene = []
induction_abnormal_gene = []
ns_gene = []

for gene_name in var_names:
    gene_expression_data = WT_Mu_time_df[gene_name]

    # Mann-Kendall test
    trend_WT, h_WT, p_WT, z_WT, Tau_WT, s_WT, var_s_WT, slope_WT, intercept_WT = mk.original_test(
        WT_Mu_time_df[gene_name])
    trend_KO, h_WT, p_KO, z_WT, Tau_WT, s_WT, var_s_WT, slope_KO, intercept_KO = mk.original_test(
        KO_Mu_time_df[gene_name])

    Bind_Mann_Kendall_result.loc[gene_name] = [trend_WT, h_WT, p_WT, z_WT, Tau_WT, s_WT, var_s_WT,
                                               slope_WT, intercept_WT, trend_KO, h_WT, p_KO, z_WT, Tau_WT, s_WT, var_s_WT, slope_KO, intercept_KO]

Bind_Mann_Kendall_result.to_csv(
    "./Driver_binding_gene_Mann_Kendall_test.csv", header=True)

#####################################
#                                   #
#          Visualization            #
#                                   #
#####################################

Abnormal_repression_gene = ["Msh4", "Hat1", "Meiob", "Sycp1", "Tesmin"]
scv.pl.scatter(adata_WT,  x='latent_time', y=Abnormal_repression_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./Abnormal_repression_gene_scatter_in_WT.pdf", xlabel="Inferred Cell Time", ylabel="")
scv.pl.scatter(adata_KO,  x='latent_time', y=Abnormal_repression_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./Abnormal_repression_gene_scatter_in_KO.pdf", xlabel="Inferred Cell Time", ylabel="")

Abnormal_induction_gene = ["Scaper", "Piwil2", "Hpse2", "Setx", "Spag9"]
scv.pl.scatter(adata_WT,  x='latent_time', y=Abnormal_induction_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./Abnormal_induction_gene_scatter_in_WT.pdf", xlabel="Inferred Cell Time", ylabel="")
scv.pl.scatter(adata_KO,  x='latent_time', y=Abnormal_induction_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./Abnormal_induction_gene_scatter_in_KO.pdf", xlabel="Inferred Cell Time", ylabel="")

Induction_driver_gene = ["Larp7", "Cep57", "Immt", "Cgrrf1", "Tprn"]
scv.pl.scatter(adata_WT,  x='latent_time', y=Induction_driver_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./WT_Induction_driver_gene_time_unsplied.pdf", xlabel="Inferred Cell Time", ylabel="")
scv.pl.scatter(adata_KO,  x='latent_time', y=Induction_driver_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./KO_Induction_driver_gene_time_unsplied.pdf", xlabel="Inferred Cell Time", ylabel="")

Repression_driver_gene = ["Msh4", "Hat1", "Meiob", "Sycp1", "Tesmin"]
scv.pl.scatter(adata_WT,  x='latent_time', y=Repression_driver_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./WT_Repression_driver_gene_time_unsplied.pdf", xlabel="Inferred Cell Time", ylabel="")
scv.pl.scatter(adata_KO,  x='latent_time', y=Repression_driver_gene, layer="Mu", ncols=6, color="celltype", frameon=True, alpha=.9, size=60,
               show=True, fontsize=14, figsize=(3, 1.85), save="./KO_Repression_driver_gene_time_unsplied.pdf", xlabel="Inferred Cell Time", ylabel="")


var_names = Binding_list
var_names = [
    name for name in var_names if name in adata_WT.var_names and name in adata_KO.var_names]
xkey = "Mu"
tkey = "latent_time"
time = adata.obs[tkey].values
time = time[np.isfinite(time)]
X = (
    adata[:, var_names].layers[xkey]
    if xkey in adata.layers.keys()
    else adata[:, var_names].X
)
df = pd.DataFrame(X[np.argsort(time)], columns=var_names)
n_convolve = 30
weights = np.ones(n_convolve) / n_convolve
for gene in var_names:
    try:
        df[gene] = np.convolve(df[gene].values, weights, mode="same")
    except ValueError as e:
        pass  # e.g. all-zero counts or nans cannot be convolved

sort_top_genes = df.columns.tolist()

palette = "viridis"
col_color = 'celltype'
col_colors = to_list(col_color)
col_color = []
for _, col in enumerate(col_colors):
    if not is_categorical(adata, col):
        obs_col = adata.obs[col]
        cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
        adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
        col += "_categorical"
        set_colors_for_categorical_obs(adata, col, palette)
    col_color.append(interpret_colorkey(adata, col)[np.argsort(time)])
col_color_WT = col_color

adata = adata_KO
time = adata.obs[tkey].values
time = time[np.isfinite(time)]
X = (
    adata[:, var_names].layers[xkey]
    if xkey in adata.layers.keys()
    else adata[:, var_names].X
)
df_KO = pd.DataFrame(X[np.argsort(time)], columns=var_names)
weights = np.ones(n_convolve) / n_convolve
for gene in var_names:
    try:
        df_KO[gene] = np.convolve(df_KO[gene].values, weights, mode="same")
    except ValueError as e:
        pass  # e.g. all-zero counts or nans cannot be convolved
max_sort = np.argsort(np.argmax(df_KO.values, axis=0))
df_KO = pd.DataFrame(df_KO.values[:, max_sort],
                     columns=df_KO.columns[max_sort])

palette = "viridis"
col_color = 'celltype'
col_colors = to_list(col_color)
col_color = []
for _, col in enumerate(col_colors):
    if not is_categorical(adata, col):
        obs_col = adata.obs[col]
        cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
        adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
        col += "_categorical"
        set_colors_for_categorical_obs(adata, col, palette)
    col_color.append(interpret_colorkey(adata, col)[np.argsort(time)])

col_color_KO = col_color
cm = sns.clustermap(pd.merge(df.T, df_KO.T.reindex(index=df.T.index), left_index=True, right_index=True),
                    standard_scale=0,
                    cmap="viridis",
                    col_cluster=False,
                    row_cluster=False,
                    yticklabels=True,
                    xticklabels=False,
                    figsize=(6.8, 20),
                    col_colors=np.concatenate(
                        (col_color_WT, col_color_KO), axis=1),
                    cbar_pos=None,
                    dendrogram_ratio=(0, 0)
                    )

plt.savefig("Binding_driver_gene.png", dpi=600)
