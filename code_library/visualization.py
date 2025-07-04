import numpy as np
import scipy
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gc
import colorsys
import seaborn as sns
import scanpy as sc

gene_highlight_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'gene_highlight_cmap', [
        colorsys.hsv_to_rgb(0.8, 0, 0.9),
        colorsys.hsv_to_rgb(0.8, 0.15, 0.75),
        colorsys.hsv_to_rgb(0.8, 0.3, 0.6),
        colorsys.hsv_to_rgb(0.8, 0.45, 0.45),
        colorsys.hsv_to_rgb(0.8, 0.6, 0.3),
        colorsys.hsv_to_rgb(0.8, 0.75, 0.15),
        colorsys.hsv_to_rgb(0.8, 0.9, 0),
        colorsys.hsv_to_rgb(0.75, 0.91, 0.2),
        colorsys.hsv_to_rgb(0.7, 0.92, 0.35),
        colorsys.hsv_to_rgb(0.65, 0.93, 0.45),
        colorsys.hsv_to_rgb(0.6, 0.94, 0.55),
        colorsys.hsv_to_rgb(0.55, 0.95, 0.6),
        colorsys.hsv_to_rgb(0.5, 0.96, 0.65),
        colorsys.hsv_to_rgb(0.45, 0.97, 0.7),
        colorsys.hsv_to_rgb(0.4, 0.98, 0.75),
        colorsys.hsv_to_rgb(0.35, 0.99, 0.8),
        colorsys.hsv_to_rgb(0.3, 1, 0.85),
        colorsys.hsv_to_rgb(0.25, 1, 0.9),
        colorsys.hsv_to_rgb(0.2, 1, 0.95),
        colorsys.hsv_to_rgb(0.15, 1, 1)
    ])
    

def str_int(s):
    try:
        _ = np.array(s, dtype=int)
        return True
    except:
        return False


def preprocessing_categories(adata, this_factor, mask=None):
    used_adata = adata[:, [0]].copy()
    if mask is None:
        all_categories = np.unique(used_adata.obs[this_factor].values.astype(str))
    else:
        all_categories = np.unique(used_adata[mask, :].obs[this_factor].values.astype(str))
    if (str_int(all_categories)):
        all_categories = all_categories[np.argsort(all_categories.astype(int))]
    else:
        all_categories = all_categories[np.argsort(all_categories)]
    return used_adata, all_categories


def clustering_plot(adata, clustering, basis="X_umap", mask=None, value=None, **kwargs):
    tmp, all_categories = preprocessing_categories(adata, clustering, mask=mask)
    for ii in all_categories:
        if mask is None:
            this_mask = tmp.obs[ii] = tmp.obs[clustering].astype(str).isin([ii])
        else:
            this_mask = np.logical_and(tmp.obs[clustering].astype(str).isin([ii]), mask)
        if value is None:
            tmp.obs[ii] = this_mask.astype(int)
        else:
            tmp.obs[ii] = np.where(this_mask, tmp.obs[value], np.nan)
    sc.pl.embedding(tmp,
                    basis,
                    color=[x for x in all_categories],
                    sort_order=kwargs.pop("sort_order", True),
                    color_map=kwargs.pop("color_map", matplotlib.colors.LinearSegmentedColormap.from_list("", ["#e0e0e0", "#FF0000"]) if value is None else None),
                    frameon=kwargs.pop("frameon", False),
                    colorbar_loc=kwargs.pop("colorbar_loc", None),
                    **kwargs)
    del tmp, all_categories
    plt.close('all')
    gc.collect()


def factor_distribution(adata, plot_factor, x_axis, path=None, ncol=None):
    tmp, x_label = preprocessing_categories(adata, x_axis)
    num_plots = len(plot_factor)
    if ncol is None:
        ncol = np.ceil(np.sqrt(num_plots)).astype(int)
    nrow = np.ceil(num_plots / ncol).astype(int)
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(len(x_label) * ncol * 0.25 + 0.5, 1 + 3 * nrow))
    this_col = 0
    this_row = 0
    plt.rcParams["axes.grid"] = False
    for index, ii in enumerate(plot_factor):
        _, y_stack = preprocessing_categories(tmp, ii)
        this_df = pd.DataFrame(
            [[np.sum(tmp.obs[ii][tmp.obs[x_axis].astype(str) == x].astype(str) == y) / tmp.obs[ii][tmp.obs[x_axis].astype(str) == x].size for y in y_stack]
             for x in x_label], index=x_label, columns=y_stack)
        if np.logical_and(nrow > 1, ncol > 1):
            this_ax=axs[this_row, this_col]
        elif nrow > 1:
            this_ax=axs[this_row]
        elif ncol > 1:
            this_ax=axs[this_col]
        else:
            this_ax=axs
        this_df.plot(kind="bar", stacked=True, title=ii, ax=this_ax)
        this_ax.legend(bbox_to_anchor=(1.02, 1))
        if this_col < (ncol - 1):
            this_col += 1
        else:
            this_row += 1
            this_col = 0
    fig.subplots_adjust(hspace=0.3)
    fig.show()
    if path is not None:
        fig.savefig(path)


def factor_distribution_per_method(adata, plot_factor, x_axis, precent="y", path=None, ncol=None):
    tmp, x_label = preprocessing_categories(adata, x_axis)
    num_plots = len(plot_factor)
    if ncol is None:
        ncol = np.ceil(np.sqrt(num_plots)).astype(int)
    nrow = np.ceil(num_plots / ncol).astype(int)
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(len(x_label) * ncol * 0.25, 1 + 3 * nrow))
    this_col = 0
    this_row = 0
    plt.rcParams["axes.grid"] = False
    if precent == "y":
        denominator_str = "np.sum(tmp.obs[ii].astype(str) == y)"
    elif precent == "x":
        denominator_str = "np.sum(tmp.obs[x_axis].astype(str) == x)"
    elif precent == "no":
        denominator_str = "1.0"
    else:
        raise NameError('"precent" should be one of "y", "x" or "no"!')
    for index, ii in enumerate(plot_factor):
        _, y_stack = preprocessing_categories(tmp, ii)
        this_df = pd.DataFrame(
            [[np.sum(tmp.obs[ii][tmp.obs[x_axis].astype(str) == x].astype(str) == y) / eval(denominator_str) for y in y_stack] for x in
             x_label], index=x_label, columns=y_stack)
        if np.logical_and(nrow > 1, ncol > 1):
            this_df.plot(kind="bar", stacked=False, title=ii, ax=axs[this_row, this_col])
        elif nrow > 1:
            this_df.plot(kind="bar", stacked=False, title=ii, ax=axs[this_row])
        elif ncol > 1:
            this_df.plot(kind="bar", stacked=False, title=ii, ax=axs[this_col])
        else:
            this_df.plot(kind="bar", stacked=False, title=ii, ax=axs)
        if this_col < (ncol - 1):
            this_col += 1
        else:
            this_row += 1
            this_col = 0
    plt.tight_layout()
    fig.show()
    if path is not None:
        fig.savefig(path)


def distribution_matrix(this_df, title=None, path=None):
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(1 + this_df.shape[0] * 0.5, 1 + 0.5 * this_df.shape[1]))
    plt.rcParams["axes.grid"] = False
    this_df.plot(kind="bar", stacked=False, title=title, ax=axs)
    plt.tight_layout()
    fig.show()
    if path is not None:
        fig.savefig(path)

def _wedge_dir(vertices, direction):
    """
    Args:
      <vertices>  The vertices from matplotlib.collections.PolyCollection
      <direction> Direction must be 'horizontal' or 'vertical' according to how
                   your plot is laid out.
    Returns:
      - a string in ['top', 'bottom', 'left', 'right'] that determines where the
         half of the violinplot is relative to the center.
    """
    if direction == 'horizontal':
        result = (direction, len(set(vertices[1:5,1])) == 1)
    elif direction == 'vertical':
        result = (direction, len(set(vertices[-3:-1,0])) == 1)
    outcome_key = {('horizontal', True): 'bottom',
                   ('horizontal', False): 'top',
                   ('vertical', True): 'left',
                   ('vertical', False): 'right'}
    # if the first couple x/y values after the start are the same, it
    #  is the input direction. If not, it is the opposite
    return outcome_key[result]


# https://stackoverflow.com/questions/43357274/separate-halves-of-split-violinplot-to-compare-tail-data
def offset_violinplot_halves(ax, delta, width, inner, direction):
    """
    This function offsets the halves of a violinplot to compare tails
    or to plot something else in between them. This is specifically designed
    for violinplots by Seaborn that use the option `split=True`.

    For lines, this works on the assumption that Seaborn plots everything with
     integers as the center.

    Args:
     <ax>    The axis that contains the violinplots.
     <delta> The amount of space to put between the two halves of the violinplot
     <width> The total width of the violinplot, as passed to sns.violinplot()
     <inner> The type of inner in the seaborn
     <direction> Orientation of violinplot. 'hotizontal' or 'vertical'.

    Returns:
     - NA, modifies the <ax> directly
    """
    # offset stuff
    if inner == 'sticks':
        lines = ax.get_lines()
        for line in lines:
            if direction == 'horizontal':
                data = line.get_ydata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is top, move neg, direction backwards for horizontal
                    data -= delta
                else:
                    # type is bottom, move pos, direction backward for hori
                    data += delta
                line.set_ydata(data)
            elif direction == 'vertical':
                data = line.get_xdata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is left, move neg
                    data -= delta
                else:
                    # type is left, move pos
                    data += delta
                line.set_xdata(data)
    for ii, item in enumerate(ax.collections):
        # axis contains PolyCollections and PathCollections
        if isinstance(item, matplotlib.collections.PolyCollection):
            # get path
            path, = item.get_paths()
            vertices = path.vertices
            half_type = _wedge_dir(vertices, direction)
            # shift x-coordinates of path
            if half_type in ['top','bottom']:
               if inner in ["sticks", None]:
                    if half_type == 'top': # -> up
                        vertices[:,1] -= delta
                    elif half_type == 'bottom': # -> down
                        vertices[:,1] += delta
            elif half_type in ['left', 'right']:
                if inner in ["sticks", None]:
                    if half_type == 'left': # -> left
                        vertices[:,0] -= delta
                    elif half_type == 'right': # -> down
                        vertices[:,0] += delta


def plot_dotplot_with_annotations(adata, var_names=None, groupby=None, dendrogram=False, add_totals=True, n_genes=10,
                                  swap_axes=False,
                                  highlight_genes=None,
                                  highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5},
                                  highlight_cluster_line=False,
                                  cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},
                                  **kwargs):
    """
    Plot a dotplot with additional visual annotations for differential gene expression analysis.

    Args:
    adata (AnnData): AnnData object containing the expression data.
    var_names (list or str): Variable names (genes) to plot or 'deg' to use differentially expressed genes.
    groupby (str): The key for the categorical annotation to group by.
    return_fig (bool): If True, returns the figure object for further customization.
    dendrogram (bool): If True, adds a dendrogram based on the data.
    n_genes (int): Number of top genes to show if `var_names` is 'deg'.

    Returns:
    matplotlib.figure.Figure or None: The dotplot figure if `return_fig` is True, otherwise None.
    """
    # Generate a dotplot based on the type of genes specified in `var_names`
    if var_names == 'deg':
        # Create a dotplot for differentially expressed genes
        fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, groupby=groupby, return_fig=True, show=False, dendrogram=dendrogram, swap_axes=swap_axes, **kwargs)
    else:
        # Create a standard dotplot
        fig = sc.pl.dotplot(adata, var_names=var_names, groupby=groupby, return_fig=True, show=False, dendrogram=dendrogram, swap_axes=swap_axes, **kwargs)
    ax = fig.add_totals().get_axes()['mainplot_ax']
    # Add horizontal spans in the background for better visual separation
    if swap_axes:
        [ax.axvspan(y, y + 1, facecolor='lightgray', alpha=0.3, zorder=0.5) for y in range(1, len(ax.get_xticks()), 2)]
    else:
        [ax.axhspan(y, y + 1, facecolor='lightgray', alpha=0.3, zorder=0.5) for y in range(1, len(ax.get_yticks()), 2)]
    # Add highlighting for the cluster that is currently being visualized
    if highlight_cluster_line:
        if dendrogram:
            cluster_order = adata.uns[f'dendrogram_{groupby}']['categories_ordered']
        else:
            cluster_order = var_names.keys()
        total = ax.get_xlim()[1]
        start = 0 + abs(ax.get_xlim()[0] / total)
        end = 0
        for ix, cl in enumerate(cluster_order):
            span_length = len(var_names[cl]) / total
            end += span_length
            ax.axhspan(ix, ix + 1, xmin=start, xmax=end, **cluster_line_highlight_params)
            start += span_length
    # Add vertical lines to visually separate groups of genes or ranked genes
    plot_lines = True
    if isinstance(var_names, dict):
        # If var_names is a dictionary, add vertical lines after each group of genes
        if dendrogram:
            cluster_order = adata.uns[f'dendrogram_{groupby}']['categories_ordered']
            if set(var_names.keys()) == set(cluster_order):
                cluster_order = adata.uns[f'dendrogram_{groupby}']['categories_ordered']
                annotation_range = np.cumsum([len(var_names[cl]) for cl in cluster_order])[:-1]
            else:
                annotation_range = np.cumsum([len(genes) for genes in var_names.values()])[:-1]
        else:
            annotation_range = np.cumsum([len(genes) for genes in var_names.values()])[:-1]
    elif var_names == 'deg':
        # If plotting differentially expressed genes, add vertical lines every `n_genes` units
        annotation_range = range(n_genes, int(ax.get_xlim()[1]), n_genes)
    else:
        plot_lines = False
    if plot_lines:
        if swap_axes:
            [ax.axhline(i, xmin=0, ls=':', lw=0.5, c='dimgray') for i in annotation_range]
        else:
            [ax.axvline(i, ymin=0, ls=':', lw=0.5, c='dimgray') for i in annotation_range]

    if highlight_genes:
        for label in ax.get_xticklabels():
            if label.get_text() in highlight_genes:
                label.set_bbox(highlight_params)
    fig = plt.gcf()
    return fig
