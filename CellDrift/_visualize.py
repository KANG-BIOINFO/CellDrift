import pandas as pd 
import numpy as np
import scanpy as sc
from itertools import combinations
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import pearsonr, spearmanr
import sklearn.metrics
import random
from adjustText import adjust_text
from matplotlib import rcParams
import matplotlib.patches as patches
import warnings

# --- default settings ----
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

warnings.filterwarnings("ignore")
sc.settings.figdir = "figures_celldrift"
output_folder = "figures_celldrift/"
sc.settings.set_figure_params(dpi = 120, facecolor = "white")

rcParams['font.family'] = 'Arial'

# create default cmap
def create_cmap(cmap_name, cmap_list):
    cmap_o = plt.get_cmap(cmap_name)
    for i in range(cmap_o.N):
        rgba = cmap_o(i)
        cmap_list.append(matplotlib.colors.rgb2hex(rgba))
    return cmap_list

cmap = []
for name in ["Paired", "tab20b", "tab20c"]:
    cmap = create_cmap(name, cmap)

def plot_gene_QC(adata, min_cells_perGene):
    """
    Plot Quality control process for genes
    """
    y_percentile = np.percentile(adata.var["n_cells"], 90)

    fig, ax = plt.subplots()

    # violin plot
    sns.violinplot(y = adata.var["n_cells"], ax = ax)

    # scatter plot
    n_genes = adata.var.shape[0]
    xa = np.random.normal(0, 0.05, n_genes)
    ya = adata.var["n_cells"]
    ax.scatter(xa, ya, s = 0.05, c = "black")
    ax.set_ylim(0, y_percentile)

    # threshold
    ax.axhline(y = min_cells_perGene, color = "#CC0000", )
    ax.set_title("Number of cells with non-zero expression for genes")
    ax.set_xticks([])
    fig.show()

    filename = (output_folder + "plot_quality_control_genes.pdf")
    fig.savefig(filename, bbox_inches = "tight")


def plot_coefficients(adata, genes, cmap = cmap, max_ci_interval = 5):
    """
    Visualization of direct estimated coefficients for a gene list.

    Parameters
    ----------
    adata
        AnnData
    genes
        a list of genes for visualization
    cmap
        color map (default = "Paired")

    Examples
    --------
    ct.plot_coefficients(adata, ["ISG15", "ISG20", "ISG20L2"])
    """
    df_mu = pd.DataFrame(data = adata.varm["pred_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
    df_sigma = pd.DataFrame(data = adata.varm["pred_sigma"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])

    n_coeff = df_mu.shape[1]
    n_half = int(n_coeff / 2)
    if isinstance(genes, str):
        genes = [genes]
    
    # mean
    fig, ax = plt.subplots()

    for gene, color in zip(genes, cmap):
        mean_arr = df_mu.loc[gene,:]
        # errorbar
        conf_int = np.clip(1.96 * df_sigma.loc[gene,:], 0, max_ci_interval)
        ax.errorbar(x = np.arange(n_coeff)[:n_half], y = mean_arr[:n_half], yerr = conf_int[:n_half], c = color)
        ax.errorbar(x = np.arange(n_coeff)[n_half:], y = mean_arr[n_half:], yerr = conf_int[n_half:], c = color)

        ax.plot(np.arange(n_coeff)[:n_half], mean_arr[:n_half], marker = 'o', linestyle = "", c = color, label = gene)
        ax.plot(np.arange(n_coeff)[n_half:], mean_arr[n_half:], marker = 'v', linestyle = "", c = color)

    ax.axhline(y=0.0, color='#E0E0E0', linestyle='--',)
    ax.set_xticks(np.arange(df_mu.shape[1]))
    ax.set_xticklabels(df_mu.columns, rotation = 45, ha = "right", fontsize = 8)
    ax.set_xlabel("Covariate Name")
    ax.set_ylabel("Coefficients")
    ax.set_title("GLM coefficients for " + gene) if len(genes) == 1 else ax.set_title("GLM coefficients for " + genes[0] + ", etc.")
    ax.grid(False)
    plt.legend()

    filename = ("coefficients_" + gene + ".pdf") if len(genes) == 1 else ("coefficients_" + genes[0] + " etc.pdf")
    fig.savefig(output_folder + filename, bbox_inches = "tight")
    plt.draw()
    plt.pause(0.001)


def plot_pairwise_coefficients(adata, genes, cmap = cmap, ci_min = -30, ci_max = 30):
    """
    Visualization of estimated pairwise coefficients for a gene list.

    Parameters
    ----------
    adata
        AnnData
    genes
        a list of genes for visualization
    cmap
        color map (default = "Paired")
    max_ci_interval
        Longest confidence interval bar, default is 5
    Examples
    --------
    ct.plot_coefficients(adata, ["ISG15", "ISG20", "ISG20L2"])
    """
    df_mu = pd.DataFrame(data = adata.varm["multcomp_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).copy()
    df_pval = pd.DataFrame(data = adata.varm["multcomp_pvals"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).copy()
    df_lci = pd.DataFrame(data = adata.varm["multcomp_lci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).copy()
    df_uci = pd.DataFrame(data = adata.varm["multcomp_uci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).copy()
    df_interval = abs(df_uci - df_lci) / 2

    n_coeff = df_mu.shape[1]
    if isinstance(genes, str):
        genes = [genes]
    
    # mean
    fig, ax = plt.subplots()

    for gene, color in zip(genes, cmap):

        mean_arr = df_mu.loc[gene, :]
        lci_arr = df_lci.loc[gene,:]
        uci_arr = df_uci.loc[gene,:]
        pval = df_mu.loc[gene, :]

        # errorbar
        conf_int = df_interval.loc[gene, :]
        noise = np.clip(np.random.normal(0, 0.05, n_coeff), -0.03, 0.03)

        # remove nodes with extreme values
        extreme_idx = np.unique(np.append(np.arange(n_coeff)[uci_arr < ci_min], np.arange(n_coeff)[lci_arr > ci_max]))
        mean_arr[extreme_idx] = np.nan
        conf_int[extreme_idx] = np.nan

        ax.errorbar(x = np.arange(n_coeff) + noise, y = mean_arr, yerr = conf_int, c = color, elinewidth = 0.75)
        ax.plot(np.arange(n_coeff), mean_arr, marker = 'o', linestyle = "", c = color, label = gene)

    ax.axhline(y=0.0, color='#E0E0E0', linestyle='--',)
    ax.set_xticks(np.arange(df_mu.shape[1]))
    ax.set_xticklabels(df_mu.columns, rotation = 45, ha = "right", fontsize = 8)
    ax.set_xlabel("Covariate Name")
    ax.set_ylabel("Coefficients")
    ax.set_title("GLM coefficients for " + gene) if len(genes) == 1 else ax.set_title("GLM coefficients for " + genes[0] + ", etc.")

    ax.grid(False)

    # plt.legend()
    # ax.legend(bbox_to_anchor=(1.1, 1.05))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    filename = ("pairwise_coefficients_" + gene + ".pdf") if len(genes) == 1 else ("pairwise_coefficients_" + genes[0] + " etc.pdf")
    fig.savefig(output_folder + filename, bbox_inches = "tight")
    plt.draw()
    plt.pause(0.001)


def plot_correlation(adata, genes, n_permute = 20, figsize = (10,4.5)):
    """
    Plot correlation between given genes.
    """
    df_mu = pd.DataFrame(data = adata.varm["pred_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
    df_sigma = pd.DataFrame(data = adata.varm["pred_sigma"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])

    n_coeff = df_mu.shape[1]
    n_half = int(n_coeff / 2)
    if isinstance(genes, str):
        genes = [genes]
    
    # compute ensemble correlation (cell type / perturbation separately)
    corr_arr_cell = np.zeros(shape = (n_permute, len(genes), len(genes)))
    corr_arr_pert = np.zeros(shape = (n_permute, len(genes), len(genes)))
                             
    for idx in range(n_permute):
        sampled_vals = np.random.normal(loc = df_mu.loc[genes,:], scale = df_sigma.loc[genes,:])
        sampled_vals_cell = sampled_vals[:, : n_half]
        sampled_vals_pert = sampled_vals[:, n_half : ]
        for j in range(len(genes)):
            for k in range(j, len(genes)):
                corr_arr_cell[idx,j,k] = corr_arr_cell[idx, k, j] = pearsonr(sampled_vals_cell[j,:], sampled_vals_cell[k,:])[0]
                corr_arr_pert[idx,j,k] = corr_arr_pert[idx, k, j] = pearsonr(sampled_vals_pert[j,:], sampled_vals_pert[k,:])[0]

    corr_arr_cell = np.median(corr_arr_cell, axis = 0)
    corr_arr_pert = np.median(corr_arr_pert, axis = 0)
    
    #sns.clustermap(corr_arr_cell, cmap = 'bwr', figsize = figsize, xticklabels = genes, yticklabels = genes)
    #sns.clustermap(corr_arr_pert, cmap = 'bwr', figsize = figsize, xticklabels = genes, yticklabels = genes)
    fig, (ax1, ax2) = plt.subplots(figsize = figsize, nrows = 1, ncols = 2)
    sns.heatmap(corr_arr_cell, cmap = 'bwr', xticklabels = genes, yticklabels = genes, ax = ax1)
    sns.heatmap(corr_arr_pert, cmap = 'bwr', xticklabels = genes, yticklabels = genes, ax = ax2)

    filename = ("correlation_of_coeff_" + gene + ".pdf") if len(genes) == 1 else ("correlation_of_coeff_" + genes[0] + " etc.pdf")
    output_filename = adata.uns["celldrift"]["figure_dir"] + filename
    fig.savefig(output_filename, bbox_inches = "tight")


def plot_top_perturbed_genes(adata, cell_type, n_genes = 10, p_threshold = 0.05):
    """
    scatter plots for coefficients in specific cell type (or perturbation) vs. the reference cell type (or perturbation).
    Parameters
    ----------
    adata
        anndata
    cell_type
        selected cell type for plotting pairwise coefficients
    n_genes
        number of genes to show on each side (up and down genes)
    """
    coeff = cell_type + "_" + adata.uns["celldrift"]["group_perturb"] + "-" + cell_type + "_" + adata.uns["celldrift"]["group_control"]

    # initialization of 3 input data frames
    df_pvals = pd.DataFrame(data = adata.varm["multcomp_pvals"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
    sign_genes = df_pvals.loc[df_pvals[coeff] < p_threshold, :].index.values

    df_mu = pd.DataFrame(data = adata.varm["multcomp_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).loc[sign_genes,:]
    df_lci = pd.DataFrame(data = adata.varm["multcomp_lci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).loc[sign_genes,:]
    df_uci = pd.DataFrame(data = adata.varm["multcomp_uci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"]).loc[sign_genes,:]
    df_interval = (df_uci - df_lci) / 2

    # order and select top / bottom gene
    df_mu = df_mu.sort_values(by = coeff, axis = 0, ascending = False)
    df_interval = df_interval.loc[list(df_mu.index.values),:]
    df_pvals = df_pvals.loc[list(df_mu.index.values),:]
    
    # retrieve information
    up_genes = df_mu.index.values[:n_genes]
    down_genes = df_mu.index.values[-n_genes:]
    
    up_params = df_mu[coeff][:n_genes]
    down_params = df_mu[coeff][-n_genes:]
    
    up_err = df_interval[coeff][:n_genes]
    down_err = df_interval[coeff][-n_genes:]
    
    up_pvals = df_pvals[coeff][:n_genes]
    down_pvals = df_pvals[coeff][-n_genes:]
    
    # draw the plot
    fig, ax = plt.subplots(figsize = (9,4.5))
    ax.errorbar(x = np.arange(n_genes * 2), 
                y = np.concatenate((up_params, down_params)), 
                yerr = list(up_err.values)+ list(down_err.values),
                linestyle='',
                ecolor = ["#FD8484"] * n_genes + ["#3399FF"] * n_genes)
    
    # three cases
    ax.scatter(x = np.arange(n_genes * 2), 
               y = np.concatenate((up_params, down_params)),
               marker = "D",
               s = 80,
               facecolors = ["#FD8484"] * n_genes + ["#3399FF"] * n_genes)
    
    ax.axhline(y=0.5, color='#004C99', linestyle='--')
    
    ax.set_xticklabels(labels = list(up_genes) + list(down_genes), rotation = 45, ha = "right", fontsize = 15)
    ax.tick_params(axis='y', which='major', labelsize=15)
    
    ax.set_xticks(np.arange(n_genes*2))
    
    ax.set_xlabel("Genes", fontsize = 20)
    ax.set_ylabel("diff coefficients", fontsize = 20)
    
    filename = "Top_perturbed_genes_" + cell_type + ".pdf"
    output_filename = output_folder + filename
    fig.savefig(output_filename, bbox_inches = "tight")

def categorize_genes_pairwiseCoeff(df_mu, df_pvals, genes, pairwise = True, p_threshold = 0.05):
    """
    Separate genes into different categories: non-significant in both coefficients; significant only in one coefficient (A or B); significant in both.
    """
    n_genes = len(genes)
    
    # create a coefficient matrix
    df_type = pd.DataFrame(index = genes)
    df_mu = df_mu.loc[genes, :]
    df_pvals = df_pvals.loc[genes, :]
    
    # categorize
    significance = [list(df_pvals.iloc[i,:] < p_threshold) for i in range(n_genes)]
    categories = []
    for term in significance:
        if term == [True, True]:
            categories.append("Both")
        elif term == [True, False]:
            categories.append("Coeff 1 only")
        elif term == [False, True]:
            categories.append("Coeff 2 only")
        else:
            categories.append("Neither")

    df_type["category"] = categories
    return df_type


def plot_coefficient_comparison(adata, coeffs, pairwise = True, errorbar = True, genes = "All", genes_highlight = None, arrows = False, figsize = (14,6)):
    """
    Parameters
    ----------
    adata
        AnnData
    coeffs
        A pair of coefficients to be visualized
    pairwise
        Whether coefficients to compare are from direct GLM estimated coefficients or from emmeans pooled pairwise coefficients
    genes
        Selection of genes
    genes_highlight
        Genes to highlight with text labels
    figsize
        Figure size
    """    
    if pairwise:
        df_mu = pd.DataFrame(data = adata.varm["multcomp_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
        df_pvals = pd.DataFrame(data = adata.varm["multcomp_pvals"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
        df_lci = pd.DataFrame(data = adata.varm["multcomp_lci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
        df_uci = pd.DataFrame(data = adata.varm["multcomp_uci"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
        df_interval = (df_uci - df_lci) / 2

    else:
        df_mu = pd.DataFrame(data = adata.varm["pred_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
        df_sigma = pd.DataFrame(data = adata.varm["pred_sigma"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
        df_pvals = pd.DataFrame(data = adata.varm["pred_pvals"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
    
    df_mu = df_mu[coeffs]
    df_pvals = df_pvals[coeffs]
    if genes == "All":
        # genes = df_mu.index.values
        genes = df_mu.loc[(np.max(abs(df_mu), axis = 1)) < 10, :].index.values # apply a hard criteria
        df_mu = df_mu.loc[genes, :]
        df_pvals = df_pvals.loc[genes, :]
        if pairwise:
            df_interval = df_interval.loc[genes, :]
        else:
            df_sigma = df_sigma.loc[genes, :]
    elif isinstance(genes, str): # single gene
        genes = [genes]
    
    coeffs_1 = coeffs[0]
    coeffs_2 = coeffs[1]
    
    # separate genes into 4 categories
    df_type = categorize_genes_pairwiseCoeff(df_mu, df_pvals, genes = genes, pairwise = pairwise)
    categories = ["Neither", "Coeff 1 only", "Coeff 2 only", "Both"]
    colors = ["#A0A0A0", "#66CC00", "#3399FF", "#FF6666"]
    
    fig, (ax1, ax2) = plt.subplots(figsize = figsize, nrows = 1, ncols = 2)#, sharex = True, sharey = True)
    
    for category, color in zip(categories, colors):
        df_sub = df_type.loc[df_type["category"] == category, :]
        genes_sub = df_sub.index.values
        
        x = df_mu.loc[genes_sub, coeffs_1]
        y = df_mu.loc[genes_sub, coeffs_2]

        ax1.scatter(x, y, c = color, s = 2)
        if errorbar:
            if pairwise:
                ax1.errorbar(x, y, xerr = df_interval.loc[genes_sub, coeffs_1], linestyle = "", 
                            color = color, elinewidth = 0.15)
                ax1.errorbar(x, y, yerr = df_interval.loc[genes_sub, coeffs_2], linestyle = "", 
                            color = color, elinewidth = 0.15) 
            else:
                ax1.errorbar(x, y, xerr = 1.96 * df_sigma.loc[genes_sub, coeffs_1], linestyle = "", 
                            color = color, elinewidth = 0.15)
                ax1.errorbar(x, y, yerr = 1.96 * df_sigma.loc[genes_sub, coeffs_2], linestyle = "", 
                            color = color, elinewidth = 0.15) 
        
    ax1.axhline(y = 0, color = "#808080", linestyle = "-", linewidth = .2)
    ax1.axvline(x = 0, color = "#808080", linestyle = "-", linewidth = .2)
        
    ax1.set_xlabel(coeffs_1, weight = "bold", size = 13)
    ax1.set_ylabel(coeffs_2, weight = "bold", size = 13)
    
    # set xlim and ylim
    x = np.array(df_mu[coeffs_1])
    y = np.array(df_mu[coeffs_2])

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    ax1.set_xlim(xmin - 0.1, xmax + 0.1)
    ax1.set_ylim(ymin - 0.1, ymax + 0.1)
    df_coor = pd.DataFrame(data = np.array([x, y]).T, index = genes, columns = ["x", "y"])

    # build a linear model between coeff1 and 2
    X = x.reshape(-1, 1)
    reg = sklearn.linear_model.LinearRegression().fit(X, y)
    y_pred = reg.predict(X)

    # add r2 score
    r2_score = np.round(sklearn.metrics.r2_score(y, y_pred), 3)
    ax1.text(x = 0.95 * xmin + 0.05 * xmax, y = 0.95 * ymax + 0.05 * ymin, s = "R$^{2}$ = " + str(r2_score),
            fontsize = 12, weight = "bold", color = "#000066")
    
    # randomly highlight genes (2 + 2 + 4 + 4)
    if genes_highlight == None:
        genes_coeff1 = random.sample(list(df_type.loc[df_type["category"] == "Coeff 1 only", :].index.values), 2)
        genes_coeff2 = random.sample(list(df_type.loc[df_type["category"] == "Coeff 2 only", :].index.values), 2)
        genes_both = random.sample(list(df_type.loc[df_type["category"] == "Both", :].index.values), 4)
        genes_highlight = genes_coeff1 + genes_coeff2 + genes_both
    
    texts = [ax1.text(x = df_coor.loc[gene, "x"], y = df_coor.loc[gene, "y"], s = gene, fontsize = 10, fontfamily = "Helvetica") for gene in genes_highlight]
    if arrows:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax = ax1)
    else:
        adjust_text(texts, ax = ax1)

    ### add contour plot    
    # perform kde
    xmin, xmax = np.percentile(x, 1), np.percentile(x, 99) # set x lim and y lim b2y percentile
    ymin, ymax = np.percentile(y, 1), np.percentile(y, 99)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    #ax2.set_xlim(0.8 * xmin + 0.2 * xmax, 0.8 * xmax + 0.2 * xmin)
    #ax2.set_ylim(0.8 * ymin + 0.2 * ymax, 0.8 * ymax + 0.2 * ymin)
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlabel(coeffs_1, weight = "bold", size = 13)
    ax2.axhline(y = 0, color = "#808080", linestyle = "-", linewidth = .2)
    ax2.axvline(x = 0, color = "#808080", linestyle = "-", linewidth = .2)
    
    # Create a Rectangle patch
    rect = patches.Rectangle((xmin, ymin), xmax-xmin, ymax-ymin, linewidth=0.5, edgecolor=None, facecolor='#66B2FF', alpha = 0.08, fill = True)
    # Add the patch to the Axes
    ax1.add_patch(rect)
    
    # contour plot
    cfset = ax2.contourf(xx, yy, f, cmap='Blues',alpha = 1)
    ## Or kernel density estimate plot instead of the contourf plot
    #ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
    # Contour plot
    cset = ax2.contour(xx, yy, f, colors='k')
    # Label plot
    ax2.clabel(cset, inline=1, fontsize=10)

    filename = "coefficients_distributions" + coeffs_1 + "_" + coeffs_2 + ".pdf"
    output_filename = output_folder + filename
    fig.savefig(output_filename, bbox_inches = "tight")
    

def plot_patternScore(adata, gene_set, reduction = "umap", vmax = 10, vmin = 0):
    """
    Draw accumulated expression values for genes in a gene pattern on reduced dimensions.
    """
    pert_key = adata.uns["celldrift"]["perturb_key"]
    ctrl_level = adata.uns["celldrift"]["group_control"]
    pert_level = adata.uns["celldrift"]["group_perturb"]

    gene_set_name = gene_set[0] + " et al"
    if type(adata.X) == np.ndarray:
        adata.obs["pattern_score"] = np.mean(adata[:, gene_set].X, axis = 1)
    else:
        adata.obs["pattern_score"] = np.mean(adata[:, gene_set].X.toarray(), axis = 1)

    sc.pl.umap(adata[adata.obs[pert_key] == ctrl_level], color = "pattern_score", vmax = vmax, vmin = vmin, 
                use_raw = False, show = False, save = "umap_patternScore_" + ctrl_level + "_" + gene_set_name + ".pdf")
    sc.pl.umap(adata[adata.obs[pert_key] == pert_level], color = "pattern_score", vmax = vmax, vmin = vmin,
                use_raw = False, show = False, save = "umap_patternScore_" + pert_level + "_" + gene_set_name + ".pdf")


def plot_patternScore_violinplot(adata, gene_set):
    """
    plot gene pattern score using violin plot (revision required)
    """
    meta = adata.obs.copy()

    fig, ax = plt.subplots()
    sns.violinplot(x = "cell_stim", y = "pattern_score", data = meta, ax = ax,
                scale = "width")
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = "right")
    ax.set_ylim(0,15)
    ax.set_xlabel("cell_stim", fontsize = 14)
    ax.set_ylabel("pattern_score", fontsize = 14)
    fig.suptitle("gene pattern score", fontsize = 14)

    plt.show()
    return 1

def create_linkage_matrix(model):
    """
    Create linkage matrix for dendrogram generation.
    """
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    return linkage_matrix


def get_dendrogram_orders(model, **kwargs):
    """
    get x labels of a dendrogram
    """
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    linkage_matrix = create_linkage_matrix(model)

    # Plot the corresponding dendrogram
    fig, ax = plt.subplots()
    dendrogram(linkage_matrix, ax = ax, **kwargs)
    x_labels = ax.get_xlabel()
    plt.close()

    return ax


def plot_dendrogram(model, ax, **kwargs):
    """
    Plot a dendrogram with ax input 
    """
    # Create linkage matrix and then plot the dendrogram
    linkage_matrix = create_linkage_matrix(model)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, ax = ax, **kwargs)

    return ax


def plot_cell_type_correlation_raw(adata):
    """
    plot correlation matrix of cell type coefficients and perturbation coefficients using direct estimation of GLM.
    """

    df_mu = pd.DataFrame(data = adata.varm["pred_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["coefficient_names"])
    n_coeffs_side = int(df_mu.shape[1] / 2)
    df_mu_cell = df_mu.iloc[:, :n_coeffs_side]
    df_mu_pert = df_mu.iloc[:, n_coeffs_side:]

    df_mu_cell.rename({"intercept": "B cells"}, axis = 1, inplace = True)
    cell_types = df_mu_cell.columns
    n_coeffs = len(cell_types)

    ### -- cell type correlation -- ###
    corr_cell = np.eye(n_coeffs, n_coeffs)
    for idx_i, idx_j in combinations(range(n_coeffs), 2):
        corr_cell[idx_i, idx_j] = corr_cell[idx_j, idx_i] = pearsonr(df_mu_cell.iloc[:, idx_i], df_mu_cell.iloc[:, idx_j])[0]
    df_corr_cell = pd.DataFrame(data = corr_cell, index = cell_types, columns = cell_types)

    # get the dendrogram (cell type)
    model_c = AgglomerativeClustering(n_clusters = None, distance_threshold = 0)
    model_c = model_c.fit(df_corr_cell)
    ax_dendro = get_dendrogram_orders(model_c) # dendrogram
    orders_c = [int(item.get_text()) for item in ax_dendro.get_xticklabels()]
    df_corr_cell = df_corr_cell.iloc[orders_c, orders_c]

    ### -- perturbation correlation -- ###
    corr_pert = np.eye(n_coeffs, n_coeffs)
    for idx_i, idx_j in combinations(range(n_coeffs), 2):
        corr_pert[idx_i, idx_j] = corr_pert[idx_j, idx_i] = pearsonr(df_mu_pert.iloc[:, idx_i], df_mu_pert.iloc[:, idx_j])[0]
    df_corr_pert = pd.DataFrame(data = corr_pert, index = cell_types, columns = cell_types)
    
    # get the dendrogram (perturbation)
    model_p = AgglomerativeClustering(n_clusters = None, distance_threshold = 0)
    model_p = model_p.fit(df_corr_pert)
    ax_dendro = get_dendrogram_orders(model_p) # dendrogram
    orders_p = [int(item.get_text()) for item in ax_dendro.get_xticklabels()]
    df_corr_pert = df_corr_pert.iloc[orders_c, orders_p]
    
    # vmin and vmax
    vmin = min(np.min(corr_cell), np.min(corr_pert))
    vmax = max(np.max(corr_cell), np.max(corr_pert))

    # initial figure
    fig = plt.figure(figsize = (13.5,9))
    gs = GridSpec(nrows = 2, ncols = 2, height_ratios=(0.15,0.85), width_ratios = (0.50, 0.50))
    gs.update(wspace = 0.075, hspace = 0.025)
    
    # dendrogram upleft
    ax0 = fig.add_subplot(gs[0,0])
    plot_dendrogram(model_c, ax0)
    ax0.xaxis.set_ticklabels([])

    # dendrogram upright
    ax1 = fig.add_subplot(gs[0,1])
    plot_dendrogram(model_p, ax1)
    ax1.xaxis.set_ticklabels([])
    
    # heatmap 1
    ax3 = fig.add_subplot(gs[1, 0])
    sns.heatmap(df_corr_cell, cmap = "bwr", ax = ax3, vmin = vmin, vmax = vmax, cbar = False)
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation = 45, ha = "right")
    ax3.set_xlabel("correlation - cell type coefficients", fontsize = 16)
    
    # heatmap 2
    ax4 = fig.add_subplot(gs[1, 1])
    sns.heatmap(df_corr_pert, cmap = "bwr", ax = ax4, vmin = vmin, vmax = vmax, cbar = False)
    ax4.yaxis.set_ticklabels([])
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation = 45, ha = "right")
    ax4.set_xlabel("correlation - perturbation coefficients", fontsize = 16)
    
    filename = "heatmap_correlation_of_correlations.pdf"
    output_filename = adata.uns["figure_dir"] + filename
    fig.savefig(output_filename, bbox_inches = "tight")


def plot_cell_type_correlation(adata, method = "pearson", vmin = None, vmax = None):
    """
    plot correlation matrix of pairwise coefficients between cell populations.
    """
    if method not in ["pearson", "spearman"]:
        raise Exception("Correlation method should be one of pearson and spearman.")

    df_mu = pd.DataFrame(data = adata.varm["multcomp_mu"], index = adata.var_names, columns = adata.uns["celldrift"]["multicomp_terms"])
    n_coeffs_side = int(df_mu.shape[1])

    cell_types = [i.split("_")[0] for i in df_mu.columns]
    n_coeffs = len(cell_types)

    ### -- cell type correlation -- ###
    corr_cell = np.eye(n_coeffs, n_coeffs)
    for idx_i, idx_j in combinations(range(n_coeffs), 2):
        if method == "pearson":
            corr_cell[idx_i, idx_j] = corr_cell[idx_j, idx_i] = pearsonr(df_mu.iloc[:, idx_i], df_mu.iloc[:, idx_j])[0]
        elif method == "spearman":
            corr_cell[idx_i, idx_j] = corr_cell[idx_j, idx_i] = spearmanr(df_mu.iloc[:, idx_i], df_mu.iloc[:, idx_j])[0]
    df_corr_cell = pd.DataFrame(data = corr_cell, index = cell_types, columns = cell_types)

    # get the dendrogram (cell type)
    model_c = AgglomerativeClustering(n_clusters = None, distance_threshold = 0)
    model_c = model_c.fit(df_corr_cell)
    ax_dendro = get_dendrogram_orders(model_c) # dendrogram
    orders_c = [int(item.get_text()) for item in ax_dendro.get_xticklabels()]
    df_corr_cell = df_corr_cell.iloc[orders_c, orders_c]

    # vmin and vmax
    vmin = np.min(corr_cell) if vmin == None else vmin
    vmax = np.max(corr_cell) if vmax == None else vmax

    fig = sns.clustermap(df_corr_cell, cmap = "bwr", vmin = vmin, vmax = vmax)
    filename = "heatmap_correlation_of_correlations.pdf"
    output_filename = output_folder + filename
    fig.savefig(output_filename, bbox_inches = "tight")

'''
visualization for factor analysis
'''
def draw_components_scatter(model, X_transformed, observation_names, output_name, figsize = (12, 5.5)):
    
    fig, axes = plt.subplots(ncols = 2, figsize = figsize)
    
    axes[0].scatter(x = X_transformed[:, 0], y = X_transformed[:, 1])#, c = data.target)
    axes[0].set_title('components for samples')
    for idx, text in enumerate(observation_names):
        axes[0].text(x = X_transformed[idx, 0], y = X_transformed[idx, 1], s = text, fontsize = 4)
    
    axes[1].scatter(x = model.components_[0, :], y = model.components_[1, :], c = '#606060')
    axes[1].set_title('components for features')

    # for feature_idx, feature_name in enumerate(feature_names):
    #     axes[1].text(x = model.components_[0, feature_idx], y = model.components_[1, feature_idx], s = feature_name)

    fig.savefig('figures_celldrift/' + output_name + '.png')

def draw_correlation(corr_mtx, output_name, feature_names, figsize = (6,6)):
    df_corr = pd.DataFrame(data = corr_mtx,
                            index = feature_names,
                            columns = feature_names)

    fig = sns.clustermap(df_corr, cmap = 'bwr', figsize = figsize)
    fig.savefig('figures_celldrift/' + output_name + '.png')

def draw_components_heatmap(model, X_transformed, output_name, figsize = (12, 5.5)):
    fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = figsize)
    sns.heatmap(X_transformed, cmap = 'bwr', ax = ax1)
    sns.heatmap(model.components_.T, cmap = 'bwr', ax = ax2)
    fig.savefig('figures_celldrift/' + output_name + '.png')

def draw_components_clustermap(model, X_transformed, observation_names, feature_names, output_name, figsize = (12, 5.5)):
    # (perturbation + cell type) x factor 
    df1 = pd.DataFrame(data = X_transformed, index = observation_names, columns = [('Factor_' + str(i)) for i in range(X_transformed.shape[1])])
    g1 = sns.clustermap(df1, cmap = 'bwr')
    g1.savefig('figures_celldrift/' + output_name + '_pert_type_x_factor.png')

    # feature x factor
    df2 = pd.DataFrame(data = model.components_.T, index = feature_names, columns = [('Factor_' + str(i)) for i in range(X_transformed.shape[1])])
    g2 = sns.clustermap(df2, cmap = 'bwr')
    g2.savefig('figures_celldrift/' + output_name + '_feature_x_factor.png')
