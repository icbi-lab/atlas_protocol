# Libraries for visualization
import warnings
from collections.abc import Sequence
from itertools import zip_longest
from math import ceil

import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns
import statsmodels.stats.multitest


def plot_paired(
    adata,
    groupby,
    *,
    paired_by=None,
    var_names=None,
    show=True,
    return_fig=False,
    n_cols=4,
    panel_size=(3, 4),
    show_legend=False,
    hue=None,
    size=10,
    ylabel="expression",
    pvalues: Sequence[float] = None,
    pvalue_template=lambda x: f"unadj. p={x:.2f}, t-test",
    adjust_fdr=False,
    boxplot_properties=None,
):
    """Pairwise expression plot.Makes on panel with a paired scatterplot for each variable.

    Parameters
    ----------
    adata
        adata matrix (usually pseudobulk).
    groupby
        Column containing the grouping. Must contain exactely two different values.
    paired_by
        Column indicating the pairing (e.g. "patient")
    var_names
        Only plot these variables. Default is to plot all.
    adjust_fdr
        Adjust p-values for multiple testing using the Benjamini-Hochberg procedure.
    boxplot_properties
        Properties to pass to the boxplot function.
    hue
        Column indicating the hue.
    n_cols
        Number of columns in the figure.
    panel_size
        Size of each panel.
    pvalue_template
        Template for the p-value annotation. Must contain a single placeholder for the p-value.
    pvalues
        P-values to annotate. Must be the same length as var_names.
    return_fig
        Return the figure object.
    show
        Show the figure.
    show_legend
        Show the legend.
    size
        Size of the points.
    ylabel
        Label for the y-axis.
    """
    if boxplot_properties is None:
        boxplot_properties = {}
    groups = adata.obs[groupby].unique()
    if len(groups) != 2:
        raise ValueError("The number of groups in the group_by column must be exactely 2")

    if var_names is None:
        var_names = adata.var_names
        if len(var_names) > 20:
            warnings.warn(
                "You are plotting more than 20 variables which may be slow. "
                "Explicitly set the `var_names` paraloeter to turn this off. ",
                stacklevel=2,
            )

    X = adata[:, var_names].X
    try:
        X = X.toarray()
    except AttributeError:
        pass

    groupby_cols = [groupby]
    if paired_by is not None:
        groupby_cols.insert(0, paired_by)
    if hue is not None:
        groupby_cols.insert(0, hue)

    df = adata.obs.loc[:, groupby_cols].join(pd.DataFrame(X, index=adata.obs_names, columns=var_names))

    if paired_by is not None:
        # remove unpaired samples
        df[paired_by] = df[paired_by].astype(str)
        df.set_index(paired_by, inplace=True)
        has_matching_samples = df.groupby(paired_by).apply(lambda x: sorted(x[groupby]) == sorted(groups))
        has_matching_samples = has_matching_samples.index[has_matching_samples].values
        removed_samples = adata.obs[paired_by].nunique() - len(has_matching_samples)
        if removed_samples:
            warnings.warn(f"{removed_samples} unpaired samples removed", stacklevel=2)

        # perform statistics (paired ttest)
        if pvalues is None:
            _, pvalues = scipy.stats.ttest_rel(
                df.loc[
                    df[groupby] == groups[0],
                    var_names,
                ].loc[has_matching_samples, :],
                df.loc[
                    df[groupby] == groups[1],
                    var_names,
                ].loc[has_matching_samples],
            )

        df = df.loc[has_matching_samples, :]
        df.reset_index(drop=False, inplace=True)

    else:
        if pvalues is None:
            _, pvalues = scipy.stats.ttest_ind(
                df.loc[
                    df[groupby] == groups[0],
                    var_names,
                ],
                df.loc[
                    df[groupby] == groups[1],
                    var_names,
                ],
            )

    if adjust_fdr:
        pvalues = statsmodels.stats.multitest.fdrcorrection(pvalues)[1]

    # transform data for seaborn
    df_melt = df.melt(
        id_vars=groupby_cols,
        var_name="var",
        value_name="val",
    )

    # start plotting
    n_panels = len(var_names)
    nrows = ceil(n_panels / n_cols)
    ncols = min(n_cols, n_panels)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * panel_size[0], nrows * panel_size[1]),
        tight_layout=True,
        squeeze=False,
    )
    axes = axes.flatten()
    if hue is None:
        hue = paired_by
    for i, (var, ax) in enumerate(zip_longest(var_names, axes)):
        if var is not None:
            sns.stripplot(
                x=groupby,
                data=df_melt.loc[lambda x: x["var"] == var],  # noqa: B023
                y="val",
                ax=ax,
                hue=hue,
                size=size,
                linewidth=1,
            )
            # sns.lineplot(
            #    x=groupby,
            #    data=df_melt.loc[lambda x: x["var"] == var],
            #    hue=hue,
            #    y="val",
            #    ax=ax,
            #    legend=False,
            #    ci=None,
            # )
            sns.boxplot(
                x=groupby,
                data=df_melt.loc[lambda x: x["var"] == var],  # noqa: B023
                y="val",
                ax=ax,
                color="white",
                fliersize=0,
                **boxplot_properties,
            )

            ax.set_xlabel("")
            ax.tick_params(
                axis="x",
                # rotation=0,
                labelsize=9,
            )
            ax.legend().set_visible(False)
            ax.set_ylabel(ylabel)
            ax.set_title(var + "\n" + pvalue_template(pvalues[i]))
        else:
            ax.set_visible(False)
    fig.tight_layout()

    if show_legend is True:
        axes[n_panels - 1].legend().set_visible(True)
        axes[n_panels - 1].legend(bbox_to_anchor=(1.1, 1.05))

    if show:
        plt.show()

    if return_fig:
        return fig
