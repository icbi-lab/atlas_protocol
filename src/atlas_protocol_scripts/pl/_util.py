"""Helper functions for creating plots."""

import os
from pathlib import Path
from typing import Literal

import altair as alt
import matplotlib as mpl
import seaborn


def reshape_clustermap(cmap: seaborn.matrix.ClusterGrid, cell_width: float = 0.02, cell_height: float = 0.02):
    """Resizes the components of a seaborn clustermap object to match the specified cell widths and heights.

    Parameters
    ----------
    cmap
        The seaborn clustermap object to reshape.
    cell_width
        The width of each cell in the heatmap. Default is 0.02.
    cell_height
        The height of each cell in the heatmap. Default is 0.02.
    """
    # Get the number of rows and columns in the heatmap data
    ny, nx = cmap.data2d.shape

    # Calculate the new width and height of the heatmap based on the number of rows and columns and the cell sizes
    hmap_width = nx * cell_width
    hmap_height = ny * cell_height

    # Set the position of the heatmap axes to the new width and height
    hmap_orig_pos = cmap.ax_heatmap.get_position()
    cmap.ax_heatmap.set_position([hmap_orig_pos.x0, hmap_orig_pos.y0, hmap_width, hmap_height])

    # Set the position of the column dendrogram axes to the new width and height plus the height of the heatmap axes
    top_dg_pos = cmap.ax_col_dendrogram.get_position()
    cmap.ax_col_dendrogram.set_position(
        [hmap_orig_pos.x0, hmap_orig_pos.y0 + hmap_height, hmap_width, top_dg_pos.height]
    )

    # Set the position of the row dendrogram axes to the new width and height plus the width of the heatmap axes
    left_dg_pos = cmap.ax_row_dendrogram.get_position()
    cmap.ax_row_dendrogram.set_position([left_dg_pos.x0, left_dg_pos.y0, left_dg_pos.width, hmap_height])

    # If the heatmap has a colorbar, set its position to be below the heatmap axes
    if cmap.ax_cbar:
        cbar_pos = cmap.ax_cbar.get_position()
        hmap_pos = cmap.ax_heatmap.get_position()
        cmap.ax_cbar.set_position([cbar_pos.x0, hmap_pos.y1, cbar_pos.width, cbar_pos.height])


def save_fig_mfmt(
    fig,
    res_dir: str,
    filename: str,
    fmt: Literal["all", "pdf", "png", "svg"] = "all",
    **kwargs,
) -> None:
    """Save a matplotlib or altair figure in the specified format(s) to the given directory with the given filename.

    Parameters
    ----------
    fig
        The matplotlib or altair figure to save.
    res_dir
        The directory in which to save the figure.
    filename
        The name to use for the figure file.
    fmt
        The format in which to save the figure.
        Default is "all", which saves the figure in all supported formats.
        Supported formats are "pdf", "svg", "png".
    **kwargs
        Additional keyword arguments to be passed to the savefig method of matplotlib or the save method of altair.

    Returns
    -------
    None

    Raises
    ------
    AssertionError: If the specified plot_provider is not supported or if the specified fmt is not supported.
    """
    supported_formats = ["pdf", "svg", "png", "all"]

    # Determine the plot provider
    plot_provider = "unsupported"

    if type(fig) in [mpl.figure.Figure, seaborn.matrix.ClusterGrid]:
        plot_provider = "mpl"
    if type(fig) is alt.vegalite.v4.api.LayerChart:
        plot_provider = "altair"

    assert plot_provider in ["mpl", "altair"], "Plot provider not supported: " + plot_provider
    assert fmt in supported_formats, (
        "Error: format " + fmt + " not supported formats [" + ",".join(supported_formats) + "]"
    )

    # Remove the "all" option from the supported_formats list, since it is only used for determining if all formats should be saved
    supported_formats.remove("all")

    # If fmt is "all", save the figure in all supported formats
    if fmt == "all":
        for f in supported_formats:
            # Construct and create the full path for the figure file in the current format
            d = Path(res_dir, f)
            os.makedirs(d, mode=0o750, exist_ok=True)

            # Construct the full path and filename for the figure file in the current format
            fn = Path(d, f"{filename}.{f}")

            # Save the figure in the current format using the specified plot_provider library
            if plot_provider == "mpl":
                fig.savefig(fn, **kwargs)
            else:
                fig.save(fn, fmt=f, mode="vega-lite", method="node")

    # If fmt is not "all", save the figure in the specified format
    else:
        # Construct and create the full path for the figure file in the specified format
        d = Path(res_dir, fmt)
        os.makedirs(d, mode=0o750, exist_ok=True)

        # Construct the full path and filename for the figure file in the specified format
        fn = Path(d, f"{filename}.{fmt}")

        # Save the figure in the specified format using the specified plot_provider library
        if plot_provider == "mpl":
            fig.savefig(fn, **kwargs)
        else:
            fig.save(fn, fmt=fmt, mode="vega-lite", method="node")
