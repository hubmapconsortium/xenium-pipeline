#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
#import manhole
import matplotlib.pyplot as plt
import squidpy as sq
import spatialdata as sd

from common import Assay
from plot_utils import new_plot
import spatialdata_plot

from spatialdata.transformations import (
    Affine,
    MapAxis,
    Scale,
    Sequence,
    Translation,
    get_transformation,
    set_transformation,
)

def main(assay: Assay, h5ad_file: Path, sdata_zarr: Path):
    if assay in {Assay.XENIUM}:
        sdata = sd.read_zarr(sdata_zarr)
        adata = anndata.read(h5ad_file)
        sdata.tables['table'] = adata


        with new_plot():
            sdata.pl.render_images("morphology_focus").pl.render_shapes(
                "cell_circles",
            ).pl.show(
                title=f"leiden cluster over Morphology image",
                figsize=(10, 5),
            )
            plt.savefig("spatial_scatter.pdf", bbox_inches="tight")

        adata.obsm["spatial"] = adata.obsm["X_spatial"]

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.nhood_enrichment(adata, cluster_key="leiden")
            plt.savefig("neighborhood_enrichment.pdf", bbox_inches="tight")

        sq.gr.co_occurrence(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.co_occurrence(adata, cluster_key="leiden")
            plt.savefig("co_occurrence.pdf", bbox_inches="tight")

        sq.gr.centrality_scores(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.centrality_scores(adata, cluster_key="leiden")
            plt.savefig("centrality_scores.pdf", bbox_inches="tight")

        sq.gr.interaction_matrix(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.interaction_matrix(adata, cluster_key="leiden")
            plt.savefig("interaction_matrix.pdf", bbox_inches="tight")

        #        sq.gr.ripley(adata, cluster_key="leiden")

        #        with new_plot():
        #            sq.pl.ripley(adata, cluster_key="leiden")
        #            plt.savefig("ripley.pdf", bbox_inches="tight")

        output_file = Path("squidpy_annotated.h5ad")
        print("Saving output to", output_file.absolute())
        # Save normalized/etc. data
        adata.write_h5ad(output_file)


if __name__ == "__main__":
#    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("sdata_zarr", type=Path)

    args = p.parse_args()

    main(args.assay, args.h5ad_file, args.sdata_zarr)
