#!/usr/bin/env python3

from common import Assay

import spatialdata as sd
#from spatialdata_io import xenium, cosmx, cosmx_proteomics
from spatialdata_io import xenium, cosmx, visium_hd, cosmx_proteomics
import anndata
import annsel as an

import re
import os
import shutil
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from os import walk, fspath
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
from contextlib import contextmanager


debug_out_dir = Path("crop-debug")
schema_url_pattern = re.compile(r"\{(.+)\}OME")
ome_tiff_pattern = re.compile(r"(?P<basename>.*)\.ome\.tif(f?)$")
nanostring_counts_file_pattern = re.compile(r"(?P<basename>.*)_exprMat_file\.csv\.gz$")
nanostring_meta_file_pattern = re.compile(r"(?P<basename>.*)_metadata_file\.csv\.gz$")
nanostring_fov_file_pattern = re.compile(r"(?P<basename>.*)_fov_positions_file.csv$")

XENIUM_ZARR_PATH = "Xenium.zarr"


def rearrange_data(data_directory: Path, proteomics:bool=False)->Path:
    mask_dir = 'ProteinMasks' if proteomics else 'CellLabels'
    image_dir = 'ProteinImages' if proteomics else 'CellComposite'
    data_directory = Path(data_directory) / 'raw/'
    directory = Path('cosmx')
    os.makedirs(directory)
    os.makedirs(directory / mask_dir)
    os.makedirs(directory / image_dir)
    for f in data_directory.iterdir():
        if f.is_file():
            shutil.copy(f, directory / f.name)
    for d in (data_directory / 'images').iterdir():
        if not d.is_dir():
            continue
        for f in d.iterdir():
            if 'CellLabels' in f.name:
                shutil.copy(f, directory / mask_dir)
            elif '.tif' in f.name or '.TIF' in f.name:
                shutil.copy(f, directory / f'{image_dir}/_{f.stem.replace('OV', '')}.tif')
    return directory

@contextmanager
def new_plot():
    """
    When used in a `with` block, clears matplotlib internal state
    after plotting and saving things. Probably not necessary to be this
    thorough in clearing everything, but extra calls to `plt.clf()` and
    `plf.close()` don't *hurt*

    Intended usage:
        ```
        with new_plot():
            do_matplotlib_things()

            plt.savefig(path)
            # or
            fig.savefig(path)
        ```
    """
    plt.clf()
    try:
        yield
    finally:
        plt.clf()
        plt.close()
def find_geojson(directory: Path) -> Optional[Path]:
    geojson_files = list(directory.glob("**/*.geojson"))

    if len(geojson_files) > 1:
        raise ValueError(f"Found multiple GeoJSON files in {directory}")
    elif len(geojson_files) == 1:
        return geojson_files[0]
    else:
        return None

def crop_sdata(sdata, geojson_path):
    with open(geojson_path) as f:
        crop_geometry = shapely.from_geojson(f.read())
        if not isinstance(crop_geometry, shapely.GeometryCollection):
            crop_geometry = shapely.geometrycollections([crop_geometry])

        closed_geometry = shapely.GeometryCollection(
            [shapely.Polygon(poly.exterior.coords) for poly in crop_geometry.geoms]
        )

    sdata_crop = sd.polygon_query(sdata, shapely.MultiPolygon(list(closed_geometry.geoms)), 'global', filter_table=True, clip=False)


    #Workaround issue writing out points
    sdata.points['transcripts'].compute()
    sdata_crop.points['transcripts'] = sd.models.PointsModel.parse(sdata_crop.points['transcripts'].compute())

#    sdata = sdata.query.polygon_query(shapely.MultiPolygon(list(closed_geometry.geoms)), 'global', filter_table=True, clip=False)
    return sdata_crop

def find_files(directory: Path, pattern) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def get_schema_url(ome_xml_root_node: ET.Element) -> str:
    if m := schema_url_pattern.match(ome_xml_root_node.tag):
        return m.group(1)
    raise ValueError(f"Couldn't extract schema URL from tag name {ome_xml_root_node.tag}")

def find_ome_tiffs(input_dir: Path) -> Iterable[Path]:
    """
    Yields 2-tuples:
     [0] full Path to source file
     [1] output file Path (source file relative to input_dir)
    """
    for dirpath_str, _, filenames in walk(input_dir):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            if ome_tiff_pattern.match(filename):
                src_filepath = dirpath / filename
                yield src_filepath

def subset_spatialdata(sdata: sd.SpatialData):
    fovs = [k.split('_')[0] for k in sdata.images]
    table_name = list(sdata.tables)[0]
    sdata_subset = sd.filter_by_table_query(sdata, table_name=table_name, obs_expr=an.col("fov").is_in(fovs))
    sdata_subset.images = sdata.images
    return sdata_subset


def main(assay: Assay, data_directory: Path):
    if assay == assay.XENIUM:
        sdata = xenium(data_directory / "lab_processed/xenium_bundle")

        maybe_geojson = find_geojson(data_directory)
        if maybe_geojson:
            sdata = crop_sdata(sdata, maybe_geojson)
        sdata.write(XENIUM_ZARR_PATH)
        sdata = sd.read_zarr(XENIUM_ZARR_PATH)
        adata = sdata.tables["table"]

    elif assay == assay.COSMX:

        directory = rearrange_data(data_directory)
        sdata = cosmx(directory)
        sdata = subset_spatialdata(sdata)
        sdata.write('CosMx.zarr')
        adata = sdata.tables["table"]

    elif assay == assay.COSMX_PROTEOMICS:
        data_directory = rearrange_data(data_directory, proteomics=True)
        sdata = cosmx_proteomics(data_directory)
        sdata.write('CosMx.zarr')
        adata = sdata.tables["table"]

    elif assay == assay.VISIUM_HD:
        directory = data_directory / "lab_processed/spaceranger/"
        sdata = visium_hd(directory)

        maybe_geojson = find_geojson(data_directory)
        if maybe_geojson:
            sdata = crop_sdata(sdata, maybe_geojson)
        sdata.write('Visium_HD.zarr')
        sdata = sd.read_zarr('Visium_HD.zarr')
        adata = sdata.tables["square_008um"]

    adata.obsm["X_spatial"] = adata.obsm["spatial"]

    adata.write("expr.h5ad")


if __name__ == "__main__":
#    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("data_directory", type=Path)
    args = p.parse_args()

    main(args.assay, args.data_directory)
