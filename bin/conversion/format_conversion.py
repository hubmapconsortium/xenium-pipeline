#!/usr/bin/env python3

from common import Assay

import spatialdata as sd
from spatialdata_io import xenium
import anndata

import re
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from os import walk, fspath
from pathlib import Path
from typing import Iterable, List, Tuple

import aicsimageio
import manhole
import tifffile as tf
from pint import Quantity, UnitRegistry

schema_url_pattern = re.compile(r"\{(.+)\}OME")
ome_tiff_pattern = re.compile(r"(?P<basename>.*)\.ome\.tif(f?)$")
nanostring_counts_file_pattern = re.compile(r"(?P<basename>.*)_exprMat_file\.csv\.gz$")
nanostring_meta_file_pattern = re.compile(r"(?P<basename>.*)_metadata_file\.csv\.gz$")
nanostring_fov_file_pattern = re.compile(r"(?P<basename>.*)_fov_positions_file.csv$")

XENIUM_ZARR_PATH = "Xenium.zarr"

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


def physical_dimension_func(img: aicsimageio.AICSImage) -> Tuple[List[float], List[str]]:
    """
    Returns lists of physical dimensions of pixels and corresponding units
    read from OME-XML metadata of input image
    """

    # aicsimageio parses the OME-XML metadata when loading an image,
    # and uses that metadata to populate various data structures in
    # the AICSImage object. The AICSImage.metadata.to_xml() function
    # constructs a new OME-XML string from that metadata, so anything
    # ignored by aicsimageio won't be present in that XML document.
    # Unfortunately, current aicsimageio ignores physical size units,
    # so we have to parse the original XML ourselves:
    root = ET.fromstring(img.xarray_dask_data.unprocessed[270])
    schema_url = get_schema_url(root)
    pixel_node_attrib = root.findall(f".//{{{schema_url}}}Pixels")[0].attrib

    values = []
    units = []
    for dimension in "XY":
        unit = pixel_node_attrib[f"PhysicalSize{dimension}Unit"]
        value = float(pixel_node_attrib[f"PhysicalSize{dimension}"])
        values.append(value)
        units.append(unit)

    return values, units

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

def main(assay: Assay, data_directory: Path):
    if assay == assay.XENIUM:
        sdata = xenium(data_directory)
        sdata.write(XENIUM_ZARR_PATH)
        sdata = sd.read_zarr(XENIUM_ZARR_PATH)
        adata = sdata.tables["table"]

        tiff_file = list(find_ome_tiffs(input_dir=data_directory))[0]
        img = tf.imread(fspath(tiff_file))
        library_id = list(adata.uns["spatial"].keys())[0]
        adata.uns["spatial"][library_id]["images"] = {"hires": img}
        adata.uns["spatial"][library_id]["scalefactors"] = {
            "tissue_hires_scalef": 1.0,
        }
        img = aicsimageio.AICSImage(tiff_file)
        values, units = physical_dimension_func(img)
        ureg = UnitRegistry()
        Q_ = ureg.Quantity

    elif assay == assay.COSMX:
        counts_file = find_files(data_directory, nanostring_counts_file_pattern)
        metadata_file = find_files(data_directory, nanostring_meta_file_pattern)
        fov_file = find_files(data_directory, nanostring_fov_file_pattern)
        adata = anndata.read.nanostring(path=data_directory, counts_file=counts_file, metadata_file=metadata_file, \
                                        fov_file=fov_file)

    adata.obsm["X_spatial"] = adata.obsm["spatial"]

    adata.write("expr.h5ad")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("data_directory", type=Path)
    args = p.parse_args()

    main(args.assay, args.data_directory)
