import h5py
import numpy as np

def hdf5_structure(filename):
    """
    Return nested dictionary describing groups and datasets in an HDF5 file.
    Dataset data are not loaded.
    """

    def explore(obj):
        out = {}

        for name, item in obj.items():
            if isinstance(item, h5py.Group):
                out[name] = {
                    "type": "Group",
                    "children": explore(item)
                }

            elif isinstance(item, h5py.Dataset):
                out[name] = {
                    "type": "Dataset",
                    "shape": item.shape,
                    "dtype": str(item.dtype),
                    "size": item.size
                }

        return out

    with h5py.File(filename, "r") as f:
        return explore(f)

def hdf5_structure_detailed(filename):
    """
    Return nested dictionary describing groups, datasets, and attached
    HDF5 attributes in an HDF5 file.

    Dataset data are not loaded.
    Attribute values are loaded because they are usually small.
    """

    def convert_attr_value(value):
        """Convert HDF5/numpy attribute values to readable Python objects."""
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")

        if isinstance(value, np.ndarray):
            if value.dtype.kind == "S":  # byte strings
                return [
                    v.decode("utf-8", errors="replace")
                    for v in value
                ]
            return value.tolist()

        if isinstance(value, np.generic):
            return value.item()

        return value

    def read_attrs(obj):
        """Read attached attributes of an HDF5 object."""
        return {
            name: convert_attr_value(value)
            for name, value in obj.attrs.items()
        }

    def explore(obj):
        out = {}

        for name, item in obj.items():

            if isinstance(item, h5py.Group):
                out[name] = {
                    "type": "Group",
                    "attributes": read_attrs(item),
                    "children": explore(item),
                }

            elif isinstance(item, h5py.Dataset):
                out[name] = {
                    "type": "Dataset",
                    "shape": item.shape,
                    "dtype": str(item.dtype),
                    "size": item.size,
                    "attributes": read_attrs(item),
                }

        return out

    with h5py.File(filename, "r") as f:
        return {
            "type": "File",
            "attributes": read_attrs(f),
            "children": explore(f),
        }