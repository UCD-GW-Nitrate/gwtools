import h5py

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