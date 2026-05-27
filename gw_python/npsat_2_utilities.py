from pathlib import Path
import numpy as np

def write_scatter_interpolant(prefix, node_xy, tri_ids, data,HOR_type,VER_type, return_filenames=False,):
    """
        Write scattered interpolant geometry/connectivity and values files.

        Files written
        -------------
        {prefix}_scatter.dat
            Nnodes Ntri
            HOR_type VER_type
            node_xy
            tri_ids

        {prefix}_values.dat
            data

        Returns
        -------
        Optional tuple:
            spatial_filename, data_filename
        """

    prefix = Path(prefix)

    spatial_filename = prefix.with_name(prefix.name + "_scatter.dat")
    data_filename = prefix.with_name(prefix.name + "_values.dat")

    node_xy = np.asarray(node_xy, dtype=float)
    tri_ids = np.asarray(tri_ids, dtype=int)
    data = np.asarray(data)

    if node_xy.ndim != 2 or node_xy.shape[1] != 2:
        raise ValueError("node_xy must have shape (Nnodes, 2)")

    if tri_ids.ndim != 2 or tri_ids.shape[1] != 3:
        raise ValueError("tri_ids must have shape (Ntri, 3)")

    if data.ndim not in (1, 2):
        raise ValueError("data must be either 1D or 2D")

    Nnodes = node_xy.shape[0]
    Ntri = tri_ids.shape[0]

    with open(spatial_filename, "w") as f:
        f.write(f"{Nnodes} {Ntri}\n")
        f.write(f"{HOR_type} {VER_type}\n")

        np.savetxt(f, node_xy, fmt="%.10g")

        np.savetxt(f, tri_ids, fmt="%d")

    with open(data_filename, "w") as f:
        if data.ndim == 1:
            np.savetxt(f, data.reshape(-1, 1), fmt="%.10g")
        else:
            np.savetxt(f, data, fmt="%.10g")

    if return_filenames:
        return str(spatial_filename), str(data_filename)
