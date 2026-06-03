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


def write_scatter_tri_interpolant(
    prefix,
    node_xy,
    tri,
    IDs,
    ELEV=None,
    data=None,
    HOR_type="LINEAR",
    VER_type="LINEAR",
    return_filenames=False,
):
    prefix = Path(prefix)

    scatter_filename = prefix.with_name(prefix.name + "_scatter_tri.dat")
    data_filename = prefix.with_name(prefix.name + "_scatter_data.dat")

    node_xy = np.asarray(node_xy, dtype=float)
    tri = np.asarray(tri, dtype=int)
    IDs = np.asarray(IDs, dtype=int)

    if node_xy.ndim != 2 or node_xy.shape[1] != 2:
        raise ValueError("node_xy must have shape (Nnodes, 2)")

    if tri.ndim != 2 or tri.shape[1] != 3:
        raise ValueError("tri must have shape (Ntri, 3)")

    n_nodes = node_xy.shape[0]
    n_tri = tri.shape[0]

    if IDs.ndim == 1:
        IDs = IDs.reshape(-1, 1)

    if IDs.shape[0] != n_nodes:
        raise ValueError("IDs must have the same number of rows as node_xy")

    if ELEV is not None:
        ELEV = np.asarray(ELEV, dtype=float)

        if ELEV.ndim == 1:
            ELEV = ELEV.reshape(-1, 1)

        if ELEV.shape[0] != n_nodes:
            raise ValueError("ELEV must have the same number of rows as node_xy")

        if IDs.shape[1] != ELEV.shape[1] + 1:
            raise ValueError("IDs must have one more column than ELEV")
    else:
        if IDs.shape[1] != 1:
            raise ValueError("If ELEV is None, IDs must have only one column")

    with open(scatter_filename, "w") as f:
        f.write(f"{n_nodes:d} {n_tri:d}\n")
        f.write(f"{HOR_type} {VER_type}\n")

        for i in range(n_nodes):
            row = [
                f"{node_xy[i, 0]:.12g}",
                f"{node_xy[i, 1]:.12g}",
            ]

            if ELEV is None:
                row.append(f"{IDs[i, 0]:d}")
            else:
                for j in range(ELEV.shape[1]):
                    row.append(f"{IDs[i, j]:d}")
                    row.append(f"{ELEV[i, j]:.12g}")

                row.append(f"{IDs[i, -1]:d}")

            f.write(" ".join(row) + "\n")

        for t in tri:
            f.write(f"{int(t[0])} {int(t[1])} {int(t[2])}\n")

    if data is not None:
        data = np.asarray(data, dtype=float)
        np.savetxt(data_filename, data, fmt="%.12g")

    if return_filenames:
        return str(scatter_filename.name), str(data_filename.name)

def write_interpolation_master(prefix, regions):
    """
       Write an interpolation master file describing one or more interpolation regions.

       Parameters
       ----------
       prefix : str or Path
           Output file prefix. The master file is written as:

               <prefix>_master.dat

       regions : list of dict
           List of interpolation region definitions. Each dictionary must contain:

           Type : str
               Interpolation type. Supported values depend on the interpolation
               framework (e.g. "CONST", "SCATTER", "GRIDDED").

           values_file : str
               File containing the interpolation values.

           Optional keys
           -------------

           spatial_file : str
               Spatial interpolation file. Required for non-CONST interpolants.

           ll_pnt : array_like, shape (2,)
               Lower-left corner [x, y] of a rectangular region.

           ur_pnt : array_like, shape (2,)
               Upper-right corner [x, y] of a rectangular region.

           region : array_like, shape (nv, 2)
               Polygon defining the interpolation region. Each row contains:

                   [x, y]

               Vertices should be ordered around the polygon perimeter
               (clockwise or counter-clockwise). The polygon does not need
               to be explicitly closed; the first and last vertices will be
               connected automatically.

       Region specification
       --------------------
       Each region can be defined in one of two ways:

       1. Bounding box region
          Provide both:

              ll_pnt = [xmin, ymin]
              ur_pnt = [xmax, ymax]

          In this case nv = 2 is written to the master file and the two
          corner coordinates are printed.

       2. Polygon region
          Provide:

              region = np.array([
                  [x1, y1],
                  [x2, y2],
                  ...
                  [xnv, ynv]
              ])

          In this case nv is the number of polygon vertices and all
          vertices are written to the master file.

       Notes
       -----
       For CONST interpolants the spatial file is ignored and a dummy
       filename is written to the master file.
       """

    prefix = Path(prefix)
    master_filename = prefix.with_name(prefix.name + "_master.dat")

    def _is_empty(x):
        return x is None or np.asarray(x).size == 0

    with open(master_filename, "w") as f:

        f.write(f"{len(regions)}\n")

        for item in regions:

            interp_type = item["Type"]
            spatial_file = item.get("spatial_file", "")
            values_file = item["values_file"]

            ll_pnt = item.get("ll_pnt", None)
            ur_pnt = item.get("ur_pnt", None)
            region = item.get("region", None)

            has_bbox = not _is_empty(ll_pnt) and not _is_empty(ur_pnt)

            # determine nv
            if has_bbox:
                nv = 2
            else:
                region = np.asarray(region, dtype=float)

                if region.ndim != 2 or region.shape[1] != 2:
                    raise ValueError("region must have shape (nv,2)")

                nv = region.shape[0]

            # print header line
            if interp_type.upper() == "CONST":
                f.write(f"{nv} {interp_type} {"dummy_file.dat"} {values_file}\n")
            else:
                f.write(f"{nv} {interp_type} {spatial_file} {values_file}\n")

            # print bbox or polygon
            if has_bbox:

                ll_pnt = np.asarray(ll_pnt, dtype=float).ravel()
                ur_pnt = np.asarray(ur_pnt, dtype=float).ravel()

                f.write(f"{ll_pnt[0]:.12g} {ll_pnt[1]:.12g}\n")
                f.write(f"{ur_pnt[0]:.12g} {ur_pnt[1]:.12g}\n")

            else:

                for x, y in region:
                    f.write(f"{x:.12g} {y:.12g}\n")

def split_prop_elev_matrix(A, split_output=False):
    """
    Input columns:
    [HK1, ELEV12, HK2, ELEV23, ..., HKn]

    Parameters
    ----------
    A : ndarray
        Input matrix.
    split_output : bool, default=False
        If False:
            returns (id_elev_mat, hk_linear_mat)

        If True:
            returns (id_mat, elev_mat, hk_linear_mat)

    Returns
    -------
    If split_output=False:
        id_elev_mat : ndarray
            [ID1, ELEV12, ID2, ELEV23, ..., IDn]
        hk_linear_mat : ndarray

    If split_output=True:
        id_mat : ndarray
            IDs only, shape (nrows, nlayers)
        elev_mat : ndarray
            Elevations only, shape (nrows, nlayers-1)
        hk_linear_mat : ndarray
    """

    A = np.asarray(A)

    if A.ndim != 2:
        raise ValueError("A must be 2D")

    if A.shape[1] % 2 == 0:
        raise ValueError(
            "A must have columns [HK1, ELEV12, HK2, ELEV23, ..., HKn]"
        )

    n_rows = A.shape[0]

    hk_cols = np.arange(0, A.shape[1], 2)
    elev_cols = np.arange(1, A.shape[1], 2)

    # HK values flattened layer-by-layer
    hk_linear_mat = A[:, hk_cols].T.ravel()

    # IDs:
    # layer 1 -> 0..n_rows-1
    # layer 2 -> n_rows..2*n_rows-1
    # ...
    id_mat = np.arange(len(hk_cols) * n_rows).reshape(len(hk_cols), n_rows).T

    if split_output:
        elev_mat = A[:, elev_cols]
        return id_mat, elev_mat, hk_linear_mat

    id_elev_mat = A.copy()
    id_elev_mat[:, hk_cols] = id_mat

    return id_elev_mat, hk_linear_mat
