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
                f.write(f"{nv} {interp_type} {'dummy_file.dat'} {values_file}\n")
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

def write_const_interpolation(prefix, region, arr):
    """
    Write files for a spatially constant interpolation.

    The value is constant in space over the selected region, but `arr` may
    contain multiple time values.

    This writes:

        prefix.dat
        prefix_master.dat

    Parameters
    ----------
    prefix : str or Path
        Output file prefix.

    region : dict
        Dictionary defining exactly one region.

        Option 1: polygon region
            {
                "region": array-like of shape (N, 2)
            }

        Option 2: bounding box region
            {
                "ll_pnt": [xmin, ymin],
                "ur_pnt": [xmax, ymax]
            }

        The dictionary may contain both, but if `ll_pnt` and `ur_pnt`
        are provided, the bounding-box definition is used.

    arr : array-like
        Constant interpolation values. Can be scalar, 1D, or 2D.
        Written to `prefix.dat`.

    Returns
    -------
    tuple[str, str], optional
        `(master_filename, values_filename)` if return_filenames=True.
    """

    prefix = Path(prefix)

    values_filename = prefix.with_name(prefix.name + "_values.dat")

    arr = np.asarray(arr, dtype=float)
    np.savetxt(values_filename, np.atleast_2d(arr), fmt="%.12g")

    ll_pnt = region.get("ll_pnt", None)
    ur_pnt = region.get("ur_pnt", None)
    polygon = region.get("region", None)

    def _is_empty(x):
        return x is None or np.asarray(x).size == 0

    has_bbox = not _is_empty(ll_pnt) and not _is_empty(ur_pnt)
    has_polygon = not _is_empty(polygon)

    if not has_bbox and not has_polygon:
        raise ValueError(
            "region must define either ll_pnt/ur_pnt or a polygon using key 'region'"
        )

    if has_bbox:
        ll_pnt = np.asarray(ll_pnt, dtype=float).ravel()
        ur_pnt = np.asarray(ur_pnt, dtype=float).ravel()

        if ll_pnt.size != 2 or ur_pnt.size != 2:
            raise ValueError("ll_pnt and ur_pnt must each be [x, y]")

        master_region = {
            "region": None,
            "Type": "CONST",
            "spatial_file": "dummy_file.dat",
            "values_file": values_filename.name,
            "ll_pnt": ll_pnt,
            "ur_pnt": ur_pnt,
        }

    else:
        polygon = np.asarray(polygon, dtype=float)

        if polygon.ndim != 2 or polygon.shape[1] != 2:
            raise ValueError("region['region'] must have shape (N, 2)")

        if polygon.shape[0] <= 2:
            raise ValueError("polygon region must have more than 2 points")

        master_region = {
            "region": polygon,
            "Type": "CONST",
            "spatial_file": "dummy_file.dat",
            "values_file": values_filename.name,
            "ll_pnt": None,
            "ur_pnt": None,
        }

    write_interpolation_master( prefix,[master_region])


def write_npsat_mesh(filename, xy, ids, cw=False):
    """
    Write an NPSAT mesh file.

    Parameters
    ----------
    filename : str or Path
        Output file name.

    xy : ndarray
        Node coordinates with shape (n_xy, 2).

    ids : ndarray
        Element connectivity with shape (n_elem, 4).
        IDs are assumed to be zero-based.

    cw : bool, default=False
        If False, elements are written counter-clockwise.
        If True, elements are written clockwise.
    """

    xy = np.asarray(xy, dtype=float)
    ids = np.asarray(ids, dtype=int).copy()

    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError("xy must have shape (n_xy, 2).")

    if ids.ndim != 2 or ids.shape[1] != 4:
        raise ValueError("ids must have shape (n_elem, 4).")

    if ids.min() < 0 or ids.max() >= xy.shape[0]:
        raise ValueError("ids contain node indices outside the xy array.")

    # Coordinates of each element: shape (n_elem, 4, 2)
    elem_xy = xy[ids]

    # Signed polygon area using shoelace formula
    x = elem_xy[:, :, 0]
    y = elem_xy[:, :, 1]

    signed_area = 0.5 * np.sum(
        x * np.roll(y, -1, axis=1) -
        y * np.roll(x, -1, axis=1),
        axis=1
    )

    # signed_area > 0  => CCW
    # signed_area < 0  => CW
    is_ccw = signed_area > 0

    desired_ccw = not cw

    flip_mask = is_ccw != desired_ccw

    # Reverse node order where orientation is wrong
    ids[flip_mask, :] = ids[flip_mask, ::-1]

    n_xy = xy.shape[0]
    n_elem = ids.shape[0]

    filename = Path(filename)

    with open(filename, "w") as f:
        f.write(f"{n_xy} {n_elem}\n")
        np.savetxt(f, xy, fmt="%.2f")
        np.savetxt(f, ids, fmt="%d")

def read_npsat_mesh(filename):
    """
    Read an NPSAT mesh file.

    Parameters
    ----------
    filename : str or Path
        Input mesh file.

    Returns
    -------
    xy : ndarray
        Node coordinates with shape (n_xy, 2).

    ids : ndarray
        Element connectivity with shape (n_elem, 4).
        IDs are returned exactly as stored in the file.
    """

    with open(filename, "r") as f:

        # Header
        n_xy, n_elem = map(int, f.readline().split())

        # Coordinates
        xy = np.loadtxt(
            [next(f) for _ in range(n_xy)],
            dtype=float
        )

        # Connectivity
        ids = np.loadtxt(
            [next(f) for _ in range(n_elem)],
            dtype=int
        )

    xy = np.atleast_2d(xy)
    ids = np.atleast_2d(ids)

    if xy.shape != (n_xy, 2):
        raise ValueError(
            f"Expected coordinates shape ({n_xy}, 2), got {xy.shape}"
        )

    if ids.shape != (n_elem, 4):
        raise ValueError(
            f"Expected connectivity shape ({n_elem}, 4), got {ids.shape}"
        )

    return xy, ids


def _values_to_4d(Values):
    """
    Return Values as (ntime, nrow, ncol, nlay).
    """
    Values = np.asarray(Values)

    if Values.ndim == 2:
        # (nrow, ncol)
        return Values[np.newaxis, :, :, np.newaxis]

    if Values.ndim == 3:
        # (nrow, ncol, nlay)
        return Values[np.newaxis, :, :, :]

    if Values.ndim == 4:
        # (ntime, nrow, ncol, nlay)
        return Values

    raise ValueError(
        "Values must be 2D (nrow, ncol), "
        "3D (nrow, ncol, nlay), or "
        "4D (ntime, nrow, ncol, nlay)"
    )


def _ids_to_3d(IDs):
    """
    Return IDs as (nlay, nrow, ncol).
    """
    IDs = np.asarray(IDs)

    if IDs.ndim == 2:
        return IDs[np.newaxis, :, :]

    if IDs.ndim == 3:
        return IDs

    raise ValueError("IDs must be 2D (nrow, ncol) or 3D (nlay, nrow, ncol)")


def _elev_to_3d(Elev, nrow, ncol):
    """
    Return Elev as (nelev, nrow, ncol), or None.
    """
    elev_is_empty = Elev is None or np.asarray(Elev).size == 0

    if elev_is_empty:
        return None

    Elev = np.asarray(Elev)

    if Elev.ndim == 2:
        Elev_3d = Elev[np.newaxis, :, :]

    elif Elev.ndim == 3:
        # Accept either (nelev, nrow, ncol) or (nrow, ncol, nelev)
        if Elev.shape[1:] == (nrow, ncol):
            Elev_3d = Elev
        elif Elev.shape[:2] == (nrow, ncol):
            Elev_3d = np.moveaxis(Elev, 2, 0)
        else:
            raise ValueError(
                f"Elev shape {Elev.shape} is incompatible with "
                f"(nrow, ncol)=({nrow}, {ncol})"
            )

    else:
        raise ValueError("Elev must be None, empty, 2D, or 3D")

    if Elev_3d.shape[1:] != (nrow, ncol):
        raise ValueError(
            f"Elev grid shape {Elev_3d.shape[1:]} does not match "
            f"grid shape {(nrow, ncol)}"
        )

    return Elev_3d


def write_gridded_interpolant(prefix, IDs, Elev, grid, values,
                              return_filenames=False,
                              raster_order=True,):
    """
    Write gridded interpolant files:
      prefix_gridded_grid.dat
      prefix_gridded_data.dat

    IDs
        2D: (nrow, ncol)
        3D: (nlay, nrow, ncol)

    Elev
        None, empty, 2D, or 3D.
        3D may be either:
          (nelev, nrow, ncol)
        or:
          (nrow, ncol, nelev)

    values
        Usually 2D:
          (n_id, ntime)

        A 1D array is written as one column.
    """

    prefix = Path(prefix)

    grid_filename = prefix.with_name(prefix.name + "_gridded_grid.dat")
    data_filename = prefix.with_name(prefix.name + "_gridded_data.dat")

    IDs_3d = _ids_to_3d(IDs)

    nlay, nrow, ncol = IDs_3d.shape

    Elev_3d = _elev_to_3d(Elev, nrow, ncol)

    if raster_order:
        IDs_3d = IDs_3d[:, ::-1, :]
        if Elev_3d is not None:
            Elev_3d = Elev_3d[:, ::-1, :]

    values = np.asarray(values)

    if values.ndim == 1:
        values = values[:, np.newaxis]

    if values.ndim != 2:
        raise ValueError("values must be 1D or 2D, usually (n_id, ntime)")

    required_keys = ["xorig", "yorig", "dx", "dy"]
    for key in required_keys:
        if key not in grid:
            raise KeyError(f"grid dictionary is missing key: {key}")

    with open(grid_filename, "w") as f:
        f.write(f"{float(grid['xorig']):.6g} {ncol}\n")
        f.write(f"{float(grid['yorig']):.6g} {nrow}\n")
        f.write(f"{float(grid['dx']):.6g} {float(grid['dy']):.6g}\n")
        f.write(f"{nlay}\n")

        for ilay in range(nlay):
            np.savetxt(f, IDs_3d[ilay], fmt="%d")

            if Elev_3d is not None:
                elev_idx = min(ilay, Elev_3d.shape[0] - 1)
                np.savetxt(f, Elev_3d[elev_idx], fmt="%.12g")

    np.savetxt(data_filename, values, fmt="%.12g")

    if return_filenames:
        return str(grid_filename.name), str(data_filename.name)


def write_partitioned_gridded_interpolant(
    prefix,
    Values,
    Elev,
    xorig,
    yorig,
    dx,
    dy,
    NR,
    NC,
    return_filenames=False,
    zero_tol=0.0,
):
    """
    Write partitioned gridded interpolants, storing data only for nonzero cells.

    Values accepted shapes:
      (nrow, ncol)
      (nrow, ncol, nlay)
      (ntime, nrow, ncol, nlay)

    Internally treated as:
      (ntime, nrow, ncol, nlay)

    For each subregion:
      - inactive / zero cells receive ID = -1
      - active cells receive IDs 0, 1, ..., n_active-1
      - data file has only n_active rows
    """

    prefix = Path(prefix)
    master_filename = prefix.with_name(prefix.name + "_master.dat")

    Values_4d = _values_to_4d(Values)
    ntime, nrow, ncol, nlay = Values_4d.shape

    Elev_3d = _elev_to_3d(Elev, nrow, ncol)

    row_blocks = np.array_split(np.arange(nrow), NR)
    col_blocks = np.array_split(np.arange(ncol), NC)

    x_edges = xorig + np.arange(ncol + 1) * dx
    y_edges = yorig + np.arange(nrow + 1) * dy

    region_records = []
    master_lines = []

    for ir, rows in enumerate(row_blocks):
        for ic, cols in enumerate(col_blocks):

            r0, r1 = rows[0], rows[-1] + 1
            c0, c1 = cols[0], cols[-1] + 1

            sub_values = Values_4d[:, r0:r1, c0:c1, :]
            _, sub_nrow, sub_ncol, sub_nlay = sub_values.shape

            # Active if any timestep is nonzero for that row/col/layer
            active_rc_lay = np.any(
                np.isfinite(sub_values) & (np.abs(sub_values) > zero_tol),
                axis=0
            )  # shape: (sub_nrow, sub_ncol, sub_nlay)

            n_active = int(np.count_nonzero(active_rc_lay))

            # Skip completely zero subregions
            if n_active == 0:
                continue

            # IDs in row-col-layer layout first
            sub_ids_rc_lay = np.full(
                (sub_nrow, sub_ncol, sub_nlay),
                -1,
                dtype=int
            )

            sub_ids_rc_lay[active_rc_lay] = np.arange(n_active, dtype=int)

            # Convert to layer-row-col layout expected by write_gridded_interpolant
            sub_ids = np.moveaxis(sub_ids_rc_lay, 2, 0)

            # Data table: n_active x ntime
            # This ordering is consistent with active_rc_lay boolean indexing.
            sub_data = np.zeros((n_active, ntime), dtype=float)

            for itime in range(ntime):
                sub_data[:, itime] = sub_values[itime][active_rc_lay]


            if Elev_3d is None:
                sub_elev = None
            else:
                sub_elev = Elev_3d[:, r0:r1, c0:c1]

            sub_prefix = prefix.with_name(
                f"{prefix.name}_r{ir:03d}_c{ic:03d}"
            )

            sub_grid = {
                "xorig": xorig + c0 * dx,
                "yorig": yorig + (nrow - r1)*dy, #yorig + r0 * dy,
                "dx": dx,
                "dy": dy,
            }

            grid_file, data_file = write_gridded_interpolant(
                sub_prefix,
                sub_ids,
                sub_elev,
                sub_grid,
                sub_data,
                return_filenames=True,
                raster_order=True
            )

            xmin = min(x_edges[c0], x_edges[c1])
            xmax = max(x_edges[c0], x_edges[c1])
            ymin = yorig + (nrow - r1)*dy #min(y_edges[r0], y_edges[r1])
            ymax = yorig + (nrow - r0)*dy #max(y_edges[r0], y_edges[r1])

            master_lines.append(
                f"2 GRIDDED {grid_file} {data_file}\n"
                f"{xmin:.6f} {ymin:.6f}\n"
                f"{xmax:.6f} {ymax:.6f}\n"
            )

            region_records.append({
                "ir": ir,
                "ic": ic,
                "row_start": r0,
                "row_end": r1,
                "col_start": c0,
                "col_end": c1,
                "ntime": ntime,
                "nlay": sub_nlay,
                "n_active": n_active,
                "grid_file": grid_file,
                "data_file": data_file,
                "xmin": xmin,
                "ymin": ymin,
                "xmax": xmax,
                "ymax": ymax,
            })

    with open(master_filename, "w") as fmaster:
        fmaster.write(f"{len(master_lines)}\n")
        for line in master_lines:
            fmaster.write(line)

    if return_filenames:
        return str(master_filename.name), region_records


def calc_relative_layer_positions(elev, top_elev):
    """
    Calculate relative vertical positions of intermediate layer elevations
    between top_elev and the bottom elevation elev[:, -1].

    Returns
    -------
    vert_rel : np.ndarray
        Array with shape (n_nodes, n_elev - 2).

        Relative position is:
            0 at top_elev
            1 at bottom elevation

        Example:
            top = 100, bottom = 0, intermediates = 75, 50, 25
            returns [0.25, 0.50, 0.75]

    Notes
    -----
    If one or more intermediate elevations are above top_elev, those
    positions are redistributed evenly between top_elev and the first
    intermediate elevation below top_elev.
    """

    elev = np.asarray(elev, dtype=float)
    top_elev = np.asarray(top_elev, dtype=float)

    if elev.ndim != 2:
        raise ValueError("elev must be a 2D array with shape (n_nodes, n_elev)")

    n_nodes, n_elev = elev.shape

    if n_elev < 3:
        raise ValueError("elev must have at least 3 columns: top, intermediate, bottom")

    if top_elev.shape[0] != n_nodes:
        raise ValueError("top_elev must have length n_nodes")

    bot_elev = elev[:, -1]
    mid_elev = elev[:, 1:-1]

    vert_rel = np.full((n_nodes, n_elev - 2), np.nan, dtype=float)

    for i in range(n_nodes):

        top = top_elev[i]
        bot = bot_elev[i]
        mids = mid_elev[i, :].copy()

        den = top - bot

        if not np.isfinite(den) or den == 0:
            continue

        # Standard relative positions
        rel = (top - mids) / den

        # Find intermediates that are above top_elev
        above = mids > top

        if np.any(above):

            # First layer that is below or equal to top_elev
            below_idx = np.where(mids <= top)[0]

            if below_idx.size == 0:
                # All intermediate layers are above top_elev.
                # Distribute them evenly between top and bottom.
                rel = np.linspace(
                    0,
                    1,
                    n_elev
                )[1:-1]

            else:
                k = below_idx[0]

                # Relative position of first valid below-top layer
                rel_k = (top - mids[k]) / den

                # Distribute layers 0..k-1 evenly between 0 and rel_k
                rel[:k] = np.linspace(
                    0,
                    rel_k,
                    k + 2
                )[1:-1]

                # Keep normal relative positions from k onward
                rel[k:] = (top - mids[k:]) / den

        vert_rel[i, :] = rel

    return vert_rel
