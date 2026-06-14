import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, Point, LineString

def mesh_to_gdf(nodes, elements, crs):
    node_lookup = nodes.set_index("ID")[["X", "Y"]]
    geometries = []
    ies = []
    irges = []

    for _, row in elements.iterrows():
        nids = [row["ND1"], row["ND2"], row["ND3"]]
        if row["ND4"] != 0:
            nids.append(row["ND4"])

        # Convert node IDs to coordinates
        coords = [(node_lookup.loc[nid, "X"], node_lookup.loc[nid, "Y"]) for nid in nids]

        # Create polygon
        poly = Polygon(coords)

        geometries.append(poly)
        ies.append(row["IE"])
        irges.append(row["IRGE"])

    gdf = gpd.GeoDataFrame(
        {"IE": ies, "IRGE": irges},
        geometry=geometries,
        crs=crs
    )

    return gdf

def get_date_df(start_date, end_date,freq="MS"):
    c2vs_dates = pd.date_range(start=start_date, end=end_date, freq=freq)
    c2vs_dates_df = pd.DataFrame({
        'date': c2vs_dates,
        'year': c2vs_dates.year,
        'month': c2vs_dates.month,
        'days_in_month': c2vs_dates.days_in_month
    })
    return c2vs_dates_df

def stratigraphy_to_elevations(
    strat_df,
    conversion_factor=0.3048,
    id_col="ID",
    gse_col="GSE",
    remove_zero_layers=True,
    verbose=True,
):
    """
        Convert IWFM-style stratigraphy thickness columns to elevation surfaces.

        Parameters
        ----------
        strat_df : DataFrame
            Stratigraphy table with fields:
            ID, GSE, A1, L1, A2, L2, ...
        conversion_factor : float, default=0.3048
            Unit conversion factor. Default converts feet to meters.
            Use 1.0 if the input is already in desired units.
        id_col : str, default="ID"
            Node/point ID field.
        gse_col : str, default="GSE"
            Ground surface elevation field.
        remove_zero_layers : bool, default=True
            If True, remove elevation surfaces where the corresponding
            layer thickness is zero for all rows.
        verbose : bool, default=True
            If True, print reminder messages.

        Returns
        -------
        ids : ndarray
            Original ID values.
        elev : ndarray
            Elevation matrix with shape (n_points, n_surfaces).
            First column is GSE.
        elev_cols : list[str]
            Names of elevation surfaces.
        """

    ids = strat_df[id_col].to_numpy()

    # Detect A/L columns
    a_cols = sorted(
        [c for c in strat_df.columns if c.startswith("A")],
        key=lambda x: int(x[1:])
    )
    l_cols = sorted(
        [c for c in strat_df.columns if c.startswith("L")],
        key=lambda x: int(x[1:])
    )

    thickness_cols = []
    for a, l in zip(a_cols, l_cols):
        thickness_cols.extend([a, l])

    # Thickness matrix
    thickness = strat_df[thickness_cols].to_numpy(dtype=float)

    # Convert units
    gse = strat_df[gse_col].to_numpy(dtype=float) * conversion_factor
    thickness = thickness * conversion_factor

    # Build elevation surfaces: GSE, GSE-A1, GSE-A1-L1, ...
    elev = np.column_stack([
        gse,
        gse[:, None] - np.cumsum(thickness, axis=1)
    ])

    elev_cols = [gse_col] + [f"BOT_{c}" for c in thickness_cols]

    # Remove surfaces associated with layers that are zero everywhere
    if remove_zero_layers:
        keep = [True]  # always keep GSE

        for j in range(thickness.shape[1]):
            keep.append(not np.all(thickness[:, j] == 0))

        keep = np.array(keep, dtype=bool)

        elev = elev[:, keep]
        elev_cols = [c for c, k in zip(elev_cols, keep) if k]

    # Replace zero-thickness elevations by average of above and below
    corrected_count = 0

    for j in range(1, elev.shape[1] - 1):
        zero_mask = np.isclose(elev[:, j], elev[:, j - 1])

        if np.any(zero_mask):
            elev[zero_mask, j] = 0.5 * (
                    elev[zero_mask, j - 1] +
                    elev[zero_mask, j + 1]
            )
            corrected_count += int(np.sum(zero_mask))


    if verbose and corrected_count > 0:
        print(
            f"Reminder: {corrected_count} zero-thickness elevation values "
            "were replaced by the mean elevation of the surfaces above and below."
        )

    return ids, elev, elev_cols


def embed_hk_anomaly_to_nodes(
    HVK, # Horizontal vertical K properties
    HKanomaly,
    elem_df,
    id_col="ID",
    hk_col="PKH",
    elem_id_col="IE",
    elem_node_cols=("ND1", "ND2", "ND3", "ND4"),
    anomaly_elem_col=1,
    anomaly_first_layer_col=2,
):
    """
    Embed element-based HK anomaly values into node-based HKV dataframes.

    Parameters
    ----------
    HVK : list of DataFrame
        One dataframe per layer. Each dataframe must contain node IDs and
        hydraulic conductivity values, e.g. ID and PKH.
    HKanomaly : ndarray
        Array where:
        column 0 = unused ID,
        column 1 = element ID,
        columns 2 onward = anomaly values per layer.
    elem_df : DataFrame
        Element connectivity table with fields IE, ND1, ND2, ND3, ND4.
        If ND4 == 0, the element is treated as a triangle.

    Returns
    -------
    HKV_new : list of DataFrame
        Updated copy of HKV with anomaly values assigned to nodes touching
        anomalous elements.
    """

    nlay = len(HVK)

    if HKanomaly.shape[1] - anomaly_first_layer_col != nlay:
        raise ValueError(
            "Number of anomaly layer columns does not match len(HKV). "
            f"Got {HKanomaly.shape[1] - anomaly_first_layer_col} anomaly layers "
            f"and {nlay} HKV layers."
        )

    elem_node_cols = [c for c in elem_node_cols if c in elem_df.columns]

    # Element ID -> node IDs
    elem_conn = elem_df.set_index(elem_id_col)[elem_node_cols]

    # Work on copies so original HKV is not modified
    HKV_new = [df.copy() for df in HVK]

    # For fast node lookup in each layer dataframe
    node_index_maps = [
        pd.Series(df.index.to_numpy(), index=df[id_col]).to_dict()
        for df in HKV_new
    ]

    for row in HKanomaly:
        ie = row[anomaly_elem_col]

        if ie not in elem_conn.index:
            continue

        nodes = elem_conn.loc[ie].to_numpy()
        nodes = nodes[(~pd.isna(nodes)) & (nodes != 0)]

        for ilay in range(nlay):
            anomaly_value = row[anomaly_first_layer_col + ilay]

            idx_map = node_index_maps[ilay]

            for node_id in nodes:
                if node_id in idx_map:
                    HKV_new[ilay].at[idx_map[node_id], hk_col] = anomaly_value

    return HKV_new


def hvk_to_arrays(HVK, props):
    """
    Convert a layer-wise HVK list into node x layer numpy arrays.

    Parameters
    ----------
    HVK : list of DataFrame
        One dataframe per layer.
    props : str or list[str]
        Property name(s) to extract.

    Returns
    -------
    ndarray or tuple(ndarray)
        If props is a string:
            arr : (nnodes, nlayers)

        If props is a list:
            tuple of arrays in the same order as props.
    """

    single_prop = isinstance(props, str)

    if single_prop:
        props = [props]

    nlay = len(HVK)

    # sanity checks
    for p in props:
        if p not in HVK[0].columns:
            raise KeyError(f"Property '{p}' not found in HVK.")

    nnodes = len(HVK[0])

    out = {}

    for prop in props:
        arr = np.empty((nnodes, nlay), dtype=float)

        for ilay, df in enumerate(HVK):
            arr[:, ilay] = df[prop].to_numpy()

        out[prop] = arr

    if single_prop:
        return out[props[0]]

    return tuple(out[p] for p in props)


def watersheds_to_gdfs(watersheds, nodes_df, crs):
    """
    Convert IWFM small watershed connectivity to GeoDataFrames.

    Parameters
    ----------
    watersheds : list[dict]
        Output from read_iwfm_smallwatersheds().

    nodes_df : pandas.DataFrame
        DataFrame with fields:
        ID, X, Y

    crs : str, dict, pyproj.CRS, or int
        CRS passed directly to GeoDataFrame.

    Returns
    -------
    line_gdf : GeoDataFrame
    point_gdf : GeoDataFrame
    """

    required_cols = {"ID", "X", "Y"}
    missing = required_cols - set(nodes_df.columns)
    if missing:
        raise ValueError(f"nodes_df is missing required columns: {missing}")

    node_lookup = nodes_df.set_index("ID")[["X", "Y"]].to_dict("index")

    line_records = []
    point_records = []

    for ws in watersheds:
        ws_id = ws["ID"]
        iwb_list = ws["IWB"]
        qmaxwb_list = ws["QMAXWB"]

        if len(iwb_list) != len(qmaxwb_list):
            raise ValueError(
                f"Watershed ID {ws_id}: IWB and QMAXWB have different lengths"
            )

        base_attrs = {
            key: value
            for key, value in ws.items()
            if key not in {"IWB", "QMAXWB"}
        }

        coords = []

        for iwb, qmaxwb in zip(iwb_list, qmaxwb_list):
            if iwb not in node_lookup:
                raise KeyError(
                    f"Watershed ID {ws_id}: IWB node {iwb} not found in nodes_df"
                )

            x = node_lookup[iwb]["X"]
            y = node_lookup[iwb]["Y"]

            coords.append((x, y))

            point_records.append({
                **base_attrs,
                "IWB": iwb,
                "QMAXWB": qmaxwb,
                "geometry": Point(x, y),
            })

        # Skip watersheds with only one IWB
        if len(coords) > 1:
            line_records.append({
                **base_attrs,
                "IWB": iwb_list,
                "geometry": LineString(coords),
            })

    line_gdf = gpd.GeoDataFrame(line_records, geometry="geometry",crs=crs)
    point_gdf = gpd.GeoDataFrame(point_records, geometry="geometry",crs=crs)

    return line_gdf, point_gdf

