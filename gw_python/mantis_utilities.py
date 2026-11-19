import os
import subprocess
import numpy as np
import h5py
from dataclasses import dataclass
from scipy.stats import lognorm
from tqdm import tqdm

def mantis_inputs():
    """Return client settings and selectable Scenario.h inputs.

    Optional values of None leave the server defaults in effect. Scenario/map
    names must exist in the server configuration. Existing example selections
    are retained for compatibility.
    """

    scenario = {
        # -----------------
        # Mantis client options
        # -----------------
        'client' : None,
        'infile' : 'incomingMSG.dat',
        'outfile' : 'testClientResults.dat',
        'descr' : 'This is a description of the simulation input. It will be ignored',

        # -----------------
        # GUI Scenario Options
        # -----------------
        'modelArea' : 'CentralValley',
        'bMap' : 'Townships',
        'Regions' : ['39M02S10E', '39M02S09E'],
        'flowScen' : 'CVHM2_MAR24_Padj',
        'wellType' : 'VI',
        'unsatScen' : 'CVHM2_MAR24',
        'unsatWC' : 0.01,
        'endSimYear' : 2150,

        # -----------------
        # Loading options
        # -----------------
        'loadScen' : 'NONE',
        'startRed' : 0,
        'endRed' : 0,
        'Crops' : None,

        # -----------------
        # options that take default values
        # -----------------
        'minRch' : None, # 10;
        'rchMap' : None,
        'maxConc' : None, #250;
        'maxAge' : None, #400;
        'constRed' : None, #1.0;
        # Removal applies where the map value exceeds the threshold.
        # Set all three options to enable; denitMap=None disables denitrification.
        'denitMap' : None, # Configured map name, not a filename.
        'denitThreshold' : None, # 0..1
        'denitFractionRemoved' : None, # 0..1 (0.30 removes 30%).
        'LoadTransitionName' : None, #'NONE';
        'LoadTransitionStart' : None, #2005;
        'LoadTransitionEnd' : None, #2015;
        'LoadTransitionSubScen' : None, # Named historic RASTER loading subset.
        'startSimYear' : None, #1945;
        'por' : None, #10;
        'urfType' : None, #'LGNRM' or 'ADE';
        'ADELR' : None, # [lambda, retardation]
        'initConcParam' : None, #mean std min max
        'surfConc' : None, # 0
        'surfRivDist' : None, # 50
        'surfRivInfl' : None, # 200
        

        # -----------------
        # Filter options
        # -----------------
        'RadSelect' : None,
        'RectSelect' : None,
        'DepthRange' : None,
        'UnsatRange' : None,
        'Wt2tRange' : None,
        'ScreenLenRange' : None,
        'SourceArea' : None,

        # -----------------
        # Raster loading Options
        # -----------------
        'loadSubScen' : None,
        'modifierName' : None,
        'modifierType' : None, # 'Multiply' or 'Replace'
        'modifierUnit' : None, # 'CONC' or 'MASS'

        # -----------------
        # Misc Options
        # -----------------
        'getids' : None,
        'getOtherInfo' : None,
        'getotherinfo' : None, # Server spelling; overrides getOtherInfo if set.
        'DebugID' : None,
        'printLF' : None,
        'printURF' : None,
        'printBTC' : None,
        'printWellBTC' : None,
        'PixelRadius' : None,
        'maxSourceCells' : None
    }

    return scenario


def run_mantis(scenario, return_metadata=False):
    """Run the client and return (btc, success).

    With return_metadata=True, return (btc, well_info, success), separating
    metadata according to the scenario's getids and getOtherInfo settings.
    The default preserves the legacy array, including any metadata columns.
    """
    btc = []
    tf = 0
    if scenario['client'] is None:
        msg = "No client was specified. Returning empty lists."
        print(f"\033[93mWARNING:\033[0m {msg}")  # 93 = yellow
        if return_metadata:
            return np.empty((0, 0)), None, False
        return btc, tf

    # Delete the output file if it already exists
    outfile = scenario["outfile"]
    if os.path.exists(outfile):
        os.remove(outfile)

    write_mantis_input(scenario)

    cmd = f"{scenario['client']} {scenario['infile']} {scenario['outfile']}"
    # cross-platform null device
    with open(os.devnull, 'w') as null:
        subprocess.call(cmd, stdout=null, stderr=null)

    other_info = scenario.get("getotherinfo")
    if other_info is None:
        other_info = scenario.get("getOtherInfo")
    return read_mantis_output(
        scenario["outfile"],
        getids=bool(scenario.get("getids")),
        getOtherInfo=bool(other_info),
        return_metadata=return_metadata,
    )

def write_mantis_input(scenario: dict):
    """Write client input using Scenario.h option names without mutating scenario.

    Set denitMap, denitThreshold, and
    denitFractionRemoved together. These serialize to the server's
    single denit option. Setting denitMap=None disables denitrification,
    ignoring retained threshold/fraction values and any legacy denit list.
    The legacy denit list [map, threshold, fraction] is accepted only when
    the denitMap key is absent.
    The historic raster
    subset uses LoadTransitionSubScen (loadTransitionSubScen is also accepted).
    Boolean output flags are serialized as integers for the C++ parser.
    """
    denit = scenario.get("denit")
    if "denitMap" in scenario:
        denit = None
        if scenario["denitMap"] is not None:
            denit = [scenario.get(key) for key in (
                "denitMap", "denitThreshold", "denitFractionRemoved",
            )]
            if any(value is None for value in denit):
                raise ValueError(
                    "Set denitThreshold and denitFractionRemoved when denitMap is enabled."
                )
    if denit is not None:
        if isinstance(denit, (str, bytes)) or len(denit) != 3:
            raise ValueError("denit must contain [map_name, threshold, fraction_removed]")
        if not isinstance(denit[0], str) or not denit[0] or any(c.isspace() for c in denit[0]):
            raise ValueError("denit map_name must be a nonempty name without whitespace")
        if not all(0 <= float(value) <= 1 for value in denit[1:]):
            raise ValueError("denit threshold and fraction_removed must be between 0 and 1")

    with open(scenario["infile"], "w") as f:

        # Description comments
        if scenario.get("descr"):
            f.write(f"# { scenario['descr']}\n")

        # -----------------
        # GUI Scenario Options
        # -----------------
        f.write(f"modelArea {scenario['modelArea']}\n")
        f.write(f"bMap {scenario['bMap']}\n")
        f.write(f"Nregions {len(scenario['Regions'])}")
        for region in scenario["Regions"]:
            f.write(f" {region}")
        f.write("\n")
        f.write(f"flowScen {scenario['flowScen']}\n")
        f.write(f"wellType {scenario['wellType']}\n")
        f.write(f"unsatScen {scenario['unsatScen']}\n")
        f.write(f"unsatWC {scenario['unsatWC']}\n")
        f.write(f"endSimYear {scenario['endSimYear']}\n")

        # -----------------
        # loading options
        # -----------------
        f.write(f"loadScen {scenario['loadScen']}\n")
        f.write(f"startRed {scenario['startRed']}\n")
        f.write(f"endRed {scenario['endRed']}\n")
        crops = scenario.get("Crops")
        if crops is None or len(crops) == 0:
            crops = [[-9, 1]]
        f.write(f"Ncrops {len(crops)}\n")
        for row in crops:
            f.write(f"{row[0]} {row[1]:.3f}\n")

        # -----------------
        # Options with defaults
        # -----------------
        if scenario.get("minRch") is not None:
            f.write(f"minRch {scenario['minRch']:.2f}\n")
        if scenario.get("rchMap"):
            f.write(f"rchMap {scenario['rchMap']}\n")
        if scenario.get("maxConc") is not None:
            f.write(f"maxConc {scenario['maxConc']:.2f}\n")
        if scenario.get("maxAge") is not None:
            f.write(f"maxAge {scenario['maxAge']:.2f}\n")
        if scenario.get("constRed") is not None:
            f.write(f"constRed {scenario['constRed']:.3f}\n")
        if scenario.get("LoadTransitionName"):
            transition_start = scenario.get("LoadTransitionStart")
            transition_end = scenario.get("LoadTransitionEnd")
            if transition_start is None:
                transition_start = 2005
            if transition_end is None:
                transition_end = 2015
            f.write(
                f"loadTrans {scenario['LoadTransitionName']} "
                f"{transition_start} {transition_end}\n"
            )
        else:
            f.write(
                f"loadTrans {'NONE'} "
                f"{0} {0}\n"
            )

        if denit is not None:
            f.write(f"denit {denit[0]} {float(denit[1])} {float(denit[2])}\n")
        transition_subset = scenario.get("LoadTransitionSubScen")
        if transition_subset is None:
            transition_subset = scenario.get("loadTransitionSubScen")
        if transition_subset:
            f.write(f"LoadTransitionSubScen {transition_subset}\n")

        if scenario.get("startSimYear") is not None:
            f.write(f"startSimYear {scenario['startSimYear']}\n")
        if scenario.get("por") is not None:
            f.write(f"por {scenario['por']}\n")
        if scenario.get("urfType"):
            f.write(f"urfType {scenario['urfType']}\n")

        if scenario.get("ADELR") is not None:
            f.write(f"ADELR {scenario['ADELR'][0]:.5e} {scenario['ADELR'][1]:.5e}\n")

        if scenario.get("initConcParam") is not None:
            vals = scenario["initConcParam"]
            f.write(f"initConcParam {vals[0]:.5e} {vals[1]:.5e} {vals[2]:.5e} {vals[3]:.5e}\n")

        if scenario.get("surfConc") is not None:
            f.write(f"surfConc {scenario['surfConc']:.2f}\n")

        if scenario.get("surfRivDist") is not None:
            f.write(f"surfRivDist {scenario['surfRivDist']:.2f}\n")

        if scenario.get("surfRivInfl") is not None:
            f.write(f"surfRivInfl {scenario['surfRivInfl']:.2f}\n")

        # -----------------
        # Filter Options
        # -----------------
        for key, fmt in [
            ("RadSelect", "RadSelect {:.2f} {:.2f} {:.2f}\n"),
            ("RectSelect", "RectSelect {:.2f} {:.2f} {:.2f} {:.2f}\n"),
            ("DepthRange", "DepthRange {:.2f} {:.2f}\n"),
            ("UnsatRange", "UnsatRange {:.2f} {:.2f}\n"),
            ("Wt2tRange", "Wt2tRange {:.2f} {:.2f}\n"),
            ("ScreenLenRange", "ScreenLenRange {:.2f} {:.2f}\n"),
        ]:
            if scenario.get(key) is not None:
                f.write(fmt.format(*scenario[key]))

        if scenario.get("SourceArea") is not None:
            f.write(
                f"SourceArea {scenario['SourceArea'][0]} {scenario['SourceArea'][1]} "
                f"{scenario['SourceArea'][2]} {scenario['SourceArea'][3]:.2f}\n"
            )

        # -----------------
        # Raster loading Options
        # -----------------
        for key in ["loadSubScen", "modifierName", "modifierType", "modifierUnit"]:
            if scenario.get(key):
                f.write(f"{key} {scenario[key]}\n")

        # -----------------
        # Misc Options
        # -----------------
        if scenario.get("getids") is not None:
            f.write(f"getids {int(bool(scenario['getids']))}\n")
        other_info = scenario.get("getotherinfo")
        if other_info is None:
            other_info = scenario.get("getOtherInfo")
        if other_info is not None:
            f.write(f"getotherinfo {int(bool(other_info))}\n")

        if scenario.get("DebugID"):
            f.write(f"DebugID {scenario['DebugID']}\n")
        for key in ("printLF", "printURF", "printBTC", "printWellBTC"):
            if scenario.get(key) is not None:
                f.write(f"{key} {int(bool(scenario[key]))}\n")

        if scenario.get("PixelRadius") is not None:
            f.write(f"PixelRadius {scenario['PixelRadius']}\n")
        if scenario.get("maxSourceCells") is not None:
            f.write(f"maxSourceCells {scenario['maxSourceCells']}\n")


def read_mantis_output(filename, *, getids=False, getOtherInfo=False,
                       return_metadata=False):
    """Read a Mantis client output file.

    By default return (data, success), preserving all columns for existing
    callers. With return_metadata=True return (btc, well_info, success): btc
    contains only concentrations and well_info is a pandas DataFrame with
    matching row order and integer columns well_id, count_rivers, count_max_age
    (the counts are included only with getOtherInfo). Without getids,
    well_info is None and getOtherInfo has no effect, matching the server.

    Pass the same flags used to generate the file; the header reports total
    columns and cannot distinguish metadata from simulation years. For an
    empty result, metadata mode returns an empty 2D array, None, and False.
    """
    with open(filename, "r") as f:
        # first line: Nbtc and Nyrs
        first_line = f.readline().split()
        Nbtc, Nyrs = int(first_line[0]), int(first_line[1])

        if Nbtc == 0:
            if return_metadata:
                return np.empty((0, 0)), None, False
            return np.array([]), False

        # read remaining numbers
        data = np.fromfile(f, sep=" ", count=Nbtc * Nyrs)

        # reshape to Nyrs x Nbtc then transpose -> Nbtc x Nyrs
        btc = data.reshape((Nbtc,Nyrs))

        if return_metadata:
            well_info = None
            if getids:
                import pandas as pd

                columns = ["well_id"]
                if getOtherInfo:
                    columns += ["count_rivers", "count_max_age"]
                n_metadata = len(columns)
                if Nyrs < n_metadata:
                    raise ValueError("Output has fewer columns than the requested metadata")
                metadata = btc[:, :n_metadata]
                if not np.all(np.isfinite(metadata) & (metadata == np.floor(metadata))):
                    raise ValueError("Well IDs and streamline counts must be finite integers")
                well_info = pd.DataFrame(metadata.astype(np.int64), columns=columns)
                btc = btc[:, n_metadata:]
            return btc, well_info, True

        return btc, True

@dataclass
class gridSpec:
    cornerX: float
    cornerY: float
    cellSize: int
    Nrows: int
    Ncols: int

def mantis_grid_spec(name):
    gs = gridSpec(cornerX = np.nan,
          cornerY = np.nan,
          cellSize = -9,
          Nrows = -9,
          Ncols = -9)

    if name == "CentralValley":
        gs.cornerX = -223300
        gs.cornerY = -344600
        gs.cellSize = 50
        gs.Nrows = 12863
        gs.Ncols = 7046
    elif name == "TuleRiver":
        gs.cornerX = 41800
        gs.cornerY = -247550
        gs.cellSize = 50
        gs.Nrows = 1027
        gs.Ncols = 1068
    elif name == "Modesto":
        gs.cornerX = -111800
        gs.cornerY = -83200
        gs.cellSize = 50
        gs.Nrows = 1477
        gs.Ncols = 1428

    return gs

def calc_grid_row_col(P, gs: gridSpec):
    XG = gs.cornerX + np.concatenate((
        [0],
        np.cumsum(gs.cellSize * np.ones(gs.Ncols))
    ))
    YG = gs.cornerY + np.concatenate((
        [0],
        np.cumsum(gs.cellSize * np.ones(gs.Nrows))
    ))

    Np = len(P)
    ind = np.arange(Np)
    R = np.full(Np, np.nan)
    C = np.full(Np, np.nan)
    P_bck = P.copy()

    # Process X coordinates
    mask = P['X'] < XG[0]
    if mask.any():
        C[ind[mask],] = 1
        P = P[~mask]
        ind = ind[~mask]

    mask = P['X'] > XG[-1]
    if mask.any():
        C[ind[mask],] = len(XG)  # or 1 if matching MATLAB 1-based behavior
        P = P[~mask]
        ind = ind[~mask]

    for ii in range(len(XG) - 1):
        mask = (P['X'] >= XG[ii]) & (P['X'] <= XG[ii + 1])
        if mask.any():
            C[ind[mask],] = ii + 1  # +1 to match MATLAB 1-based index
            P = P[~mask]
            ind = ind[~mask]

    # Restore
    P = P_bck
    ind = np.arange(Np)
    Ny = len(YG) - 1

    # Process Y coordinates
    mask = P['Y'] < YG[0]
    if mask.any():
        R[ind[mask],] = Ny
        P = P[~mask]
        ind = ind[~mask]

    mask = P['Y'] > YG[-1]
    if mask.any():
        R[ind[mask],] = 1
        P = P[~mask]
        ind = ind[~mask]

    for ii in range(Ny):
        mask = (P['Y'] >= YG[ii]) & (P['Y'] <= YG[ii + 1])
        if mask.any():
            R[ind[mask],] = Ny - ii  # Match MATLAB's reversed row index
            P = P[~mask]
            ind = ind[~mask]

    # Restore
    P = P_bck
    P['urfI'] = R
    P['urfJ'] = C

    return P

def write_linear_data(names, data, filename, suppress_print=False):
    """Write (cell, scenario) values in Mantis LinearData format.

    HDF5 stores Data as (scenario, cell); ASCII stores one cell per line.
    A one-dimensional input is accepted for a single named scenario.
    Runtime loading requires HDF5 with exactly one scenario and values in
    the server raster's active-cell order.
    Set suppress_print=True to suppress the HDF5 success message.
    """
    names = list(names)
    data = np.asarray(data)

    # Validate dimensions
    n = len(names)
    if n == 0 or any(not isinstance(name, str) or not name for name in names) or len(set(names)) != n:
        raise ValueError("names must contain unique, nonempty scenario names.")
    if data.ndim == 1 and n == 1:
        data = data.reshape(-1, 1)
    if data.ndim != 2:
        raise ValueError("data must have shape (cells, scenarios).")
    if data.shape[1] != n:
        raise ValueError(f"Data must have {n} columns, but has {data.shape[1]}.")

    # Check file extension
    ext = os.path.splitext(filename)[1].lower()

    if ext in (".h5", ".hdf5"):
        # Save as HDF5
        with h5py.File(filename, "w") as f:
            f.create_dataset("Names", data=np.array(names, dtype=h5py.string_dtype()))
            f.create_dataset("Data", data=data.T)
        if not suppress_print:
            print(f"HDF5 file written: {filename}")
    else:
        # Save as ASCII
        with open(filename, "w") as f:
            f.write(f"{n}\n")
            for name in names:
                f.write(f"{name}\n")
            np.savetxt(f, data, fmt="%.3f")
        print(f"ASCII file written: {filename}")

def read_linear_data(filename):
    """Read names and (cell, scenario) values from a Mantis LinearData file.

    HDF5 uses the server's (scenario, cell) layout. Legacy Python files with
    (cell, scenario) layout are also accepted when the shape is unambiguous.
    """
    ext = os.path.splitext(filename)[1].lower()

    if ext in (".h5", ".hdf5"):
        # Read HDF5
        with h5py.File(filename, "r") as f:
            names = list(f["Names"].asstr()[()])
            data = np.array(f["Data"])
            if data.ndim != 2:
                raise ValueError("HDF5 Data must be two-dimensional.")
            if data.shape[0] == len(names):
                data = data.T
            elif data.shape[1] != len(names):
                raise ValueError("HDF5 Data dimensions do not match Names.")
    else:
        # Read ASCII
        with open(filename, "r") as f:
            # First line: number of names
            n = int(f.readline().strip())
            # Next n lines: names
            names = [f.readline().strip() for _ in range(n)]
            # Remaining: array
            data = np.loadtxt(f, ndmin=2)

        # Ensure shape consistency (convert to 2D even if 1 row)
        if data.ndim == 1:
            data = data.reshape(1, -1)

    return names, data

def read_raster(filename):
    """Read raster indices as (N, 2), accepting either HDF5 orientation.

    A (2, 2) HDF5 dataset is treated as already having two columns.
    """
    ext = os.path.splitext(filename)[1].lower()
    if ext == ".h5":

        with h5py.File(filename, "r") as f:
            if 'Raster' not in f:
                raise KeyError("Dataset 'Raster' not found in HDF5 file.")
            data = np.array(f["Raster"])
            if data.ndim != 2:
                raise ValueError("HDF5 Raster dataset must be two-dimensional.")
            if data.shape[1] != 2:
                if data.shape[0] != 2:
                    raise ValueError("HDF5 Raster dataset must have shape (N, 2) or (2, N).")
                data = data.T

    else:
        data = np.loadtxt(filename)
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError("ASCII file must contain exactly two columns.")

    return data

def linear_to_raster(linear_data, idx, nrow, ncol):
    """Place (N,) or (N, 1) values at zero-based (N, 2) row/column indices.

    Return a (nrow, ncol) array with NaN in cells not listed in idx.
    Index pairs should be unique for a lossless round trip.
    """
    linear_data = np.asarray(linear_data)
    idx = np.asarray(idx)
    if linear_data.ndim == 2 and linear_data.shape[1] == 1:
        linear_data = linear_data[:, 0]
    if linear_data.ndim != 1:
        raise ValueError("linear_data must have shape (N,) or (N, 1).")
    if idx.ndim != 2 or idx.shape[1] != 2:
        raise ValueError("idx must have shape (N, 2).")
    if len(linear_data) != len(idx):
        raise ValueError("linear_data and idx must have the same number of rows.")
    if not np.issubdtype(idx.dtype, np.integer):
        raise ValueError("idx must contain integer row/column indices.")
    if np.any(idx < 0) or np.any(idx[:, 0] >= nrow) or np.any(idx[:, 1] >= ncol):
        raise ValueError("idx contains indices outside the raster dimensions.")
    out = np.full((nrow, ncol), np.nan)
    out[idx[:, 0], idx[:, 1]] = linear_data
    return out

def raster_to_linear(raster, idx):
    """Extract raster values as (N,) in the order of the (N, 2) idx array.

    idx contains zero-based [row, column] pairs, as in linear_to_raster.
    Values (including negatives and NaNs) and their dtype are preserved.
    Use raster_to_linear(raster, idx).reshape(-1, 1) for an (N, 1) result.
    """
    raster = np.asarray(raster)
    idx = np.asarray(idx)
    if raster.ndim != 2:
        raise ValueError("raster must be a two-dimensional array.")
    if idx.ndim != 2 or idx.shape[1] != 2:
        raise ValueError("idx must have shape (N, 2).")
    if not np.issubdtype(idx.dtype, np.integer):
        raise ValueError("idx must contain integer row/column indices.")
    if (np.any(idx < 0) or np.any(idx[:, 0] >= raster.shape[0])
            or np.any(idx[:, 1] >= raster.shape[1])):
        raise ValueError("idx contains indices outside the raster dimensions.")
    return raster[idx[:, 0], idx[:, 1]]

def convert_raster_to_column(raster):
    rows, cols = np.where(raster >= 0)
    values = raster[rows, cols]

    order = np.argsort(values)

    return np.column_stack((rows[order], cols[order]))

def convert_column_to_raster(raster_column, shape, fill_value=-1):
    raster = np.full(shape, fill_value, dtype=int)
    for i, (r, c) in enumerate(raster_column):
        raster[r, c] = i
    return raster

def convolute_urf(urf, lf, init_conc=0.0):
    """
    Convolve URFs with loading functions and compute pre-existing concentration BTC.

    Parameters
    ----------
    urf : array_like, shape (Ns, Nt)
        Unit response functions.

    lf : array_like, shape (Ns, Nt) or (Nt,) or (1, Nt)
        Loading functions. If one loading function is provided, it is applied
        to all URFs.

    init_conc : float, default 0.0
        Initial/background concentration used to compute prebtc.

    Returns
    -------
    btc : ndarray, shape (Ns, Nt_lf)
        Breakthrough curve from loading functions.

    prebtc : ndarray, shape (Ns, Nt_lf)
        Contribution from initial concentration:
        prebtc = init_conc * (1 - cumulative_urf)
    """

    urf = np.asarray(urf, dtype=float).copy()
    lf = np.asarray(lf, dtype=float)

    if urf.ndim != 2:
        raise ValueError("urf must be a 2D array with shape (Ns, Nt_urf).")

    if lf.ndim == 1:
        lf = lf[None, :]

    if lf.ndim != 2:
        raise ValueError("lf must be a 1D or 2D array.")

    Ns, Nt_urf = urf.shape
    Ns_lf, Nt_lf = lf.shape

    if Ns_lf == 1:
        lf = np.repeat(lf, Ns, axis=0)
    elif Ns_lf != Ns:
        raise ValueError(
            "The number of URFs must equal the number of loading functions, "
            "unless only one loading function is provided."
        )

    urf[urf < 0] = 0.0

    btc = np.zeros((Ns, Nt_lf), dtype=float)
    prebtc = np.zeros((Ns, Nt_lf), dtype=float)

    n_conv = min(Nt_urf, Nt_lf)

    for i in range(Ns):
        # Loading contribution
        full_conv = np.convolve(lf[i, :], urf[i, :], mode="full")
        btc[i, :] = full_conv[:Nt_lf]

        # Initial concentration contribution
        cumurf = np.cumsum(urf[i, :n_conv])
        prebtc[i, :n_conv] = init_conc * (1.0 - cumurf)

        if Nt_lf > n_conv:
            prebtc[i, n_conv:] = init_conc * (1.0 - cumurf[-1])

    return btc, prebtc


def simulate_feedback_concentrations(
    VIwells, VIurfs, VImsa, scenario, viparam

):
    """
    Simulate salinity feedback for all wells.

    Parameters
    ----------
    VIwells : DataFrame
        Must contain field 'Eid'.

    VIurfs : DataFrame
        Must contain fields:
            Eid
            v_cds

        Rows correspond one-to-one with VImsa rows.

    VImsa : ndarray
        MSA matrix.
        Column 6 = mu
        Column 7 = sigma

    scenario : dict
        {
            'applied_mm': 200,
            'perc_ratio': 1/3,
            'init_conc': 200,
            'n_years': 200
        }

    viparam : dict
        {
            'm_col':6
            's_col:'7
        }

    Returns
    -------
    applied_conc_all : ndarray
        Shape (Nwells, Nyears)

    perc_conc_all : ndarray
        Shape (Nwells, Nyears)
    """

    applied_mm = scenario["applied_mm"]
    perc_ratio = scenario["perc_ratio"]
    init_conc = scenario["init_conc"]
    n_years = scenario.get("n_years", 200)
    m_col = viparam['m_col']
    s_col = viparam['s_col']
    w_field = viparam['w_field']

    perc_mm = applied_mm * perc_ratio

    x = np.arange(1, n_years + 1)

    eids = VIurfs["Eid"].to_numpy()

    # Build lookup once
    eid_to_rows = {}
    for i, eid in enumerate(eids):
        eid_to_rows.setdefault(eid, []).append(i)

    n_wells = len(VIwells)

    applied_conc_all = np.zeros((n_wells, n_years))
    perc_conc_all = np.zeros((n_wells, n_years))

    for iw, well in enumerate(tqdm(VIwells.itertuples(),
                                   total=len(VIwells),
                                   desc="Simulating wells")):

        rows = eid_to_rows.get(well.Eid)

        if rows is None:
            continue

        rows = np.asarray(rows)

        msa = VImsa[rows, :]
        urfs = VIurfs.iloc[rows]

        n_paths = len(rows)

        urfs_mat = np.zeros((n_paths, n_years))

        for ip in range(n_paths):

            mu_ln = msa[ip, m_col]
            sigma_ln = msa[ip, s_col]

            if (
                not np.isfinite(mu_ln)
                or not np.isfinite(sigma_ln)
                or mu_ln <= 0
                or sigma_ln <= 0
            ):
                urfs_mat[ip, :] = 0.0
                continue

            urfs_mat[ip, :] = lognorm.pdf(
                x,
                s=sigma_ln,
                scale=np.exp(mu_ln),
            )

        wgh = urfs[w_field].to_numpy(dtype=float)

        if wgh.sum() <= 0:
            continue

        wgh /= wgh.sum()

        lf = np.zeros(n_years)

        for iyr in range(n_years):

            btc, prebtc = convolute_urf(
                urfs_mat,
                lf[: iyr + 1],
                init_conc,
            )

            tot_conc = btc[:, -1] + prebtc[:, -1]

            well_conc = np.sum(tot_conc * wgh)

            applied_conc_all[iw, iyr] = well_conc

            perc_conc = (
                well_conc
                * applied_mm
                / perc_mm
            )

            perc_conc_all[iw, iyr] = perc_conc

            lf[iyr] = perc_conc

    return applied_conc_all, perc_conc_all


def expand_urfs(urf_df, m_field, s_field, n_years=100):
    """
    Expand URF parameters into annual lognormal response functions.

    Parameters
    ----------
    urf_df : pandas.DataFrame
        DataFrame containing the URF mean and standard deviation fields.
    m_field : str
        Column name containing the mean.
    s_field : str
        Column name containing the standard deviation.
    n_years : int, default=100
        Number of years over which to expand the URFs.

    Returns
    -------
    urfs : np.ndarray
        Array of shape (len(urf_df), n_years).
    """

    n = len(urf_df)
    urfs = np.zeros((n, n_years), dtype=float)

    x = np.arange(1, n_years + 1, dtype=float)

    # Special URF for m = -1
    urf_m1 = np.array([
        0.48443252, 0.41642340, 0.08307405, 0.01338364,
        0.00219913, 0.00039049, 0.00007584, 0.00001608,
        0.00000370, 0.00000092, 0.00000023
    ])

    # Special URF for m = -2
    urf_m2 = np.array([
        0.00095171, 0.19513243, 0.41957050, 0.25244126,
        0.09333424, 0.02803643, 0.00771540, 0.00206063,
        0.00055016, 0.00014916, 0.00004141, 0.00001182,
        0.00000347, 0.00000106, 0.00000032
    ])

    for i, (m, s) in enumerate(
        zip(urf_df[m_field].to_numpy(),
            urf_df[s_field].to_numpy())
    ):

        # Special predefined URFs
        if m == -1:
            ncopy = min(len(urf_m1), n_years)
            urfs[i, :ncopy] = urf_m1[:ncopy]

        elif m == -2:
            ncopy = min(len(urf_m2), n_years)
            urfs[i, :ncopy] = urf_m2[:ncopy]

        # Any other negative mean -> all zeros
        elif m < 0:
            continue

        # Lognormal URF
        else:
            urfs[i, :] = lognorm.pdf(
                x,
                s=s,
                scale=np.exp(m)
            )

    return urfs
