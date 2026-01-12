import os
import subprocess
import numpy as np
import h5py
from dataclasses import dataclass

def mantis_inputs():

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
        'loadScen' : 'GNLM',
        'startRed' : 2020,
        'endRed' : 2030,
        'Crops' : None,

        # -----------------
        # options that take default values
        # -----------------
        'minRch' : None, # 10;
        'rchMap' : None,
        'maxConc' : None, #250;
        'maxAge' : None, #400;
        'constRed' : None, #1.0;
        'LoadTransitionName' : None, #'GNLM';
        'LoadTransitionStart' : None, #2035;
        'LoadTransitionEnd' : None, #2025;
        'startSimYear' : None, #1945;
        'por' : None, #10;
        'urfType' : None, #'LGNRM' or 'ADE';
        'ADELR' : None,
        'initConcParam' : None, #mean std min max

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
        'modifierType' : None,
        'modifierUnit' : None,

        # -----------------
        # Misc Options
        # -----------------
        'getids' : None,
        'getOtherInfo' : None,
        'DebugID' : None,
        'printLF' : None,
        'printURF' : None,
        'printBTC' : None,
        'printWellBTC' : None,
        'PixelRadius' : None,
        'maxSourceCells' : None
    }

    return scenario


def run_mantis(scenario):
    btc = []
    tf = 0
    if scenario['client'] is None:
        msg = "No client was specified. Returning empty lists."
        print(f"\033[93mWARNING:\033[0m {msg}")  # 93 = yellow
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

    btc, tf = read_mantis_output(scenario["outfile"])

    return btc, tf

def write_mantis_input(scenario: dict):
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
        f.write(f"unsatWC {scenario['unsatWC']:.2f}\n")
        f.write(f"endSimYear {scenario['endSimYear']}\n")

        # -----------------
        # loading options
        # -----------------
        f.write(f"loadScen {scenario['loadScen']}\n")
        f.write(f"startRed {scenario['startRed']}\n")
        f.write(f"endRed {scenario['endRed']}\n")
        if not scenario.get("Crops"):
            scenario["Crops"] = [[-9, 1]]
        f.write(f"Ncrops {len(scenario['Crops'])}\n")
        for row in scenario["Crops"]:
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
            if not scenario.get("LoadTransitionStart") or not scenario.get("LoadTransitionEnd"):
                scenario["LoadTransitionStart"] = 2005
                scenario["LoadTransitionEnd"] = 2015
            f.write(
                f"loadTrans {scenario['LoadTransitionName']} "
                f"{scenario['LoadTransitionStart']} {scenario['LoadTransitionEnd']}\n"
            )
        else:
            f.write(
                f"loadTrans {'NONE'} "
                f"{0} {0}\n"
            )

        if scenario.get("startSimYear"):
            f.write(f"startSimYear {scenario['startSimYear']}\n")
        if scenario.get("por"):
            f.write(f"por {scenario['por']}\n")
        if scenario.get("urfType"):
            f.write(f"urfType {scenario['urfType']}\n")

        if scenario.get("ADELR") is not None:
            f.write(f"ADELR {scenario['ADELR'][0]:.5e} {scenario['ADELR'][1]:.5e}\n")

        if scenario.get("initConcParam") is not None:
            vals = scenario["initConcParam"]
            f.write(f"initConcParam {vals[0]:.5e} {vals[1]:.5e} {vals[2]:.5e} {vals[3]:.5e}\n")

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
        if scenario.get("getids"):
            if scenario["getids"] != 0:
                f.write("getids 1\n")
                if scenario.get("getOtherInfo") and scenario["getOtherInfo"] != 0:
                    f.write("getotherinfo 1\n")

        if scenario.get("DebugID"):
            f.write(f"DebugID {scenario['DebugID']}\n")
            if scenario.get("printLF"):
                f.write(f"printLF {scenario['printLF']}\n")
            if scenario.get("printURF"):
                f.write(f"printURF {scenario['printURF']}\n")
            if scenario.get("printBTC"):
                f.write(f"printBTC {scenario['printBTC']}\n")
            if scenario.get("printWellBTC"):
                f.write(f"printWellBTC {scenario['printWellBTC']}\n")

        if scenario.get("PixelRadius"):
            f.write(f"PixelRadius {scenario['PixelRadius']}\n")
        if scenario.get("maxSourceCells"):
            f.write(f"maxSourceCells {scenario['maxSourceCells']}\n")


def read_mantis_output(filename):
    with open(filename, "r") as f:
        # first line: Nbtc and Nyrs
        first_line = f.readline().split()
        Nbtc, Nyrs = int(first_line[0]), int(first_line[1])

        if Nbtc == 0:
            return np.array([]), False

        # read remaining numbers
        data = np.fromfile(f, sep=" ", count=Nbtc * Nyrs)

        # reshape to Nyrs x Nbtc then transpose -> Nbtc x Nyrs
        btc = data.reshape((Nbtc,Nyrs))

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

def write_linear_data(names, data, filename):
    names = list(names)
    data = np.asarray(data)

    # Validate dimensions
    n = len(names)
    if data.shape[1] != n:
        raise ValueError(f"Data must have {n} columns, but has {data.shape[1]}.")

    # Check file extension
    ext = os.path.splitext(filename)[1].lower()

    if ext == ".h5":
        # Save as HDF5
        with h5py.File(filename, "w") as f:
            f.create_dataset("Names", data=np.array(names, dtype=h5py.string_dtype()))
            f.create_dataset("Data", data=data)
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
    ext = os.path.splitext(filename)[1].lower()

    if ext == ".h5":
        # Read HDF5
        with h5py.File(filename, "r") as f:
            names = [str(n) for n in f["Names"][()]]
            data = np.array(f["Data"])
    else:
        # Read ASCII
        with open(filename, "r") as f:
            # First line: number of names
            n = int(f.readline().strip())
            # Next n lines: names
            names = [f.readline().strip() for _ in range(n)]
            # Remaining: array
            data = np.loadtxt(f)

        # Ensure shape consistency (convert to 2D even if 1 row)
        if data.ndim == 1:
            data = data.reshape(1, -1)

    return names, data