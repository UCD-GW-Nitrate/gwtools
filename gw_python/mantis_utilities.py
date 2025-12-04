import numpy as np
import h5py
from dataclasses import dataclass

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