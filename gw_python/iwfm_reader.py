import sys
import pandas as pd
import numpy as np
import re
import h5py

def _read_all_file(filename):
    # Read all lines
    try:
        with open(filename, "r") as f:
            all_lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None

    return all_lines


def read_nodes(filename):
    # Read all lines
    all_lines = _read_all_file(filename)
    if all_lines is None:
        return None

    line_index = 0
    tmp, line_index = _read_next_data_block(all_lines, line_index, 'scalar', 2, float)
    nd = int(tmp[0])
    #fact = tmp[1]

    node_data, line_index = _read_next_data_block(all_lines, line_index, 'array', nd, float)

    node_df = pd.DataFrame(node_data, columns=["ID", "X", "Y"])
    node_df["ID"] = node_df["ID"].astype(int)

    return node_df

def _read_next_data_block(all_lines, line_index, data_format, n_data, data_type, C = 'C'):
    out = []
    while line_index < len(all_lines):
        line = all_lines[line_index].strip()
        # Skip comment lines ('C') and empty lines
        if line.startswith(C) or not line:
            line_index += 1
            continue

        try:
            if data_format == 'scalar':
                for ii in range(n_data):
                    data = data_type(line.split()[0])
                    out.append(data)
                    line_index += 1
                    if line_index >= len(all_lines):
                        break
                    line = all_lines[line_index].strip()
            elif data_format == 'array':
                for ii in range(n_data):
                    data = [data_type(tok) for tok in line.split()]
                    out.append(data)
                    line_index += 1
                    if line_index >= len(all_lines):
                        break
                    line = all_lines[line_index].strip()
            break
        except (ValueError, IndexError):
            # This handles separator lines or unexpected text
            print(f"Error: WHILE READING ['{line}'].")
            return None
    return out, line_index


def read_elements(filename):
    all_lines = _read_all_file(filename)
    if all_lines is None:
        return None

    line_index = 0
    tmp, line_index = _read_next_data_block(all_lines, line_index, 'scalar', 2, int)
    ne = tmp[0]
    nregn = tmp[1]

    reg_names, line_index = _read_next_data_block(all_lines, line_index, 'array', nregn, str)

    elem_data, line_index = _read_next_data_block(all_lines, line_index, 'array', ne, int)

    elem_df = pd.DataFrame(elem_data, columns=["IE", "ND1", "ND2", "ND3", "ND4", "IRGE"])
    return elem_df


def read_stratigraphy(filename, n_nodes):
    all_lines = _read_all_file(filename)
    if all_lines is None:
        return None

    line_index = 0
    tmp, line_index = _read_next_data_block(all_lines, line_index, 'scalar', 2, float)
    nlay = int(tmp[0])
    fact = tmp[1]

    elev_data, line_index = _read_next_data_block(all_lines, line_index, 'array', n_nodes, float)

    cols = ['ID', 'GSE'] + [f"{prop}{ilay}" for ilay in range(1, nlay+1) for prop in ['A', 'L']]
    elev_df = pd.DataFrame(elev_data, columns=cols)
    elev_df["ID"] = elev_df["ID"].astype(int)
    return elev_df



def read_VELOUTFILE(filename: str, Nstep: int, Nelem: int, Nlay: int):
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()

    # Total velocities per element:
    nvals = 3 * Nlay
    # Pre-allocate arrays
    CC = np.zeros((Nelem, 2), dtype=float)
    VX = np.zeros((Nelem * Nlay, Nstep), dtype=float)
    VY = np.zeros((Nelem * Nlay, Nstep), dtype=float)
    VZ = np.zeros((Nelem * Nlay, Nstep), dtype=float)
    YMD = np.zeros((Nstep, 3), dtype=int)  # year, month, day

    section = 1
    idx = 0  # position in lines
    timestep = 0

    # ---------------------------------------------------------
    # Parse file
    # ---------------------------------------------------------
    while idx < len(lines):
        line = lines[idx].strip()
        if not line or line.startswith("*"):
            idx += 1
            continue

        # -----------------------------------------------------
        # SECTION 1 — element centroids
        # -----------------------------------------------------
        if section == 1:
            for i in range(Nelem):
                tokens = lines[idx].split()
                # last 2 are X,Y
                CC[i, 0] = float(tokens[-2])
                CC[i, 1] = float(tokens[-1])
                idx += 1
            section = 2
            continue

        # -----------------------------------------------------
        # SECTION 2 — velocity blocks
        # First line contains timestamp + first element
        # Next (Nelem-1) lines contain element lines only
        # -----------------------------------------------------

        if section == 2 and timestep < Nstep:

            # -------- Parse first line of timestep (timestamp + element 1) --------
            L = lines[idx].split()

            # timestamp format e.g. "10/31/1973_24:00"
            date_token = L[0]
            # split into date and maybe hour
            print(f"Reading time step: {date_token}")
            date_part = date_token.split("_")[0]  # "10/31/1973"
            mm, dd, yyyy = date_part.split("/")
            YMD[timestep, :] = [int(yyyy), int(mm), int(dd)]

            # Identify element number (last token that is integer)
            # MATLAB logic: tmp(1) is element id
            elm = int(float(L[-(nvals + 1)]))
            vals = np.array(L[-nvals:], dtype=float)

            # store layer-wise
            base_row = elm - 1
            for lay in range(Nlay):
                row = lay * Nelem + base_row
                VX[row, timestep] = vals[3 * lay]
                VY[row, timestep] = vals[3 * lay + 1]
                VZ[row, timestep] = vals[3 * lay + 2]

            idx += 1

            # -------------------------------
            # Parse remaining Nelem-1 lines
            # -------------------------------
            for _ in range(Nelem - 1):
                arr = np.fromstring(lines[idx], sep=' ')
                elm = int(arr[0])
                vals = arr[1:1 + nvals]

                base_row = elm - 1
                for lay in range(Nlay):
                    row = lay * Nelem + base_row
                    VX[row, timestep] = vals[3 * lay]
                    VY[row, timestep] = vals[3 * lay + 1]
                    VZ[row, timestep] = vals[3 * lay + 2]

                idx += 1

            timestep += 1
            continue
        idx += 1

    return {"CC": CC, "VX": VX, "VY": VY, "VZ": VZ}

def read_gw_head(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    # --- Identify timestamp lines (they start with MM/DD/YYYY_24:00) ---
    #ts_pattern = re.compile(r"^\s*(\d{2})/(\d{2})/(\d{4})_24:00")
    ts_pattern = re.compile(r"^\s*(\d{2})/(\d{2})/(\d{4})_24:00(.*)$")
    #ts_indices = [i for i, ln in enumerate(lines) if ts_pattern.match(ln)]
    ts_indices = []
    timestamps = []
    for i, ln in enumerate(lines):
        m = ts_pattern.match(ln)
        if m:
            ts_indices.append(i)
            mm, dd, yyyy, ss = m.groups()
            timestamps.append((int(dd), int(mm), int(yyyy)))

    Nsteps = len(ts_indices)
    if Nsteps == 0:
        raise ValueError("No valid time-step lines found.")

    timestamps = pd.DataFrame(timestamps, columns=["D", "M", "Y"])

    # Detect Nlay and Nnodes from the FIRST timestep
    first_ts = ts_indices[0]
    first_line = lines[first_ts]
    m = ts_pattern.match(first_line)
    tail = m.group(4).strip()
    if not tail:
        raise ValueError("First timestep line does not contain numeric data.")
    Nnodes = len(tail.split())

    Nlay = 1
    # Auxiliary
    def is_comment(s):
        return s.lstrip().startswith("*")
    #comment = lambda s: s.lstrip().startswith("*")

    j = first_ts + 1
    while j < len(lines):
        ln = lines[j].strip()
        if not ln or is_comment(ln):
            j += 1
            continue
        if ts_pattern.match(ln):
            break
        Nlay += 1
        j += 1

    if Nlay == 0:
        raise ValueError("Could not detect number of layers (Nlay).")

    # --- Allocate output array ---
    heads = np.zeros((Nnodes, Nsteps, Nlay), dtype=float)

    # --- Parse each timestep ---
    for step_idx, idx in enumerate(ts_indices):
        row_start = idx
        lay_counter = 0
        # read next Nlay numeric rows
        while lay_counter < Nlay and row_start < len(lines):
            ln = lines[row_start].strip()

            if ln and not is_comment(ln):
                if lay_counter == 0:
                    m = ts_pattern.match(ln)
                    if m:
                        mm, dd, yyyy, ln = m.groups()
                arr = np.fromstring(ln, sep=' ')
                #arr = np.array([float(x) for x in ln.strip().split()], dtype=float)
                if arr.size != Nnodes:
                    raise ValueError(
                        f"Row {row_start} has {arr.size} values but expected {Nnodes}"
                    )
                heads[:, step_idx, lay_counter] = arr
                lay_counter += 1

            row_start += 1

        if lay_counter != Nlay:
            raise ValueError(f"Time step {step_idx} ended prematurely.")

    return heads, timestamps


def read_iwfm_tecplot_velocity(filename, Nnodes, Nlay, Ntimes, alloc_mult=2.5):
    Ntotal = Nnodes * Nlay
    VX = np.zeros((Ntotal, Ntimes), dtype=np.float64)
    VY = np.zeros((Ntotal, Ntimes), dtype=np.float64)
    VZ = np.zeros((Ntotal, Ntimes), dtype=np.float64)

    max_elem = int(alloc_mult * Nnodes)
    elem = np.zeros((max_elem, 4), dtype=np.int32)
    nodeXY = np.zeros((Nnodes, 2), dtype=np.float64)

    nan_list = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ZONE"):
                break

        n_read = 0
        for line in f:
            if line.startswith("T") or line.startswith("Z"):
                break
                # reached next section prematurely (should not happen)

            parts = line.split()
            if len(parts) < 2:
                continue

            # Parse only X, Y (first two values)
            nodeXY[n_read, 0] = float(parts[0])
            nodeXY[n_read, 1] = float(parts[1])

            n_read += 1
            if n_read >= Nnodes:
                break  # finished reading node coordinates

        # --- Read elements until TEXT ---
        e_read = 0
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue  # skip empty line
            if stripped.startswith("TEXT"):
                break  # end of element list

            parts = stripped.split()
            if len(parts) == 4:
                try:
                    elem[e_read, 0] = int(parts[0])
                    elem[e_read, 1] = int(parts[1])
                    elem[e_read, 2] = int(parts[2])
                    elem[e_read, 3] = int(parts[3])
                    e_read += 1
                except ValueError:
                    continue  # skip malformed numeric lines

                if e_read >= max_elem:
                    raise RuntimeError(
                        f"Element count exceeded allocated size {max_elem}. "
                        f"Increase the multiplier alloc_mult {alloc_mult}."
                    )

        elem = elem[:e_read, :]

        # --- Read velocities for all timesteps ---
        timestep = 0
        n_fixed = 0
        while timestep < Ntimes:
            # Move to next ZONE line
            for line in f:
                if line.startswith("ZONE"):
                    break
            #else:
            #    raise RuntimeError("Unexpected end of file while seeking ZONE.")

            # Read exactly Ntotal*3 floats for this timestep
            values = []
            while len(values) < Ntotal * 3:
                line = f.readline()
                if not line:
                    raise EOFError("Unexpected end of file while reading velocities")
                stripped = line.strip()
                if stripped.startswith(("TEXT", "ZONE")) or not stripped:
                    continue

                #  Parse floats, handle ************
                # This is ok if there are no ************ values
                #values.extend([float(x) for x in stripped.split()])
                for x in stripped.split():
                    try:
                        values.append(float(x))
                    except:
                        values.append(np.nan)
                        # Determine inode, ilay
                        idx = len(values) - 1
                        ilay = idx // Nnodes % Nlay
                        inode = idx % Nnodes
                        nan_list.append((inode, ilay, timestep))

            # For each variable block in order:
            #   VX_L1, VY_L1, VZ_L1, VX_L2, VY_L2, ...
            #
            # We fill:
            #   VX[layer*Nnodes : (layer+1)*Nnodes, t]
            #   VY[...]
            #   VZ[...]
            #
            # BLOCK order is variable-major, not node-major.
            # ----------------------------------------------------------
            values = np.array(values[:Ntotal*3], dtype=np.float64)
            VX[:, timestep] = values[0:Ntotal]
            VY[:, timestep] = values[Ntotal:2*Ntotal]
            VZ[:, timestep] = values[2*Ntotal:3*Ntotal]

            for ilay in range(Nlay):
                start = ilay * Nnodes
                stop = start + Nnodes
                for arr in (VX, VY, VZ):
                    layer_vals = arr[start:stop, timestep]
                    max_val = np.nanmax(layer_vals)
                    nan_mask = np.isnan(layer_vals)
                    n_fix_layer = np.sum(nan_mask)
                    if n_fix_layer > 0:
                        layer_vals[nan_mask] = max_val
                        n_fixed += n_fix_layer

            timestep += 1
            # --- Simple progress bar ---
            bar_len = 40
            fraction = timestep / Ntimes
            filled_len = int(bar_len * fraction)
            bar = "=" * filled_len + "-" * (bar_len - filled_len)
            sys.stdout.write(f"\rReading velocities: [{bar}] {timestep}/{Ntimes}")
            sys.stdout.flush()

        print("\nDone reading velocities.")
        if n_fixed > 0:
            print(f"[{n_fixed}] ************ values found, replaced by layer max.")

    return nodeXY, elem, VX, VY, VZ, nan_list


def read_iwfm_stream_node_budget(filename):
    with h5py.File(filename, "r") as f:
        node_entries = []

        for key in f.keys():
            m = re.match(r"NODE\s+(\d+)$", key)
            if m:
                node_id = int(m.group(1))
                node_entries.append((node_id, key))

        node_entries.sort(key=lambda x: x[0])

        node_ids = []
        data = []
        data_by_node = {}

        for node_id, key in node_entries:
            arr = f[key][()]
            node_ids.append(node_id)
            data.append(arr)
            data_by_node[node_id] = arr


    bud = {
        "NodeIDs": node_ids,
        "Data": data,
        "DataByNode": data_by_node
    }

    return bud

def stream_spec(filename):
    """
    Read IWFM-style stream specification file.

    Skips lines starting with '#' or 'C'.

    Returns
    -------
    streams : dict
        streams[ID] = {
            "ID": ID,
            "IBUR": IBUR,
            "IBDR": IBDR,
            "NAME": NAME,
            "IRV": [...],
            "IGW": [...]
        }

    meta : dict
        {"NHR": NHR, "NRTB": NRTB}
    """

    def valid_lines():
        with open(filename, "r") as f:
            for line in f:
                s = line.strip()

                if not s:
                    continue

                if s.startswith("#") or s.startswith("C"):
                    continue

                yield s

    lines = valid_lines()

    NHR = int(next(lines).split()[0])
    NRTB = int(next(lines).split()[0])

    streams = {}

    for _ in range(NHR):
        # Header line: ID IBUR IBDR NAME
        parts = next(lines).split()

        ID = int(parts[0])
        IBUR = int(parts[1])
        IBDR = int(parts[2])

        # Everything after the first 3 fields is the name
        NAME = " ".join(parts[3:])

        IRV = []
        IGW = []

        for _ in range(IBUR):
            parts = next(lines).split()

            IRV.append(int(parts[0]))
            IGW.append(int(parts[1]))

        streams[ID] = {
            "ID": ID,
            "IBUR": IBUR,
            "IBDR": IBDR,
            "NAME": NAME,
            "IRV": IRV,
            "IGW": IGW
        }

    meta = {
        "NHR": NHR,
        "NRTB": NRTB
    }

    return streams, meta

# Function to read C2VSim Crop areas
def read_iwfm_crop_area(filename, ntimes, nelem, ncrops, nskip, print_progress=False):
    """
    Reads IWFM crop area data per element.

    Python equivalent of:
        CropAreas = readIWFM_CropArea(filename, Ntimes, Nelem, Ncrops, Nskip)

    Returns
    -------
    crop_areas : dict
        crop_areas["Data"] : list of length ntimes
            Each entry is a (nelem, ncrops + 1) numpy array.
        crop_areas["YMD"] : ndarray
            Array of shape (ntimes, 3), with columns [year, month, day].
    """

    with open(filename, "r") as f:
        lines = f.read().splitlines()

    current_line = nskip

    crop_areas = {
        "Data": [None] * ntimes,
        "YMD": np.zeros((ntimes, 3), dtype=int),
    }

    for ii in range(ntimes):
        parts = lines[current_line].split()

        # First token is assumed to be something like MM/DD/YYYY_...
        date_token = parts[0]
        date_part = date_token.split("_")[0]
        month, day, year = map(int, date_part.split("/"))

        if print_progress:
            print(date_token)

        crop_areas["YMD"][ii, :] = [year, month, day]

        # Remaining values on first line
        values = [float(x) for x in parts[1:] if _is_float(x)]

        arr = np.zeros((nelem, ncrops + 1), dtype=float)
        arr[0, :] = values

        current_line += 1

        for j in range(1, nelem):
            parts = lines[current_line].split()

            values = [float(x) for x in parts if _is_float(x)]
            arr[j, :] = values

            current_line += 1

        crop_areas["Data"][ii] = arr

    return crop_areas


def _is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

# Function to read C2VSim Crop areas
def read_iwfm_crop_area(filename, ntimes, nelem, ncrops, nskip, print_progress=False):
    """
    Reads IWFM crop area data per element.

    Python equivalent of:
        CropAreas = readIWFM_CropArea(filename, Ntimes, Nelem, Ncrops, Nskip)

    Returns
    -------
    crop_areas : dict
        crop_areas["Data"] : list of length ntimes
            Each entry is a (nelem, ncrops + 1) numpy array.
        crop_areas["YMD"] : ndarray
            Array of shape (ntimes, 3), with columns [year, month, day].
    """

    with open(filename, "r") as f:
        lines = f.read().splitlines()

    current_line = nskip

    crop_areas = {
        "Data": [None] * ntimes,
        "YMD": np.zeros((ntimes, 3), dtype=int),
    }

    for ii in range(ntimes):
        parts = lines[current_line].split()

        # First token is assumed to be something like MM/DD/YYYY_...
        date_token = parts[0]
        date_part = date_token.split("_")[0]
        month, day, year = map(int, date_part.split("/"))

        if print_progress:
            print(date_token)

        crop_areas["YMD"][ii, :] = [year, month, day]

        # Remaining values on first line
        values = [float(x) for x in parts[1:] if _is_float(x)]

        arr = np.zeros((nelem, ncrops + 1), dtype=float)
        arr[0, :] = values

        current_line += 1

        for j in range(1, nelem):
            parts = lines[current_line].split()

            values = [float(x) for x in parts if _is_float(x)]
            arr[j, :] = values

            current_line += 1

        crop_areas["Data"][ii] = arr

    return crop_areas


def _is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def read_conductivity_prop(filename, n_nodes, n_layers, n_lines_skip):
    cols = ['ID', 'PKH', 'PS', 'PN', 'PV', 'PL']
    data = []

    with open(filename, 'r') as f:
        lines = f.readlines()[n_lines_skip:]  # skip header lines

    current_id = None
    for idx in range(n_nodes):
        start_idx = idx * n_layers
        block = lines[start_idx:start_idx + n_layers]
        if len(block) != n_layers:
            break

        # First line includes the ID
        parts = block[0].split()
        node_id = int(parts[0])
        values = [float(x) for x in parts[1:]]
        data.append([node_id] + values)
        # Remaining Nlayers-1 lines (continuations)
        for l in block[1:]:
            vals = [float(x) for x in l.split()]
            data.append([node_id] + vals)

    # Build main DataFrame
    df = pd.DataFrame(data, columns=cols)

    expected_rows = n_nodes * n_layers
    if len(df) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, but got {len(df)}")

    # Split into Nlayers DataFrames (one per layer)
    dfs = []
    for layer in range(n_layers):
        layer_df = df.iloc[layer::n_layers, :].reset_index(drop=True)
        dfs.append(layer_df)

    return dfs


def read_iwfm_smallwatersheds(filename):
    """
    Read an IWFM Small Watersheds file.

    Returns
    -------
    watersheds : list[dict]
        One dictionary per small watershed with keys:
        ID, AREAS, IWBTS, NWB, IWB, QMAXWB

        IWB and QMAXWB are lists with length NWB.
    """

    def is_comment_or_blank(line):
        s = line.strip()
        return (
            s == ""
            or s.startswith("C")
            or s.startswith("c")
            or s.startswith("#")
        )

    # Keep only non-comment, non-blank lines
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if not is_comment_or_blank(line)]

    # First two non-comment lines are ignored
    lines = lines[2:]

    # Next block of 4 non-comment lines: only first line is needed for NSW
    NSW = int(lines[0].split()[0])

    # Skip the full 4-line block
    lines = lines[4:]

    watersheds = []
    i = 0

    for _ in range(NSW):
        parts = lines[i].split()
        i += 1

        if len(parts) < 6:
            raise ValueError(f"Expected watershed header line with 6 values, got: {lines[i-1]}")

        ws_id = int(parts[0])
        areas = float(parts[1])
        iwbts = int(parts[2])
        nwb = int(parts[3])

        iwb = [int(parts[4])]
        qmaxwb = [float(parts[5])]

        # Remaining NWB-1 rows contain only IWB and QMAXWB
        for _ in range(nwb - 1):
            parts = lines[i].split()
            i += 1

            if len(parts) < 2:
                raise ValueError(f"Expected IWB/QMAXWB line with 2 values, got: {lines[i-1]}")

            iwb.append(int(parts[0]))
            qmaxwb.append(float(parts[1]))

        watersheds.append({
            "ID": ws_id,
            "AREAS": areas,
            "IWBTS": iwbts,
            "NWB": nwb,
            "IWB": iwb,
            "QMAXWB": qmaxwb,
        })

    return watersheds


def read_iwfm_stream_bed_param(filename, n_lines_skip = 4885, n_nodes=4634):
    """
    Read IWFM stream bed parameter file.

    Parameters
    ----------
    filename : str
        Input filename.

    n_lines_skip : int
        Number of lines to skip at the beginning of the file.

    n_nodes : int
        Number of stream nodes (rows) to read.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:

        IR      : River reach ID
        WETPR   : Wetted perimeter
        IRGW    : Groundwater node
        CSTRM   : Streambed conductance
        DSTRM   : Streambed thickness
    """

    records = []

    with open(filename, "r") as f:

        # Skip header lines
        for _ in range(n_lines_skip):
            next(f)

        while len(records) < n_nodes:

            try:
                line = next(f)
            except StopIteration:
                raise ValueError(
                    f"End of file reached after reading {len(records)} "
                    f"nodes; expected {n_nodes}."
                )

            # Remove inline comments beginning with '/'
            line = line.split("/", 1)[0].strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 5:
                continue

            records.append({
                "IR": int(parts[0]),
                "WETPR": float(parts[1]),
                "IRGW": int(parts[2]),
                "CSTRM": float(parts[3]),
                "DSTRM": float(parts[4]),
            })

    return pd.DataFrame(records)


def read_iwfm_stream_stage(filename, n_times=None, multiplier=1.0):
    """
    Read IWFM stream stage hydrograph output.

    Parameters
    ----------
    filename : str
        IWFM stream hydrograph output file.

    n_times : int or None, optional
        Number of time steps to read. If None, read all time steps.

    multiplier : float, optional
        Multiplier applied to stage values Hs.

    Returns
    -------
    dict
        Keys:
        HydId : ndarray
        IRV   : ndarray
        YMD   : ndarray with columns [year, month, day]
        Hs    : ndarray with shape (n_nodes, n_times)
    """

    hyd_id = None
    irv = None
    ymd_records = []
    hs_records = []

    with open(filename, "r") as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith("*"):
                # Header lines
                clean = s.lstrip("*").strip()

                if clean.upper().startswith("HYDROGRAPH ID"):
                    parts = clean.split()
                    hyd_id = np.array([int(x) for x in parts[2:]], dtype=int)

                elif clean.upper().startswith("NODES"):
                    parts = clean.split()
                    irv = np.array([int(x) for x in parts[1:]], dtype=int)

                continue

            # Data lines start here
            parts = s.split()

            if len(parts) < 2:
                continue

            date_token = parts[0].split("_")[0]
            month, day, year = [int(v) for v in date_token.split("/")]

            values = np.array([float(v) for v in parts[1:]], dtype=float)

            if irv is not None and len(values) != len(irv):
                raise ValueError(
                    f"Expected {len(irv)} stage values, "
                    f"but found {len(values)} on date {date_token}."
                )

            ymd_records.append([year, month, day])
            hs_records.append(values * multiplier)

            if n_times is not None and len(hs_records) >= n_times:
                break

    if hyd_id is None:
        raise ValueError("Could not find '* HYDROGRAPH ID' header line.")

    if irv is None:
        raise ValueError("Could not find '* NODES' header line.")

    ymd = np.array(ymd_records, dtype=int)

    # Stored as time x node, then transpose to node x time
    hs = np.vstack(hs_records).T

    return {
        "HydId": hyd_id,
        "IRV": irv,
        "YMD": ymd,
        "Hs": hs,
    }


