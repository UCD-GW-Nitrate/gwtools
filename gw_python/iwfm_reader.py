import sys
import pandas as pd
import numpy as np
import re

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