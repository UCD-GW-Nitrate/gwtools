import pandas as pd

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
