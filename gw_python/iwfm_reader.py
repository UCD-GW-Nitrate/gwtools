import pandas as pd

def read_nodes(filename):
    nd = None
    fact = None
    node_data = []  # List to store coordinate tuples

    # 1. Read all lines
    try:
        with open(filename, "r") as f:
            all_lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None

    line_index = 0
    # 2. Extract ND and FACT
    while line_index < len(all_lines):
        line = all_lines[line_index].strip()
        line_index += 1

        # Skip comment lines ('C') and empty lines
        if line.startswith('C') or not line:
            continue

        # Found the ND line
        try:
            nd = int(line.split()[0])
            line = all_lines[line_index].strip()
            line_index += 1
            fact = float(line.split()[0])
            break
        except (ValueError, IndexError):
            # This handles separator lines or unexpected text
            print(f"Error: WHILE READING ['{line}' not found.")
            return None

    # 3. Preallocate DataFrame
    nodes_df = pd.DataFrame({
        'NodeID': [0] * nd,
        'X': [0.0] * nd,
        'Y': [0.0] * nd
    })

    node_count = 0
    while line_index < len(all_lines) and node_count < nd:
        line = all_lines[line_index].strip()
        line_index += 1

        # Skip comment lines ('C') and empty lines
        if line.startswith('C') or not line:
            continue

        try:
            parts = line.split()
            if len(parts) >= 3:
                node_id = int(parts[0])
                coord_x = float(parts[1])
                coord_y = float(parts[2])

                nodes_df.loc[node_count] = [node_id, coord_x, coord_y]
                node_count += 1
        except (ValueError, IndexError):
            # Skip header/separator lines between FACT and coordinates
            continue

    if node_count != nd:
        print(f"Warning: Expected {nd} nodes, but read only {node_count} nodes.")

    return nodes_df