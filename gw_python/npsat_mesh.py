import gmsh
import numpy as np
import plotly.graph_objects as go

def generate_polygon_quad_mesh(polygon_coords, target_mesh_size=0.5):
    """
    Generates a 2D quad mesh for a polygon. Safe to run multiple times in Jupyter.

    Parameters:
    - polygon_coords: List/array of (x, y) coordinates defining the polygon loop.
    - target_mesh_size: Element sizing factor.

    Returns:
    - coordinates: NumPy array of shape (N, 2) tracking X and Y coordinates.
    - node_ids: NumPy array of shape (N,) matching the rows of the coordinates.
    """
    # 1. Jupyter Safe Session Management
    if not gmsh.isInitialized():
        gmsh.initialize()
    else:
        gmsh.clear()  # Wipes previous model data so you can run it cleanly again

    gmsh.model.add("Jupyter_Quad_Mesh")
    geo = gmsh.model.geo

    try:
        # 2. Build the Geometry
        point_tags = [geo.addPoint(x, y, 0.0, target_mesh_size) for x, y in polygon_coords]

        line_tags = []
        num_points = len(point_tags)
        for i in range(num_points):
            line_tags.append(geo.addLine(point_tags[i], point_tags[(i + 1) % num_points]))

        curve_loop = geo.addCurveLoop(line_tags)
        surface = geo.addPlaneSurface([curve_loop])
        geo.synchronize()

        # 3. Apply Recombination Options (Triangle -> Quad)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # Blossom algorithm
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for Quads

        # Turn off console print spam inside your Jupyter scrolling outputs
        gmsh.option.setNumber("General.Terminal", 0)

        # 4. Mesh and Extract
        gmsh.model.mesh.generate(2)

        node_ids, coords_flat, _ = gmsh.model.mesh.getNodes(dim=-1, tag=-1)

        # 5. Format to exact matching NumPy rows
        coordinates = np.reshape(coords_flat, (-1, 3))[:, :2]  # Reshape and drop Z coordinate
        node_ids = np.array(node_ids, dtype=np.int64)

        return coordinates, node_ids

    except Exception as e:
        print(f"An error occurred during mesh generation: {e}")
        return None, None



def plot_quad_mesh(polygon_coords=None, coords=None, ids=None):
    """
    Plot a polygon boundary and/or quadrilateral mesh.

    If polygon_coords is None, only the mesh is plotted.
    If coords and ids are None, only the polygon boundary is plotted.
    """

    fig = go.Figure()
    has_mesh = False
    has_polygon = False

    # --------------------------------------------------
    # 1. Plot mesh if available
    # --------------------------------------------------
    if coords is not None and ids is not None and len(coords) > 0 and len(ids) > 0:
        coords = np.asarray(coords, dtype=float)

        try:
            elem_types, _, element_node_tags = gmsh.model.mesh.getElements(dim=2)
            quad_idx = np.where(elem_types == 3)[0]

            if len(quad_idx) > 0:
                quad_node_tags = element_node_tags[quad_idx[0]].reshape(-1, 4)
                id_to_idx = {node_id: idx for idx, node_id in enumerate(ids)}

                mesh_x = []
                mesh_y = []

                for quad in quad_node_tags:
                    try:
                        idxs = [id_to_idx[nid] for nid in quad]
                        idxs.append(idxs[0])

                        mesh_x.extend([coords[idx, 0] for idx in idxs] + [None])
                        mesh_y.extend([coords[idx, 1] for idx in idxs] + [None])
                    except KeyError:
                        continue

                fig.add_trace(go.Scatter(
                    x=mesh_x,
                    y=mesh_y,
                    mode="lines",
                    line=dict(color="rgb(120, 150, 180)", width=1),
                    name="Quad Mesh",
                    hoverinfo="skip"
                ))

                fig.add_trace(go.Scatter(
                    x=coords[:, 0],
                    y=coords[:, 1],
                    mode="markers",
                    marker=dict(size=4, color="rgb(60, 90, 120)"),
                    name="Mesh Nodes",
                    text=[f"Node ID: {nid}" for nid in ids],
                    hoverinfo="text"
                ))

                has_mesh = True

        except Exception as e:
            print(f"Mesh data supplied but active Gmsh cache could not read structure: {e}")

    # --------------------------------------------------
    # 2. Plot polygon boundary only if provided
    # --------------------------------------------------
    if polygon_coords is not None:
        poly_arr = np.asarray(polygon_coords, dtype=float)

        if poly_arr.ndim != 2 or poly_arr.shape[0] < 3 or poly_arr.shape[1] < 2:
            print("Invalid polygon coordinates provided. Expected at least 3 points with x/y coordinates.")
        else:
            poly_arr = poly_arr[:, :2]

            poly_x = np.append(poly_arr[:, 0], poly_arr[0, 0])
            poly_y = np.append(poly_arr[:, 1], poly_arr[0, 1])

            fig.add_trace(go.Scatter(
                x=poly_x,
                y=poly_y,
                mode="lines",
                line=dict(color="firebrick", width=4),
                name="Polygon Boundary"
            ))

            has_polygon = True

    # --------------------------------------------------
    # 3. Title
    # --------------------------------------------------
    if has_mesh and has_polygon:
        title_text = "Polygon & Quadrilateral Mesh"
    elif has_mesh:
        title_text = "Quadrilateral Mesh"
    elif has_polygon:
        title_text = "Polygon Boundary"
    else:
        title_text = "No Valid Polygon or Mesh Data"

    fig.update_layout(
        title=title_text,
        xaxis_title="X",
        yaxis_title="Y",
        template="plotly_white",
        showlegend=True,
        yaxis=dict(scaleanchor="x", scaleratio=1),
        width=700,
        height=600
    )

    fig.show()


def export_quad_mesh_to_file(coords, ids, filename="mesh_output.txt", cw=False):
    """
    Writes mesh coordinates and quadrilateral connectivity to a file.

    Format:
    Ncoords, Nquad
    Coords [Ncoords x 2]
    Quads [Nquad x 4] (0-based indices mapping directly to the rows of Coords)

    Parameters:
    - coords: NumPy array of shape (N, 2) from Gmsh.
    - ids: NumPy array of shape (N,) mapping node tags to rows.
    - filename: Output destination path.
    - cw: True for clockwise node ordering, False for counter-clockwise (CCW).
    """
    if coords is None or ids is None or len(coords) == 0 or len(ids) == 0:
        print("Missing mesh data. Aborting file write.")
        return

    try:
        # 1. Fetch 2D elements from active Gmsh session memory
        elem_types, _, element_node_tags = gmsh.model.mesh.getElements(dim=2)
        quad_idx = np.where(elem_types == 3)[0]

        if len(quad_idx) == 0:
            print("No quadrilateral elements found in the current Gmsh session.")
            return

        # Reshape into a structural matrix of Gmsh Node Tags (M x 4)
        quad_node_tags = element_node_tags[quad_idx[0]].reshape(-1, 4)

        # 2. Map Gmsh Node Tags to 0-based row indexes matching 'coords'
        id_to_zero_index = {node_id: idx for idx, node_id in enumerate(ids)}

        quads_zero_indexed = []
        reoriented_count = 0
        for quad in quad_node_tags:
            # Translate tags into row numbers [0, 1, 2...]
            indices = [id_to_zero_index[nid] for nid in quad]

            # Extract the literal physical coordinates for orientation parsing
            p = coords[indices]  # Shape (4, 2)

            # 3. Calculate signed area via Shoelace formula to determine orientation
            # Area = 0.5 * sum(x_i * y_{i+1} - x_{i+1} * y_i)
            signed_area = 0.0
            for i in range(4):
                x1, y1 = p[i]
                x2, y2 = p[(i + 1) % 4]
                signed_area += (x1 * y2 - x2 * y1)

            is_currently_ccw = signed_area > 0

            # Modify orientation if it doesn't match the user's intent
            if cw and is_currently_ccw:
                # CCW -> CW: Swap the ring flow order (Keep start, reverse the rest)
                indices = [indices[0], indices[3], indices[2], indices[1]]
                reoriented_count += 1
            elif not cw and not is_currently_ccw:
                # CW -> CCW: Swap the ring flow order
                indices = [indices[0], indices[3], indices[2], indices[1]]
                reoriented_count += 1

            quads_zero_indexed.append(indices)

        quads_zero_indexed = np.array(quads_zero_indexed, dtype=np.int64)

        # 4. Write to File
        n_coords = coords.shape[0]
        n_quads = quads_zero_indexed.shape[0]

        with open(filename, 'w') as f:
            # Header line: Ncoords, Nquad
            f.write(f"{n_coords} {n_quads}\n")

            # Block 1: Coords [Ncoords x 2]
            for row in coords:
                f.write(f"{row[0]:.8f} {row[1]:.8f}\n")

            # Block 2: Quads [Nquad x 4]
            for quad in quads_zero_indexed:
                f.write(f"{quad[0]} {quad[1]} {quad[2]} {quad[3]}\n")

        print(f"Successfully exported mesh layout to {filename} ({'CW' if cw else 'CCW'} ordering)")
        print(f"{reoriented_count} out of {n_quads} elements were reoriented.")

    except Exception as e:
        print(f"An error occurred while attempting file export: {e}")


def plot_mesh(polygon_coords=None, coords=None, ids=None,polyline_width = 4.0):
    """
    Plot a polygon boundary and/or mesh.

    Parameters
    ----------
    polygon_coords : array-like, optional
        Polygon boundary coordinates, shape (N, 2) or (N, 3).
        If None, no polygon outline is plotted.

    coords : array-like, optional
        Mesh node coordinates, shape (Nnodes, 2) or (Nnodes, 3).

    ids : ndarray, optional
        Mesh connectivity array, shape (Nelements, ngon).
        Indices are assumed to be 0-based.
        Can represent triangles, quads, pentagons, etc.
    """

    fig = go.Figure()
    has_mesh = False
    has_polygon = False

    # --------------------------------------------------
    # Plot mesh
    # --------------------------------------------------
    if coords is not None and ids is not None:
        coords = np.asarray(coords, dtype=float)
        ids = np.asarray(ids, dtype=np.int64)

        if coords.ndim != 2 or coords.shape[1] < 2:
            print("Invalid coords. Expected shape (Nnodes, 2) or (Nnodes, 3).")
        elif ids.ndim != 2 or ids.shape[1] < 3:
            print("Invalid ids. Expected shape (Nelements, ngon), with ngon >= 3.")
        else:
            coords_xy = coords[:, :2]

            mesh_x = []
            mesh_y = []

            for elem in ids:
                elem_closed = np.append(elem, elem[0])

                mesh_x.extend(coords_xy[elem_closed, 0].tolist() + [None])
                mesh_y.extend(coords_xy[elem_closed, 1].tolist() + [None])

            fig.add_trace(go.Scatter(
                x=mesh_x,
                y=mesh_y,
                mode="lines",
                line=dict(color="rgb(120, 150, 180)", width=1),
                name=f"{ids.shape[1]}-node mesh",
                hoverinfo="skip"
            ))

            fig.add_trace(go.Scatter(
                x=coords_xy[:, 0],
                y=coords_xy[:, 1],
                mode="markers",
                marker=dict(size=4, color="rgb(60, 90, 120)"),
                name="Mesh nodes",
                text=[f"Node index: {i}" for i in range(coords_xy.shape[0])],
                hoverinfo="text"
            ))

            has_mesh = True

    # --------------------------------------------------
    # Plot polygon boundary
    # --------------------------------------------------
    if polygon_coords is not None:
        poly_arr = np.asarray(polygon_coords, dtype=float)

        if poly_arr.ndim != 2 or poly_arr.shape[0] < 3 or poly_arr.shape[1] < 2:
            print("Invalid polygon_coords. Expected at least 3 points with x/y coordinates.")
        else:
            poly_arr = poly_arr[:, :2]

            poly_x = np.append(poly_arr[:, 0], poly_arr[0, 0])
            poly_y = np.append(poly_arr[:, 1], poly_arr[0, 1])

            fig.add_trace(go.Scatter(
                x=poly_x,
                y=poly_y,
                mode="lines",
                line=dict(color="firebrick", width=polyline_width),
                name="Polygon boundary"
            ))

            has_polygon = True

    if has_mesh and has_polygon:
        title_text = "Polygon & Mesh"
    elif has_mesh:
        title_text = "Mesh"
    elif has_polygon:
        title_text = "Polygon Boundary"
    else:
        title_text = "No Valid Polygon or Mesh Data"

    fig.update_layout(
        title=title_text,
        xaxis_title="X",
        yaxis_title="Y",
        template="plotly_white",
        showlegend=True,
        yaxis=dict(scaleanchor="x", scaleratio=1),
        width=700,
        height=600
    )

    fig.show()
