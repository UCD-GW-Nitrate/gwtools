import numpy as np
import pandas as pd
from shapely.geometry import LineString
import plotly.graph_objects as go


def read_npsat_trace(filename, simplify_tolerance=0.0):
    """
    Read NPSAT_TRACE streamline file.

    Parameters
    ----------
    filename : str
    simplify_tolerance : float
        Douglas-Peucker tolerance.
        0 -> no simplification.

    Returns
    -------
    summary_df : pandas.DataFrame

        Columns:
            Eid
            Sid
            xs ys zs
            xe ye ze
            vs
            age
            er

    streamlines : list[np.ndarray]

        List of (N,3) arrays containing simplified x,y,z coordinates.
    """

    summaries = []
    streamlines = []

    points = []
    exit_reason = -1

    current_eid = None
    current_sid = None

    def finalize():
        nonlocal points, exit_reason

        if len(points) == 0:
            return

        pts = np.asarray(points)

        xyz = pts[:, :3]
        vel = pts[:, 3]

        xs, ys, zs = xyz[0]
        xe, ye, ze = xyz[-1]
        vs = vel[0]

        # -------------------------
        # Compute travel age
        # -------------------------
        if len(xyz) == 1:
            age = 0.0

        else:
            seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)

            vavg = 0.5 * (vel[:-1] + vel[1:])

            # avoid divide by zero
            vavg[vavg <= 0] = np.nan

            age = np.nansum(seg / vavg)

        summaries.append(
            dict(
                Eid=current_eid,
                Sid=current_sid,
                xs=xs,
                ys=ys,
                zs=zs,
                xe=xe,
                ye=ye,
                ze=ze,
                vs=vs,
                age=age,
                er=exit_reason,
            )
        )

        # -------------------------
        # Douglas-Peucker simplify
        # -------------------------
        if simplify_tolerance > 0 and len(xyz) > 2:

            line = LineString(xyz)

            xyz = np.asarray(
                line.simplify(
                    simplify_tolerance,
                    preserve_topology=False,
                ).coords
            )

        streamlines.append(xyz)

        points = []
        exit_reason = -1

    with open(filename) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            vals = line.split()

            pid = int(vals[0])

            # ----------------------------------------------------------
            # termination row
            # ----------------------------------------------------------
            if pid < 0:

                # exit reason = 3rd from the end
                exit_reason = int(vals[-3])

                finalize()

                current_eid = None
                current_sid = None

                continue

            eid = int(vals[1])
            sid = int(vals[2])

            # ----------------------------------------------------------
            # Missing termination row
            # ----------------------------------------------------------
            if (
                current_eid is not None
                and (eid != current_eid or sid != current_sid)
            ):
                finalize()

            current_eid = eid
            current_sid = sid

            x = float(vals[3])
            y = float(vals[4])
            z = float(vals[5])
            v = float(vals[6])

            points.append((x, y, z, v))

    # file ended without termination
    finalize()

    summary_df = pd.DataFrame(summaries)

    return summary_df, streamlines

def plot_streamlines_3d(
    streamlines,
    fig=None,
    color="royalblue",
    width=2,
    opacity=1.0,
    show_start=False,
    show_end=False,
):
    """
    Plot streamlines as 3D lines.

    Parameters
    ----------
    streamlines : list[np.ndarray]
        List of (N,3) arrays containing x,y,z coordinates.

    fig : go.Figure or None
        Existing figure to add traces to. If None a new figure is created.

    color : str
        Line color.

    width : float
        Line width.

    opacity : float
        Line opacity.

    show_start : bool
        Plot starting points.

    show_end : bool
        Plot ending points.

    Returns
    -------
    fig : go.Figure
    """

    if fig is None:
        fig = go.Figure()

    start_x = []
    start_y = []
    start_z = []

    end_x = []
    end_y = []
    end_z = []

    for xyz in streamlines:

        if len(xyz) == 0:
            continue

        fig.add_trace(
            go.Scatter3d(
                x=xyz[:, 0],
                y=xyz[:, 1],
                z=xyz[:, 2],
                mode="lines",
                line=dict(
                    color=color,
                    width=width,
                ),
                opacity=opacity,
                showlegend=False,
                hoverinfo="skip",
            )
        )

        if show_start:
            start_x.append(xyz[0, 0])
            start_y.append(xyz[0, 1])
            start_z.append(xyz[0, 2])

        if show_end:
            end_x.append(xyz[-1, 0])
            end_y.append(xyz[-1, 1])
            end_z.append(xyz[-1, 2])

    if show_start and start_x:
        fig.add_trace(
            go.Scatter3d(
                x=start_x,
                y=start_y,
                z=start_z,
                mode="markers",
                marker=dict(
                    size=4,
                    color="green",
                ),
                name="Start",
            )
        )

    if show_end and end_x:
        fig.add_trace(
            go.Scatter3d(
                x=end_x,
                y=end_y,
                z=end_z,
                mode="markers",
                marker=dict(
                    size=4,
                    color="red",
                ),
                name="End",
            )
        )

    fig.update_layout(
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, b=0, t=30),
    )

    return fig