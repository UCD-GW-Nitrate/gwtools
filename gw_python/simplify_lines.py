from shapely.geometry import LineString
import numpy as np

def douglas_peucker_3d(points, epsilon,original_idx=None):
    """
    Simplify a 3D polyline using the Ramer-Douglas-Peucker algorithm.

    Returns
    -------
    simplified_points : ndarray
        Simplified points with shape (M, 3).

    kept_rows : ndarray
        Row indices of the kept points in the original input array.
    """

    if points is None or len(points) == 0:
        return None, None

    points = np.asarray(points, dtype=float)

    if points.ndim != 2 or points.shape[1] < 2:
        return None, None

    if np.isnan(points).any():
        return None, None

    if points.shape[0] < 2:
        return None, None

    if original_idx is None:
        original_idx = np.arange(points.shape[0])
    else:
        original_idx = np.asarray(original_idx)

    start, end = points[0], points[-1]
    line_vec = end - start
    line_len = np.linalg.norm(line_vec)

    if line_len == 0:
        distances = np.linalg.norm(points - start, axis=1)
    else:
        line_unitvec = line_vec / line_len
        start_to_points = points - start
        proj_lengths = np.dot(start_to_points, line_unitvec)
        proj_points = start + np.outer(proj_lengths, line_unitvec)
        distances = np.linalg.norm(points - proj_points, axis=1)

    idx_max = np.argmax(distances)
    dmax = distances[idx_max]

    if dmax > epsilon:
        left_pts, left_idx = douglas_peucker_3d(
            points[:idx_max + 1],
            epsilon,
            original_idx[:idx_max + 1]
        )

        right_pts, right_idx = douglas_peucker_3d(
            points[idx_max:],
            epsilon,
            original_idx[idx_max:]
        )

        return (
            np.vstack((left_pts[:-1], right_pts)),
            np.concatenate((left_idx[:-1], right_idx))
        )

    else:
        return (
            np.vstack((start, end)),
            np.array([original_idx[0], original_idx[-1]])
        )


def simplify_to_linestring_3d(P, epsilon=10.0):
    """
    Simplify 3D points and convert them to a Shapely LineString.

    Returns
    -------
    line : shapely.geometry.LineString
        Simplified 3D LineString.

    kept_rows : ndarray
        Row indices of the kept points in the original input array.
    """

    if P is None or len(P) < 2:
        return None, None

    simplified, kept_rows = douglas_peucker_3d(np.asarray(P), epsilon)

    if simplified is None or simplified.shape[0] < 2:
        return None, None

    try:
        line = LineString(simplified)

        if not line.is_valid or line.is_empty:
            return None, None

        return line, kept_rows

    except Exception:
        return None, None