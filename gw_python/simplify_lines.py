from shapely.geometry import LineString
import numpy as np

def douglas_peucker_3d(points, epsilon):
    if points is None or len(points) == 0:
        return None

    points = np.asarray(points, dtype=float)

    # Remove NaNs or invalid points
    if np.isnan(points).any() or points.shape[1] != 3:
        return None

    if points.shape[0] < 2:
        return None

    # Line from first to last
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

    # Recursive subdivision
    if dmax > epsilon:
        left = douglas_peucker_3d(points[:idx_max + 1], epsilon)
        right = douglas_peucker_3d(points[idx_max:], epsilon)
        return np.vstack((left[:-1], right))
    else:
        return np.vstack((start, end))


def simplify_to_linestring_3d(P,epsilon = 10.0):
    if P is None or len(P) < 2:
        return None

    simplified = douglas_peucker_3d(np.array(P), epsilon)
    if simplified is None or simplified.shape[0] < 2:
        return None
    try:
        line = LineString(simplified)
        if not line.is_valid or line.is_empty:
            return None
        return line
    except Exception:
        return None