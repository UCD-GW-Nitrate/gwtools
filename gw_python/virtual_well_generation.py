from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


@dataclass
class PairKDE:
    x_field: str
    y_field: str
    kde: gaussian_kde
    log10_min: np.ndarray
    log10_max: np.ndarray
    n_used: int
    xy_log: np.ndarray
    xy_norm: np.ndarray
    max_density: float

    def transform(self, x, y):
        """Convert physical x,y values to normalized log10 KDE coordinates."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        lx = np.log10(x)
        ly = np.log10(y)

        vals = np.vstack([lx, ly]).T
        vals_norm = (vals - self.log10_min) / (self.log10_max - self.log10_min)

        return vals_norm.T

    def density(self, x, y):
        """Return KDE density score at physical x,y values."""
        pts = self.transform(x, y)
        return self.kde(pts)

    def accept_probability(self, x, y):
        """Return normalized acceptance probability in [0, 1]."""
        d = self.density(x, y)
        return np.minimum(d / self.max_density, 1.0)

    def sample_candidate(self):
        """Generate one random candidate pair uniformly in log10 space."""
        lx = np.random.uniform(self.log10_min[0], self.log10_max[0])
        ly = np.random.uniform(self.log10_min[1], self.log10_max[1])

        return 10 ** lx, 10 ** ly

    def sample(self, n, max_iter=1_000_000):
        """Generate n accepted random pairs using rejection sampling."""
        accepted = []

        for _ in range(max_iter):
            x, y = self.sample_candidate()
            p = self.accept_probability(x, y)[0]

            if np.random.rand() < p:
                accepted.append((x, y))

                if len(accepted) == n:
                    return np.array(accepted)

        raise RuntimeError(
            f"Only accepted {len(accepted)} points after {max_iter} iterations."
        )

    def plot(
            self,
            ax=None,
            cmap="jet",
            show_points=True,
            nx=200,
            ny=200,
            levels=20,
            filled=True,
            point_size=10,
            point_alpha=0.35,
    ):
        """
        Plot KDE heatmap in physical space with log-scaled axes.
        KDE is evaluated in normalized log10 space.
        """

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))

        # Grid in log10 physical coordinates
        lx = np.linspace(self.log10_min[0], self.log10_max[0], nx)
        ly = np.linspace(self.log10_min[1], self.log10_max[1], ny)

        LX, LY = np.meshgrid(lx, ly)

        # Convert log10 grid to normalized KDE coordinates
        XN = (LX - self.log10_min[0]) / (
                self.log10_max[0] - self.log10_min[0]
        )
        YN = (LY - self.log10_min[1]) / (
                self.log10_max[1] - self.log10_min[1]
        )

        pts = np.vstack([XN.ravel(), YN.ravel()])
        ZZ = self.kde(pts).reshape(ny, nx)

        # Physical coordinates for plotting
        X = 10 ** LX
        Y = 10 ** LY

        if filled:
            h = ax.contourf(X, Y, ZZ, levels=levels, cmap=cmap)
        else:
            h = ax.contour(X, Y, ZZ, levels=levels, cmap=cmap)

        plt.colorbar(h, ax=ax, label="KDE density")

        if show_points:
            ax.scatter(
                10 ** self.xy_log[:, 0],
                10 ** self.xy_log[:, 1],
                s=point_size,
                c="k",
                alpha=point_alpha,
                linewidths=0,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(self.x_field)
        ax.set_ylabel(self.y_field)

        ax.set_title(
            f"KDE association: {self.x_field} vs {self.y_field}\n"
            f"N = {self.n_used}"
        )

        return ax


def fit_pair_kde(df, x_field, y_field, bandwidth=None):
    """
    Fit a 2D KDE on log10-transformed and normalized x,y pairs.

    Cleaning:
    - drops NaN
    - drops inf
    - drops zero and negative values
    """

    tmp = df[[x_field, y_field]].copy()
    tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna()

    # Keep only physically meaningful positive values
    tmp = tmp[(tmp[x_field] > 0) & (tmp[y_field] > 0)]

    if len(tmp) < 5:
        raise ValueError("Not enough valid positive x,y pairs to fit KDE.")

    # Log10 transform
    log_vals = np.log10(tmp[[x_field, y_field]].to_numpy(dtype=float))

    log10_min = log_vals.min(axis=0)
    log10_max = log_vals.max(axis=0)

    if np.any(log10_max == log10_min):
        raise ValueError("One of the fields has no variation after log10 transform.")

    # Normalize to [0, 1]
    log_vals_norm = (log_vals - log10_min) / (log10_max - log10_min)

    # scipy expects shape: ndim x nsamples
    xy = log_vals_norm.T

    kde = gaussian_kde(xy, bw_method=bandwidth)
    max_density = kde(xy).max()

    model = PairKDE(
        x_field=x_field,
        y_field=y_field,
        kde=kde,
        log10_min=log10_min,
        log10_max=log10_max,
        n_used=len(tmp),
        xy_log=log_vals,
        xy_norm=log_vals_norm,
        max_density=max_density,
    )

    # Store max density using training points for acceptance normalization
    model.max_density = kde(xy).max()

    return model


def make_raster_probability_interpolator(
    prob_grid,
    transform,
    min_prob=0.10,
    fallback_prob=0.10,
    method="linear",
):
    """
    Build an interpolator for a raster probability grid.

    Returns a callable with signature f(x, y). The returned value is clipped to
    [0, 1], uses fallback_prob outside the raster or over NaN cells, and applies
    min_prob as a floor for finite interpolated values.
    """

    prob_grid = np.asarray(prob_grid, dtype=float)

    if not np.isclose(transform.b, 0.0) or not np.isclose(transform.d, 0.0):
        raise ValueError("Rotated/sheared raster transforms are not supported.")

    n_rows, n_cols = prob_grid.shape
    cols = np.arange(n_cols)
    rows = np.arange(n_rows)

    xs = transform.c + (cols + 0.5) * transform.a
    ys = transform.f + (rows + 0.5) * transform.e

    grid = prob_grid.copy()
    if ys[0] > ys[-1]:
        ys = ys[::-1]
        grid = grid[::-1, :]
    if xs[0] > xs[-1]:
        xs = xs[::-1]
        grid = grid[:, ::-1]

    interpolator = RegularGridInterpolator(
        (ys, xs),
        grid,
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )

    def interpolate(x, y):
        p = float(interpolator((y, x)))
        if not np.isfinite(p):
            p = fallback_prob
        return float(np.clip(max(p, min_prob), 0.0, 1.0))

    return interpolate

def _is_valid_positive(v):
    return pd.notna(v) and np.isfinite(v) and float(v) > 0

def _model_values_for_field(model, field):
    if model.x_field == field:
        return 10 ** model.xy_log[:, 0]
    elif model.y_field == field:
        return 10 ** model.xy_log[:, 1]
    else:
        return None

def get_field_range(field, models):
    """
    Get the narrowest available physical-value range for a field
    from all PairKDE models that contain that field.
    """
    ranges = []

    for model in models:
        vals = _model_values_for_field(model, field)
        if vals is not None:
            ranges.append((np.nanmin(vals), np.nanmax(vals)))

    if len(ranges) == 0:
        raise ValueError(f"No KDE model contains field '{field}'.")

    # Try intersection first
    lo = max(r[0] for r in ranges)
    hi = min(r[1] for r in ranges)

    if lo < hi:
        return lo, hi

    # Fallback: use the narrowest individual range
    widths = [r[1] - r[0] for r in ranges]
    return ranges[int(np.argmin(widths))]

def sample_field_from_ecdf(field, models, rng=None):
    """
    Sample a field from the empirical observed values.

    Important:
    - The admissible range is the narrowest shared range from the KDE models.
    - Duplicate values are retained.
    - Therefore commonly observed values are sampled more frequently.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Narrowest admissible range from all models containing this field
    lo, hi = get_field_range(field, models)

    vals_all = []

    for model in models:
        vals = _model_values_for_field(model, field)
        if vals is None:
            continue
        vals = np.asarray(vals, dtype=float)
        vals = vals[
            np.isfinite(vals) &
            (vals > 0) &
            (vals >= lo) &
            (vals <= hi)
        ]

        if len(vals) > 0:
            vals_all.append(vals)

    if len(vals_all) == 0:
        # Fallback only if no empirical values are available
        return float(10 ** rng.uniform(np.log10(lo), np.log10(hi)))

    vals_all = np.concatenate(vals_all)

    if len(vals_all) == 0:
        # fallback to log-uniform
        return 10 ** rng.uniform(np.log10(lo), np.log10(hi))

    # Keep duplicates intentionally.
    # Common values therefore have higher sampling probability.
    return float(rng.choice(vals_all))

def sample_field_loguniform(field, models, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    lo, hi = get_field_range(field, models)
    return float(10 ** rng.uniform(np.log10(lo), np.log10(hi)))

def sample_Q_conditioned(
    q_field,
    models,
    conditioning_prob=None,
    rng=None,
    method="ecdf",
    low_quantile=0.05,
    high_quantile=0.95,
    spread_frac=0.15,
):
    """
    Sample Q from observed WCR values, optionally conditioned by an external
    probability/intensity value.

    conditioning_prob = None:
        sample directly from the empirical Q distribution.

    conditioning_prob = 0:
        sample from the low-Q part of the empirical distribution.

    conditioning_prob = 1:
        sample from the high-Q part of the empirical distribution.
    """
    if rng is None:
        rng = np.random.default_rng()

    vals_all = []

    for model in models:
        vals = _model_values_for_field(model, q_field)

        if vals is not None:
            vals = np.asarray(vals, dtype=float)
            vals = vals[np.isfinite(vals) & (vals > 0)]

            if len(vals) > 0:
                vals_all.append(vals)

    if len(vals_all) == 0:
        raise ValueError(f"No valid Q values found for field '{q_field}'.")

    q_vals = np.concatenate(vals_all)

    if len(q_vals) == 0:
        raise ValueError(f"No valid Q values found for field '{q_field}'.")

    # Case 1: no conditioning. Use empirical distribution directly.
    if conditioning_prob is None:
        if method == "ecdf":
            return float(rng.choice(q_vals))

        q_min = np.min(q_vals)
        q_max = np.max(q_vals)

        return float(
            10 ** rng.uniform(
                np.log10(q_min),
                np.log10(q_max)
            )
        )

    # Case 2: external conditioning probability/intensity is provided.
    p = float(np.clip(conditioning_prob, 0.0, 1.0))

    qlo = np.quantile(q_vals, low_quantile)
    qhi = np.quantile(q_vals, high_quantile)

    q_target = qlo + p * (qhi - qlo)

    spread = spread_frac * (qhi - qlo)

    q_min = max(np.min(q_vals), q_target - spread)
    q_max = min(np.max(q_vals), q_target + spread)

    q_candidates = q_vals[
        (q_vals >= q_min) &
        (q_vals <= q_max)
    ]

    if method == "ecdf" and len(q_candidates) > 0:
        return float(rng.choice(q_candidates))

    # Fallback: log-uniform within the conditioned interval
    if q_min <= 0 or q_max <= 0 or q_min >= q_max:
        return float(rng.choice(q_vals))

    return float(
        10 ** rng.uniform(
            np.log10(q_min),
            np.log10(q_max)
        )
    )


def impute_D_SL_Q(
    D,
    SL,
    Q,
    model_D_Q,
    model_D_SL,
    model_SL_Q,
    swat_irrg_prob=0.0,
    max_iter=10000,
    rng=None,
    sample_method="ecdf",
    D_fld="TOTALCOMPLETEDDEPTH",
    SL_fld="ScreenLength",
    Q_fld="WELLYIELD",
):
    """
    Impute missing D, SL, and Q using three pairwise KDE models.

    Observed positive values are kept fixed.
    Missing or invalid values are sampled.

    Parameters
    ----------
    D : float
        TOTALCOMPLETEDDEPTH.
    SL : float
        ScreenLength.
    Q : float
        WELLYIELD.
    model_D_Q, model_D_SL, model_SL_Q : PairKDE
        Pairwise KDE models.
    swat_irrg_prob : float
        Raster-derived probability/intensity value in [0, 1].
    max_iter : int
        Maximum rejection-sampling attempts.
    rng : np.random.Generator or None
        Random generator.
    sample_method : {"ecdf", "loguniform"}
        How to sample non-Q missing variables.

    Returns
    -------
    D_imp, SL_imp, Q_imp, accepted, n_iter
    """

    if rng is None:
        rng = np.random.default_rng()

    models = [model_D_Q, model_D_SL, model_SL_Q]

    D_fixed = _is_valid_positive(D)
    SL_fixed = _is_valid_positive(SL)
    Q_fixed = _is_valid_positive(Q)

    D0 = float(D) if D_fixed else np.nan
    SL0 = float(SL) if SL_fixed else np.nan
    Q0 = float(Q) if Q_fixed else np.nan

    if D_fixed and SL_fixed and Q_fixed:
        return D0, SL0, Q0, True, 0

    sampler = sample_field_from_ecdf if sample_method == "ecdf" else sample_field_loguniform

    for it in range(1, max_iter + 1):

        D_cand = D0 if D_fixed else sampler(D_fld, models, rng)
        SL_cand = SL0 if SL_fixed else sampler(SL_fld, models, rng)

        Q_cand = (
            Q0
            if Q_fixed
            else sample_Q_conditioned(
                Q_fld,
                models,
                conditioning_prob=swat_irrg_prob,
                rng=rng,
                method=sample_method,
            )
        )

        # Physical constraint: screen length cannot exceed completed depth
        if SL_cand > D_cand:
            continue

        p_D_Q = model_D_Q.accept_probability(D_cand, Q_cand)[0]
        p_D_SL = model_D_SL.accept_probability(D_cand, SL_cand)[0]
        p_SL_Q = model_SL_Q.accept_probability(SL_cand, Q_cand)[0]

        if (
            rng.random() < p_D_Q and
            rng.random() < p_D_SL and
            rng.random() < p_SL_Q
        ):
            return D_cand, SL_cand, Q_cand, True, it

    return D0, SL0, Q0, False, max_iter

def _positive_values(df, field):
    vals = pd.to_numeric(df[field], errors="coerce").to_numpy(dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    return vals


def _basic_dist_stats(vals, prefix):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]

    out = {
        f"{prefix}_n": len(vals),
        f"{prefix}_mean": np.nan,
        f"{prefix}_std": np.nan,
        f"{prefix}_min": np.nan,
        f"{prefix}_p05": np.nan,
        f"{prefix}_p25": np.nan,
        f"{prefix}_p50": np.nan,
        f"{prefix}_p75": np.nan,
        f"{prefix}_p95": np.nan,
        f"{prefix}_max": np.nan,
        f"{prefix}_ecdf_x": np.array([]),
        f"{prefix}_ecdf_y": np.array([]),
        f"{prefix}_pdf_bins": np.array([]),
        f"{prefix}_pdf_density": np.array([]),
    }

    if len(vals) == 0:
        return out

    vals_sort = np.sort(vals)

    out[f"{prefix}_mean"] = np.mean(vals)
    out[f"{prefix}_std"] = np.std(vals)
    out[f"{prefix}_min"] = np.min(vals)
    out[f"{prefix}_p05"] = np.percentile(vals, 5)
    out[f"{prefix}_p25"] = np.percentile(vals, 25)
    out[f"{prefix}_p50"] = np.percentile(vals, 50)
    out[f"{prefix}_p75"] = np.percentile(vals, 75)
    out[f"{prefix}_p95"] = np.percentile(vals, 95)
    out[f"{prefix}_max"] = np.max(vals)

    out[f"{prefix}_ecdf_x"] = vals_sort
    out[f"{prefix}_ecdf_y"] = np.arange(1, len(vals_sort) + 1) / len(vals_sort)

    pdf_density, pdf_edges = np.histogram(vals, bins=30, density=True)
    out[f"{prefix}_pdf_bins"] = pdf_edges
    out[f"{prefix}_pdf_density"] = pdf_density

    return out


def make_region_stats(ireg, group_name, wcr_obs, completed_df,
    D_fld = "TOTALCOMPLETEDDEPTH",
    SL_fld = "ScreenLength",
    Q_fld = "WELLYIELD",
):
    rec = {
        "SubRegion": ireg + 1,
        "Group": group_name,
        "n_wells": len(wcr_obs),
        "n_D_obs": len(_positive_values(wcr_obs, D_fld)),
        "n_SL_obs": len(_positive_values(wcr_obs, SL_fld)),
        "n_Q_obs": len(_positive_values(wcr_obs, Q_fld)),
        "n_D_Q_pairs_obs": len(
            wcr_obs[
                (pd.to_numeric(wcr_obs[D_fld], errors="coerce") > 0) &
                (pd.to_numeric(wcr_obs[Q_fld], errors="coerce") > 0)
            ]
        ),
        "n_D_SL_pairs_obs": len(
            wcr_obs[
                (pd.to_numeric(wcr_obs[D_fld], errors="coerce") > 0) &
                (pd.to_numeric(wcr_obs[SL_fld], errors="coerce") > 0)
            ]
        ),
        "n_SL_Q_pairs_obs": len(
            wcr_obs[
                (pd.to_numeric(wcr_obs[SL_fld], errors="coerce") > 0) &
                (pd.to_numeric(wcr_obs[Q_fld], errors="coerce") > 0)
            ]
        ),
        "n_accepted": int(completed_df["Accepted"].sum()) if len(completed_df) > 0 else 0,
        "mean_n_iter": completed_df["NIter"].mean() if len(completed_df) > 0 else np.nan,
    }

    field_map = {
        "D": (D_fld, "Depth"),
        "SL": (SL_fld, "SL"),
        "Q": (Q_fld, "Wyield"),
    }

    for short_name, (obs_field, comp_field) in field_map.items():
        obs_vals = _positive_values(wcr_obs, obs_field)
        comp_vals = pd.to_numeric(completed_df[comp_field], errors="coerce").to_numpy(dtype=float)
        comp_vals = comp_vals[np.isfinite(comp_vals) & (comp_vals > 0)]

        rec.update(_basic_dist_stats(obs_vals, f"{short_name}_obs"))
        rec.update(_basic_dist_stats(comp_vals, f"{short_name}_completed"))

    return rec
