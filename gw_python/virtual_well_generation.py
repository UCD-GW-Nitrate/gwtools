from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
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
