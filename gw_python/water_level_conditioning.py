import numpy as np
import pandas as pd
import skgstat as skg


def condition_model_to_measurements(
    meas_df,
    model_df,
    n_data=6,
    slope_threshold=2.0,
    kriging_params=None,
    meas_dict=None,
    model_dict=None,
    verbose=True,
):
    """
    Condition model heads to measured heads using:

        WL_conditioned =
            kriged_measured_field
            +
            (model_field - reconstructed_model_field_from_measurement_locations)

    Parameters
    ----------
    meas_df : pandas.DataFrame
        Measurement dataframe.

    model_df : pandas.DataFrame
        Model dataframe.

    n_data : int
        Minimum number of measurements required per measured point.

    slope_threshold : float
        Maximum absolute allowed slope. Assumes meas_df slope field has same units
        intended by the user, typically ft/year.

    kriging_params : dict
        Dictionary with optional keys:
            {
                "model": "spherical",
                "n_lags": 15,
                "maxlag_meas": 400000,
                "maxlag_model": "median",
                "maxlag_sim_on_meas": 400000,
                "normalize": False,
                "min_points": 5,
                "max_points": 15,
                "mode": "exact",
            }

    meas_dict : dict
        Dictionary defining measured-data fields:
            {
                "x": "X",
                "y": "Y",
                "v": "MEDIAN_m",
                "count": "COUNT",
                "slope": "SLOPE",
            }

    model_dict : dict
        Dictionary defining model-data fields:
            {
                "x": "X",
                "y": "Y",
                "v": "mean_head_m",
            }

    Returns
    -------
    result : dict
        Dictionary with:
            - "WL_conditioned"
            - "interp_error"
            - "WLonModel"
            - "SimValOnMeas"
            - "SimValOnModel"
            - "filtered_meas_df"
            - "model_df"
            - "stats"
            - "variograms"
            - "kriging_models"
    """

    # --------------------------------------------------
    # Defaults
    # --------------------------------------------------

    if kriging_params is None:
        kriging_params = {}

    kp = {
        "model": "spherical",
        "n_lags": 15,
        "maxlag_meas": 400000,
        "maxlag_model": "median",
        "maxlag_sim_on_meas": 400000,
        "normalize": False,
        "min_points": 5,
        "max_points": 15,
        "mode": "exact",
    }
    kp.update(kriging_params)

    if meas_dict is None:
        meas_dict = {
            "x": "X",
            "y": "Y",
            "v": "MEDIAN",
            "count": "COUNT",
            "slope": "SLOPE",
        }

    if model_dict is None:
        model_dict = {
            "x": "X",
            "y": "Y",
            "v": "mean_head",
        }

    mx = meas_dict["x"]
    my = meas_dict["y"]
    mv = meas_dict["v"]
    mcount = meas_dict["count"]
    mslope = meas_dict["slope"]

    gx = model_dict["x"]
    gy = model_dict["y"]
    gv = model_dict["v"]

    # --------------------------------------------------
    # Copy and numeric conversion
    # --------------------------------------------------

    meas = meas_df.copy()
    model = model_df.copy()

    for col in [mx, my, mv, mcount, mslope]:
        meas[col] = pd.to_numeric(meas[col], errors="coerce")

    for col in [gx, gy, gv]:
        model[col] = pd.to_numeric(model[col], errors="coerce")

    # --------------------------------------------------
    # Filter measured data
    # --------------------------------------------------

    meas_mask = (
        (meas[mcount] >= n_data) &
        (np.abs(meas[mslope]) < slope_threshold) &
        np.isfinite(meas[mx]) &
        np.isfinite(meas[my]) &
        np.isfinite(meas[mv])
    )

    meas_f = meas.loc[meas_mask].copy()

    if len(meas_f) == 0:
        raise ValueError("No measured points remain after filtering.")

    # --------------------------------------------------
    # Filter model data
    # --------------------------------------------------

    model_mask = (
        np.isfinite(model[gx]) &
        np.isfinite(model[gy]) &
        np.isfinite(model[gv])
    )

    model_f = model.loc[model_mask].copy()

    if len(model_f) == 0:
        raise ValueError("No valid model points remain after filtering.")

    # --------------------------------------------------
    # Extract coordinates and values
    # --------------------------------------------------

    coords_meas = meas_f[[mx, my]].to_numpy(dtype=np.float64)
    vals_meas = meas_f[mv].to_numpy(dtype=np.float64)

    coords_model = model_f[[gx, gy]].to_numpy(dtype=np.float64)
    vals_model = model_f[gv].to_numpy(dtype=np.float64)

    if verbose:
        print("Conditioning model heads to measured heads")
        print(f"Measured points before filtering : {len(meas_df)}")
        print(f"Measured points after filtering  : {len(meas_f)}")
        print(f"Model points used                : {len(model_f)}")
        print(f"n_data                           : {n_data}")
        print(f"slope_threshold                  : {slope_threshold}")
        print(f"Kriging model                    : {kp['model']}")

    # --------------------------------------------------
    # 1. Kriging interpolant for measured heads
    # --------------------------------------------------

    Vmeas = skg.Variogram(
        coords_meas,
        vals_meas,
        maxlag=kp["maxlag_meas"],
        n_lags=kp["n_lags"],
        normalize=kp["normalize"],
        model=kp["model"],
    )

    ok_meas = skg.OrdinaryKriging(
        Vmeas,
        min_points=kp["min_points"],
        max_points=kp["max_points"],
        mode=kp["mode"],
    )

    # --------------------------------------------------
    # 2. Kriging interpolant for model heads
    # --------------------------------------------------

    Vmodel = skg.Variogram(
        coords_model,
        vals_model,
        maxlag=kp["maxlag_model"],
        n_lags=kp["n_lags"],
        normalize=kp["normalize"],
        model=kp["model"],
    )

    ok_model = skg.OrdinaryKriging(
        Vmodel,
        min_points=kp["min_points"],
        max_points=kp["max_points"],
        mode=kp["mode"],
    )

    # --------------------------------------------------
    # 3. Simulated values at measured points
    # --------------------------------------------------

    SimValOnMeas = np.asarray(
        ok_model.transform(
            meas_f[mx].to_numpy(),
            meas_f[my].to_numpy(),
        )
    ).ravel()

    valid_sim_on_meas = np.isfinite(SimValOnMeas)

    if valid_sim_on_meas.sum() == 0:
        raise ValueError("Model kriging returned no valid values at measured points.")

    # --------------------------------------------------
    # 4. Kriging interpolant for simulated values
    #    sampled at measurement locations
    # --------------------------------------------------

    VsimOnMeas = skg.Variogram(
        coords_meas[valid_sim_on_meas],
        SimValOnMeas[valid_sim_on_meas],
        maxlag=kp["maxlag_sim_on_meas"],
        n_lags=kp["n_lags"],
        normalize=kp["normalize"],
        model=kp["model"],
    )

    ok_sim_on_meas = skg.OrdinaryKriging(
        VsimOnMeas,
        min_points=kp["min_points"],
        max_points=kp["max_points"],
        mode=kp["mode"],
    )

    # --------------------------------------------------
    # 5. Reconstruct simulated field from measurement locations
    # --------------------------------------------------

    SimValOnModel = np.asarray(
        ok_sim_on_meas.transform(
            model_f[gx].to_numpy(),
            model_f[gy].to_numpy(),
        )
    ).ravel()

    # --------------------------------------------------
    # 6. Calculate model spatial-detail / interpolation-error term
    # --------------------------------------------------

    interp_error = vals_model - SimValOnModel

    # --------------------------------------------------
    # 7. Interpolate measured water levels to model points
    # --------------------------------------------------

    WLonModel = np.asarray(
        ok_meas.transform(
            model_f[gx].to_numpy(),
            model_f[gy].to_numpy(),
        )
    ).ravel()

    # --------------------------------------------------
    # 8. Final conditioned water level
    # --------------------------------------------------

    WL_conditioned = WLonModel + interp_error

    # --------------------------------------------------
    # Diagnostics
    # --------------------------------------------------

    valid_final = (
        np.isfinite(WL_conditioned) &
        np.isfinite(WLonModel) &
        np.isfinite(interp_error) &
        np.isfinite(vals_model)
    )

    stats = {
        "n_meas_total": len(meas_df),
        "n_meas_used": len(meas_f),
        "n_model_total": len(model_df),
        "n_model_used": len(model_f),
        "n_valid_conditioned": int(valid_final.sum()),
        "meas_mean": float(np.nanmean(vals_meas)),
        "meas_std": float(np.nanstd(vals_meas, ddof=1)),
        "model_mean": float(np.nanmean(vals_model)),
        "model_std": float(np.nanstd(vals_model, ddof=1)),
        "conditioned_mean": float(np.nanmean(WL_conditioned)),
        "conditioned_std": float(np.nanstd(WL_conditioned, ddof=1)),
        "interp_error_mean": float(np.nanmean(interp_error)),
        "interp_error_std": float(np.nanstd(interp_error, ddof=1)),
        "interp_error_min": float(np.nanmin(interp_error)),
        "interp_error_max": float(np.nanmax(interp_error)),
        "mean_shift_model_to_conditioned": float(
            np.nanmean(WL_conditioned - vals_model)
        ),
        "std_shift_model_to_conditioned": float(
            np.nanstd(WL_conditioned - vals_model, ddof=1)
        ),
    }

    if verbose:
        print("Conditioning completed")
        print(f"Valid conditioned values         : {stats['n_valid_conditioned']}")
        print(f"Measured mean                    : {stats['meas_mean']:.3f}")
        print(f"Model mean                       : {stats['model_mean']:.3f}")
        print(f"Conditioned mean                 : {stats['conditioned_mean']:.3f}")
        print(f"Mean shift conditioned - model   : {stats['mean_shift_model_to_conditioned']:.3f}")
        print(f"Interpolation-error std          : {stats['interp_error_std']:.3f}")

    # --------------------------------------------------
    # Return
    # --------------------------------------------------

    result = {
        "WL_conditioned": WL_conditioned,
        "interp_error": interp_error,
        "WLonModel": WLonModel,
        "SimValOnMeas": SimValOnMeas,
        "SimValOnModel": SimValOnModel,
        "filtered_meas_df": meas_f,
        "model_df": model_f,
        "stats": stats,
        "variograms": {
            "measured": Vmeas,
            "model": Vmodel,
            "sim_on_meas": VsimOnMeas,
        },
        "kriging_models": {
            "measured": ok_meas,
            "model": ok_model,
            "sim_on_meas": ok_sim_on_meas,
        },
    }

    return result