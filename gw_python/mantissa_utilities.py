import pandas as pd
import numpy as np

def read_rs_data_merge(fname):
    tmp_c = pd.read_csv(fname + "RSC.dat.gz",compression='gzip',sep=",")
    tmp_c.columns = tmp_c.columns.str.strip()
    tmp_l = pd.read_csv(fname + "RSL.dat.gz",compression='gzip',sep=",")
    tmp_l.columns = tmp_l.columns.str.strip()
    tmp_f = pd.read_csv(fname + "RSF.dat.gz",compression='gzip',sep=",")
    tmp_f.columns = tmp_f.columns.str.strip()
    merged = (tmp_c[["Eid", "Sid", "W"]]
    .merge(tmp_l[["Eid", "Sid", "UrfI", "UrfJ", "hru", "hru_idx", "UnsatD", "UnsatR", "SrfPrc", "RivInfl"]], on=["Eid", "Sid"], how="inner")
    .merge(tmp_f[["Eid", "Sid", "Aid"]], on=["Eid", "Sid"], how="inner")
    )
    c = tmp_c.drop(columns=["W"]).to_numpy()
    l = tmp_l.drop(columns=["UrfI", "UrfJ", "hru", "hru_idx", "UnsatD", "UnsatR", "SrfPrc", "RivInfl"]).to_numpy()
    f = tmp_c.drop(columns=["W"]).to_numpy()
    return merged, c, l, f


def read_swat_output_file(filename, watershed_order = None, field_names = None):
    if watershed_order is None:
        watershed_order = ['SACV', 'SJV', 'TLB', 'DM_SUB', 'Kings']
    if field_names is None:
        field_names = ['hru_num','irrSW_mm','irrGW_mm','irrsaltSW_kgha','irrsaltGW_kgha',
                       'fertsalt_kgha','dssl_kgha','Qsalt_kgha','uptk_kgha','pGW',
                       'dSoilSalt_kgha','perc_mm','Salt_perc_ppm']

    swatTAB = pd.read_csv(filename)

    # Check for missing fields and add them
    if 'yr' not in swatTAB:
        swatTAB['yr'] = 1990

    if 'pGW' not in swatTAB:
        swatTAB['pGW'] = 1.0

    # Append watershed numeric code
    swatTAB['WatershedNum'] = np.nan

    # Assign watershed numeric values
    for i, ws in enumerate(watershed_order, start=1):
        swatTAB.loc[swatTAB['Watershed'] == ws, 'WatershedNum'] = i * 1_000_000

    # Compute hru_num
    swatTAB['hru_num'] = swatTAB['WatershedNum'] + swatTAB['gis_id']

    # Get unique years and HRUs
    swat_years = swatTAB['yr'].unique()
    swat_hrus_uniq = swatTAB['hru_num'].unique()

    num_hrus = swat_hrus_uniq.shape[0]
    num_years = swat_years.shape[0]

    swat = {name: np.full((num_hrus, num_years), np.nan) for name in field_names}


    # Map HRU rows per year to consistent row order (swat_hrus_uniq)
    hru_indexer = {h: i for i, h in enumerate(swat_hrus_uniq)}

    for j, yr in zip(range(num_years), swat_years):
        idx = swatTAB.index[swatTAB['yr'] == yr]
        tmp = swatTAB.loc[idx]

        for fname in field_names:
            if fname in tmp.columns:
                swat[fname][:, j] = tmp[fname].values

    for fname in field_names:
        swat[fname] = np.nan_to_num(swat[fname], nan=0)

    swat["hrus"] = swat.pop("hru_num").astype(int)

    return swatTAB, swat