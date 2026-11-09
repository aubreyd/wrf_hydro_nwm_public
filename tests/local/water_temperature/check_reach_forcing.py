"""Check that reach air temperature in the stream model equals the catchment mean.

For headwater reaches (no upstream reach, not a lake outlet) whose lateral inflow
is below min_flow (1e-4 m3/s), run_stream_temp_eq sets T_in = T_e and therefore
T_out = T_e = T_air_smooth + SW_net/K_sw.  At night (catchment SW == 0) this is
exactly T_air_smooth, the 24-h EMA of the reach's air temperature.

This script rebuilds T_air_smooth offline from the LDASIN forcing and the
spatialweights mapping (area-weighted catchment mean of T2D), and compares it
with the model's riverTemp at those reaches/hours.

usage: python check_reach_forcing.py RUNDIR NHOURS [START]
  RUNDIR  run directory containing DOMAIN/, FORCING/ (LDASIN) and OUTDIR/ (CHRTOUT)
  NHOURS  number of hourly outputs to check after START
  START   cold-start time, default 2019-07-01T00 (run must be a cold start,
          since the EMA is seeded from forcing at the first step)

Expected result for correct code: the flip=False lag=+0 line shows max|diff| ~ 0.
"""
import sys, glob, os
import numpy as np
import xarray as xr

run, nh = sys.argv[1], int(sys.argv[2])
start = np.datetime64(sys.argv[3] if len(sys.argv) > 3 else "2019-07-01T00", "h")
dom = os.path.join(run, "DOMAIN")
forc = os.path.join(run, "FORCING")
out = os.path.join(run, "OUTDIR")
AGG = 4

rl = xr.open_dataset(f"{dom}/Route_Link.nc")
links = rl.link.values
to = rl["to"].values
headwater = ~np.isin(links, to)

sw = xr.open_dataset(f"{dom}/spatialweights.nc")
ids = sw.IDmask.values
ii = sw.i_index.values - 1
jj = sw.j_index.values - 1
w = sw.regridweight.values

times = [start + np.timedelta64(h, "h") for h in range(nh + 1)]

def fname(t, kind):
    s = str(t.astype("datetime64[m]")).replace("-", "").replace("T", "").replace(":", "")[:12]
    return f"{s}.{kind}_DOMAIN1"

def catchment_mean(field_lsm, flip):
    """field_lsm: (y, x) as stored in the file (row 0 = first row of file)."""
    f = field_lsm[::-1, :] if flip else field_lsm
    vals = f[jj // AGG, ii // AGG]
    num = np.bincount(np.searchsorted(uniq, ids), weights=w * vals, minlength=uniq.size)
    den = np.bincount(np.searchsorted(uniq, ids), weights=w, minlength=uniq.size)
    return num / den

uniq = np.unique(ids)

def run_check(flip, lag):
    """lag: forcing file time offset (h) relative to the step's output time."""
    ema = None
    alpha = 1.0 / 24.0
    diffs = []
    for k, t in enumerate(times[1:], start=1):
        tf = t + np.timedelta64(lag, "h")
        fd = xr.open_dataset(os.path.join(forc, fname(tf, "LDASIN")))
        T = catchment_mean(fd.T2D.values[0], flip) - 273.15
        S = catchment_mean(fd.SWDOWN.values[0], flip)
        ema = T.copy() if ema is None else alpha * T + (1 - alpha) * ema
        c = xr.open_dataset(os.path.join(out, fname(t, "CHRTOUT")))
        fid = c.feature_id.values
        pos = {v: n for n, v in enumerate(uniq)}
        rt = c.riverTemp.values
        latq = c.qSfcLatRunoff.values + c.qBucket.values
        hw = dict(zip(links, headwater))
        for n, f in enumerate(fid):
            if f not in pos or not hw[f] or latq[n] > 1e-4:
                continue
            p = pos[f]
            if S[p] > 0:          # night only
                continue
            diffs.append(rt[n] - ema[p])
    d = np.array(diffs)
    return d

for flip in (False, True):
    for lag in (0, -1):
        d = run_check(flip, lag)
        if d.size:
            print(f"flip={flip!s:5} lag={lag:+d}  n={d.size:5d}  mean|diff|={np.abs(d).mean():.4f}  "
                  f"max|diff|={np.abs(d).max():.4f}  frac<0.01C={np.mean(np.abs(d)<0.01):.2f}")
        else:
            print(f"flip={flip} lag={lag}: no qualifying reach-hours")
