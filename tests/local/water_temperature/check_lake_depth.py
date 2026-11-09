"""Check that the 1-D lake model depth follows the level-pool reservoir elevation.

With sf_lake_physics=1, the Noah-MP driver shifts the lake model depth
(LAKE_DEPTH2D) on each lake's LAKE_ID cells by the level-pool water-elevation
change (RESHT) since the previous LSM step.  The LSM at time t sees RESHT from
the previous routing step, i.e. the LAKEOUT water_sfc_elev written just before t.

LAKE_DEPTH2D is fill at t0 (it is set after the first LSM step), so the first
valid LDASOUT time t1 is the reference.  With H_seen(t) = water_sfc_elev of the
LAKEOUT written just before t, for every lake cell and LDASOUT time t >= t1:
    LAKE_DEPTH2D(t) - LAKE_DEPTH2D(t1) == H_seen(t) - H_seen(t1)
(this also holds when the LSM step is shorter than the LDASOUT interval).
unless the 1 m depth floor is active.  For correct code the max difference is ~0
at any LSM timestep and MPI rank count.

usage: python check_lake_depth.py RUNDIR WRFINPUT
  RUNDIR    cold-start run directory with OUTDIR/ containing LDASOUT and LAKEOUT
  WRFINPUT  the HRLDAS setup file holding LAKE_ID (as in namelist.hrldas)
"""
import sys, glob, os
import numpy as np
import xarray as xr

run, wrfinput = sys.argv[1], sys.argv[2]
lake_id = xr.open_dataset(wrfinput)["LAKE_ID"].values.squeeze()
out = os.path.join(run, "OUTDIR")

ldas = sorted(glob.glob(os.path.join(out, "*.LDASOUT_DOMAIN1")))
lake = sorted(glob.glob(os.path.join(out, "*.LAKEOUT_DOMAIN1")))
lake_stamps = [os.path.basename(f)[:12] for f in lake]

L0 = xr.open_dataset(lake[0])
t1 = os.path.basename(ldas[1])[:12]
D0 = xr.open_dataset(ldas[1])["LAKE_DEPTH2D"].values.squeeze()        # reference depth at t1
H0 = xr.open_dataset(lake[lake_stamps.index(t1) - 1]).water_sfc_elev.values  # elevation seen at t1

worst = 0.0
n = 0
nbad = 0
for f in ldas[1:]:
    t = os.path.basename(f)[:12]
    if t not in lake_stamps:
        continue
    Lprev = xr.open_dataset(lake[lake_stamps.index(t) - 1])
    D = xr.open_dataset(f)["LAKE_DEPTH2D"].values.squeeze()
    for k, fid in enumerate(L0.feature_id.values):
        m = lake_id == fid
        if not m.any():
            continue
        dH = float(Lprev.water_sfc_elev.values[k] - H0[k])
        diff = np.abs((D[m] - D0[m]) - dH)
        if not np.all(np.isfinite(diff)):
            nbad += 1
            continue
        worst = max(worst, float(diff.max()))
        n += 1

Lend = xr.open_dataset(lake[-1])
print(f"{n} lake-times checked, {nbad} with non-finite values")
if nbad or n == 0:
    print("FAIL: missing or non-finite lake depths/elevations")
print("final elevation change by lake (m):",
      dict(zip(L0.feature_id.values.tolist(),
               np.round(Lend.water_sfc_elev.values - L0.water_sfc_elev.values, 3).tolist())))
print(f"max |lake depth change - level-pool elevation change| = {worst:.5f} m")
