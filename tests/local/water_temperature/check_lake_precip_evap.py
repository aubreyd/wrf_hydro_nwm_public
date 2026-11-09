"""Check that level-pool reservoirs receive the lake-mean precipitation and evaporation.

The Noah-MP driver sums RAINBL and QFX over each lake's LAKE_ID cells, reduces the
sums and cell counts across MPI ranks, and the coupler divides to get lake means
that are passed to the reservoirs (and written to LAKEOUT as lake_precip/lake_evap).

This script recomputes those lake means from the LDASOUT 2-D lake fields
(LAKE_PRECIP2D = RAINBL*3600/dt and LAKE_EVAP2D = QFX*3600, both mm/h) and the
wrfinput LAKE_ID, and compares with LAKEOUT lake_precip/lake_evap (mm/h).  For
correct code the two agree to float precision at every output time, independent
of MPI rank count and LSM timestep.  LAKEOUT times without a matching LDASOUT
(e.g. sub-hourly hydro output) are skipped.

usage: python check_lake_precip_evap.py RUNDIR WRFINPUT
  RUNDIR    run directory with OUTDIR/ containing LDASOUT and LAKEOUT files
  WRFINPUT  the HRLDAS setup file holding LAKE_ID (as in namelist.hrldas)
"""
import sys, glob, os
import numpy as np
import xarray as xr

run, wrfinput = sys.argv[1], sys.argv[2]
lake_id = xr.open_dataset(wrfinput)["LAKE_ID"].values.squeeze()

lakeouts = sorted(glob.glob(os.path.join(run, "OUTDIR", "*.LAKEOUT_DOMAIN1")))
worst_p = worst_e = 0.0
nbad = 0
rows = []
for lf in lakeouts[1:]:            # skip t0 (no step taken yet)
    stamp = os.path.basename(lf)[:12]
    ldas = os.path.join(run, "OUTDIR", f"{stamp}.LDASOUT_DOMAIN1")
    if not os.path.exists(ldas):
        continue
    L = xr.open_dataset(lf)
    D = xr.open_dataset(ldas)
    P2 = D["LAKE_PRECIP2D"].values.squeeze()
    E2 = D["LAKE_EVAP2D"].values.squeeze()
    for k, fid in enumerate(L.feature_id.values):
        m = lake_id == fid
        if not m.any():
            continue
        p_off = P2[m].mean()
        e_off = E2[m].mean()
        p_mod = float(L.lake_precip.values[k])
        e_mod = float(L.lake_evap.values[k])
        if not np.all(np.isfinite([p_off, e_off, p_mod, e_mod])):
            nbad += 1
            continue
        worst_p = max(worst_p, abs(p_mod - p_off))
        worst_e = max(worst_e, abs(e_mod - e_off))
        rows.append((stamp, int(fid), int(m.sum()), p_off, p_mod, e_off, e_mod))

print(f"{'time':12} {'lake':>8} {'ncell':>5} {'P_ldas':>9} {'P_lakeout':>10} {'E_ldas':>9} {'E_lakeout':>10}")
for r in rows[: 4 * 4]:
    print(f"{r[0]:12} {r[1]:8d} {r[2]:5d} {r[3]:9.4f} {r[4]:10.4f} {r[5]:9.4f} {r[6]:10.4f}")
print(f"... {len(rows)} lake-times checked, {nbad} with non-finite values")
if nbad or not rows:
    print("FAIL: missing or non-finite lake precip/evap")
print(f"max |precip diff| = {worst_p:.5f} mm   max |evap diff| = {worst_e:.5f} mm/h")
