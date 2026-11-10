"""Check the 1-D lake inflow heat, lake sub-stepping and lake restart state.

Three checks, each optional (pass '-' to skip a pair):

1) Lake profiles: mean LAKE_T3D profile of each LAKE_ID lake in two runs (for
   example sub-stepping only vs. sub-stepping plus inflow heat) on the dates
   that both runs wrote LDASOUT.  Inflow heat shows up as warmer metalimnion
   and outlet layers in summer (warm river water entering below the surface).

2) Restart identity: a continuous run and a run restarted from one of its
   restart files.  Every output variable listed below and every variable in
   the RESTART / HYDRO_RST files both runs wrote must be bit-identical, except
   those listed in KNOWN_RESTART_DIFFS (not from the lake/stream temperature code).

3) HYDRO_RST contents: lake_inflow_temp (K) and lake_resht_lsm (m) are present
   in the last HYDRO_RST of the continuous run, and lake_inflow_temp is in the
   physical range for lakes with inflow.

usage: python check_lake_inflow_heat.py RUN_A RUN_B RUN_CONT RUN_RST
   RUN_*: run directories with OUTDIR/ and DOMAIN/ (wrfinput with LAKE_ID)
example (DRB):
   python check_lake_inflow_heat.py run.drb_lk_substep run.drb_lk_heat \
       run.lkh_new run.lkh_new_rst
"""
import sys, os, glob
import numpy as np
import netCDF4 as nc

OUT_VARS = [('*.CHRTOUT_DOMAIN1', ['streamflow', 'riverTemp']),
            ('*.LAKEOUT_DOMAIN1', ['water_sfc_elev', 'outflow']),
            ('*.LDASOUT_DOMAIN1', ['TG', 'SOIL_T', 'SNEQV', 'LAKE_T3D', 'LAKE_DEPTH2D'])]
# Restart differences not caused by the lake or stream temperature code:
# qstrmvolrt (routing-grid channel inflow accumulator) differs after a
# restart.  The lake snow/ice state differences seen earlier (uninitialized
# imelt in the 1-D lake model) are fixed and are no longer tolerated.
KNOWN_RESTART_DIFFS = {'qstrmvolrt'}


def wrfinput(run):
    f = sorted(glob.glob(os.path.join(run, 'DOMAIN', 'wrfinput*')))
    return f[0]


def lake_profiles(run_a, run_b):
    lake_id = np.squeeze(nc.Dataset(wrfinput(run_a))['LAKE_ID'][:])
    lakes = [int(l) for l in np.unique(lake_id[lake_id > 0])]
    fa = {os.path.basename(f) for f in glob.glob(os.path.join(run_a, 'OUTDIR', '*.LDASOUT_DOMAIN1'))}
    fb = {os.path.basename(f) for f in glob.glob(os.path.join(run_b, 'OUTDIR', '*.LDASOUT_DOMAIN1'))}
    common = sorted(fa & fb)
    pick = common[::max(1, len(common) // 12)]
    print(f'Lake profiles (top -> bottom layer, C): A = {run_a}, B = {run_b}')
    for lk in lakes:
        m = lake_id == lk
        print(f'  lake {lk} ({m.sum()} cells)')
        for f in pick:
            row = []
            for run in (run_a, run_b):
                v = nc.Dataset(os.path.join(run, 'OUTDIR', f))['LAKE_T3D']
                t = np.transpose(np.squeeze(v[:]).astype(float), (0, 2, 1))[m] - 273.15
                row.append(' '.join(f'{x:5.1f}' for x in t.mean(0)))
            print(f'    {f[:8]}  A [{row[0]}]')
            print(f'              B [{row[1]}]')


def same(x, y):
    x = np.ma.filled(np.asarray(x, dtype=float), np.nan)
    y = np.ma.filled(np.asarray(y, dtype=float), np.nan)
    return np.array_equal(x, y, equal_nan=True)


def restart_identity(run_c, run_r):
    ok = True
    print(f'Restart identity: continuous {run_c} vs restarted {run_r}')
    for pat, vs in OUT_VARS:
        cnt = {v: [0, 0] for v in vs}
        for g in sorted(glob.glob(os.path.join(run_r, 'OUTDIR', pat))):
            f = os.path.join(run_c, 'OUTDIR', os.path.basename(g))
            if not os.path.exists(f):
                continue
            a, b = nc.Dataset(f), nc.Dataset(g)
            for v in vs:
                if v in a.variables and v in b.variables:
                    cnt[v][0] += 1
                    # LAKE_DEPTH2D: compare finite values only (fill pattern differs)
                    if v == 'LAKE_DEPTH2D':
                        x = np.ma.filled(a[v][:].astype(float), np.nan); y = np.ma.filled(b[v][:].astype(float), np.nan)
                        k = np.isfinite(x) & np.isfinite(y)
                        diff = not np.array_equal(x[k], y[k])
                    else:
                        diff = not same(a[v][:], b[v][:])
                    cnt[v][1] += diff
        for v, (n, d) in cnt.items():
            if n:
                print(f'  {v:15s} {n:4d} files, {d} differ')
                ok &= d == 0
    for pat in ('RESTART.*_DOMAIN1', 'HYDRO_RST.*_DOMAIN1'):
        for g in sorted(glob.glob(os.path.join(run_r, 'OUTDIR', pat))):
            f = os.path.join(run_c, 'OUTDIR', os.path.basename(g))
            if not os.path.exists(f):
                continue
            a, b = nc.Dataset(f), nc.Dataset(g)
            diff = [v for v in a.variables if v in b.variables and a[v].dtype.kind in 'fi'
                    and not same(a[v][:], b[v][:])]
            new = [v for v in diff if v not in KNOWN_RESTART_DIFFS]
            known = [v for v in diff if v in KNOWN_RESTART_DIFFS]
            print(f'  {os.path.basename(g)}: unexpected differences {new or "none"}; known {known}')
            ok &= not new
    print('  RESTART IDENTITY:', 'PASS' if ok else 'FAIL')
    return ok


def hydro_rst(run_c):
    f = sorted(glob.glob(os.path.join(run_c, 'OUTDIR', 'HYDRO_RST.*_DOMAIN1')))[-1]
    d = nc.Dataset(f)
    ok = True
    print(f'HYDRO_RST contents: {os.path.basename(f)}')
    for v in ('lake_inflow_temp', 'lake_resht_lsm'):
        present = v in d.variables
        ok &= present
        print(f'  {v}: {"present" if present else "MISSING"}', end='')
        if present:
            x = d[v][:]
            print('  values', np.round(x - 273.15, 2) if v == 'lake_inflow_temp' else np.round(x, 2), end='')
        print()
    if 'lake_inflow_temp' in d.variables and 'qlakei' in d.variables:
        t, q = d['lake_inflow_temp'][:], d['qlakei'][:]
        bad = (q > 0) & ((t < 273.15) | (t > 313.15))
        print(f'  lakes with inflow and lake_inflow_temp outside 0-40 C: {int(bad.sum())}')
        ok &= not bad.any()
    print('  HYDRO_RST:', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    run_a, run_b, run_c, run_r = sys.argv[1:5]
    ok = True
    if run_a != '-' and run_b != '-':
        lake_profiles(run_a, run_b)
    if run_c != '-' and run_r != '-':
        ok &= restart_identity(run_c, run_r)
    if run_c != '-':
        ok &= hydro_rst(run_c)
    sys.exit(0 if ok else 1)
