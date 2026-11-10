"""Check the stream / bucket / lake temperature output variables.

With water_temperature = 2 the model writes
  CHRTOUT  T_eq, T_airSmooth, T_upstream, T_lateral, T_gwInflow   [C]
  GWOUT    T_bucket                                               [C]
  LAKEOUT  T_inflow, T_outflow                                    [C]
(-9999 = not applicable).  This script

1) (optional) compares every variable that RUN_OLD and RUN_NEW both wrote in
   CHRTOUT, GWOUT, LAKEOUT and LDASOUT: they must be bit-identical (the new
   outputs must not change the model);
2) checks that the new variables are present in RUN_NEW and consistent:
   - at each lake pour-point reach, CHRTOUT T_upstream = LAKEOUT T_inflow
     (same flow-weighted upstream temperature) and riverTemp = T_outflow;
   - CHRTOUT T_gwInflow = GWOUT T_bucket for reaches whose bucket has a
     valid temperature (GWOUT feature_id = reach ComID in UDMP);
   - riverTemp lies between the reach inflow temperature and T_eq where the
     reach has only upstream inflow and no lateral inflow (relaxation toward
     T_eq), within a small tolerance.

usage: python check_temp_outputs.py RUN_NEW [RUN_OLD]
"""
import sys, os, glob
import numpy as np
import netCDF4 as nc

NEW_CH = ['T_eq', 'T_airSmooth', 'T_upstream', 'T_lateral', 'T_gwInflow']


def arr(v):
    return np.ma.filled(np.asarray(v[:], dtype=float), np.nan)


def compare(old, new):
    ok = True
    print(f'== common variables, {new} vs {old}')
    for pat in ('*.CHRTOUT_DOMAIN1', '*.GWOUT_DOMAIN1', '*.LAKEOUT_DOMAIN1', '*.LDASOUT_DOMAIN1'):
        nfile = ndiff = 0; diffvars = set()
        for g in sorted(glob.glob(os.path.join(new, 'OUTDIR', pat))):
            f = os.path.join(old, 'OUTDIR', os.path.basename(g))
            if not os.path.exists(f):
                continue
            a, b = nc.Dataset(f), nc.Dataset(g)
            nfile += 1
            for v in a.variables:
                if v in b.variables and a[v].dtype.kind in 'fi':
                    if not np.array_equal(arr(a[v]), arr(b[v]), equal_nan=True):
                        diffvars.add(v)
        print(f'  {pat:18s} {nfile:4d} files; differing variables: {sorted(diffvars) or "none"}')
        ok &= not diffvars
    return ok


def consistency(run):
    ok = True
    out = os.path.join(run, 'OUTDIR')
    ch = sorted(glob.glob(os.path.join(out, '*.CHRTOUT_DOMAIN1')))
    lk = {os.path.basename(f)[:12]: f for f in glob.glob(os.path.join(out, '*.LAKEOUT_DOMAIN1'))}
    gw = {os.path.basename(f)[:12]: f for f in glob.glob(os.path.join(out, '*.GWOUT_DOMAIN1'))}
    d = nc.Dataset(ch[-1])
    missing = [v for v in NEW_CH if v not in d.variables]
    print(f'== {run}: CHRTOUT new variables missing: {missing or "none"}')
    ok &= not missing
    if lk:
        missing = [v for v in ('T_inflow', 'T_outflow') if v not in nc.Dataset(next(iter(lk.values()))).variables]
        print(f'   LAKEOUT new variables missing: {missing or "none"}'); ok &= not missing
    if gw:
        missing = [v for v in ('T_bucket',) if v not in nc.Dataset(next(iter(gw.values()))).variables]
        print(f'   GWOUT new variables missing: {missing or "none"}'); ok &= not missing
    if not ok:
        return ok

    # Route_Link: lake pour-point reaches (reaches whose 'to' leaves a lake)
    rl = nc.Dataset(os.path.join(run, 'DOMAIN', 'Route_Link.nc'))
    link = rl['link'][:]; lakeid = rl['NHDWaterbodyComID'][:]
    lp = None
    worst = {'lake_in': 0.0, 'lake_out': 0.0, 'gw': 0.0}
    nlake = ngw = 0; nbound = nviol = 0
    for f in ch[1:]:
        key = os.path.basename(f)[:12]
        c = nc.Dataset(f)
        fid = c['feature_id'][:]
        T = {v: arr(c[v]) for v in NEW_CH + ['riverTemp', 'q_lateral', 'streamflow']}
        if key in lk:
            l = nc.Dataset(lk[key])
            lfid = l['feature_id'][:]
            tin, tout = arr(l['T_inflow']), arr(l['T_outflow'])
            # pour-point reach of each lake: the reach carrying the lake's ComID
            # with the 'to' outside the lake; in LAKEOUT order use the LAKEPARM
            # outlet via Route_Link: the reach with NHDWaterbodyComID == lake
            # and largest streamflow is not robust, so match on T_outflow ==
            # riverTemp candidates.
            for k, L in enumerate(lfid):
                cand = np.where(lakeid == L)[0]
                if len(cand) == 0 or tout[k] < -9000:
                    continue
                cidx = np.searchsorted(np.sort(fid), link[cand])
                order = np.argsort(fid)
                idx = order[np.clip(cidx, 0, len(fid) - 1)]
                rt = T['riverTemp'][idx]
                j = idx[np.argmin(np.abs(rt - tout[k]))]
                worst['lake_out'] = max(worst['lake_out'], abs(T['riverTemp'][j] - tout[k]))
                if tin[k] > -9000 and T['T_upstream'][j] > -9000:
                    worst['lake_in'] = max(worst['lake_in'], abs(T['T_upstream'][j] - tin[k]))
                nlake += 1
        if key in gw:
            g = nc.Dataset(gw[key])
            gfid = g['feature_id'][:]; tb = arr(g['T_bucket'])
            m = dict(zip(gfid, tb))
            tg = np.array([m.get(x, -9999.0) for x in fid])
            ok_b = (tg > -9000) & (T['T_gwInflow'] > -9000)
            if ok_b.any():
                # T_gwInflow is the bucket temperature clamped to [0, 40] C
                worst['gw'] = max(worst['gw'], float(np.max(np.abs(np.clip(tg[ok_b], 0, 40) - T['T_gwInflow'][ok_b]))))
                ngw += int(ok_b.sum())
        # relaxation bound for reaches with no lateral inflow
        sel = (T['q_lateral'] <= 0) & (T['T_upstream'] > -9000) & (T['T_eq'] > -9000)
        lo = np.minimum(T['T_upstream'][sel], T['T_eq'][sel]) - 0.05
        hi = np.maximum(T['T_upstream'][sel], T['T_eq'][sel]) + 0.05
        nbound += int(sel.sum()); nviol += int(((T['riverTemp'][sel] < lo) | (T['riverTemp'][sel] > hi)).sum())
    print(f'   lake pour points ({nlake} lake-times): max |riverTemp - T_outflow| = {worst["lake_out"]:.3g} C, '
          f'max |T_upstream - T_inflow| = {worst["lake_in"]:.3g} C')
    print(f'   buckets ({ngw} reach-times): max |T_gwInflow - clamp(T_bucket)| = {worst["gw"]:.3g} C')
    print(f'   reaches with no lateral inflow ({nbound} reach-times): riverTemp outside [T_upstream, T_eq] '
          f'by more than 0.05 C: {nviol}')
    ok &= worst['lake_out'] < 1e-3 and worst['lake_in'] < 1e-3 and worst['gw'] < 1e-3
    return ok


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ok = True
    if len(sys.argv) > 2:
        ok &= compare(sys.argv[2], sys.argv[1])
    ok &= consistency(sys.argv[1])
    print('RESULT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
