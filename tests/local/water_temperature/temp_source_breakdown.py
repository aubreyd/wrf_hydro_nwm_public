"""Basin-scale breakdown of simulated stream temperature at USGS gauges.

For each gauge reach and each hourly output in a period, the reaches upstream
of (and including) the gauge reach are found from Route_Link ('link' -> 'to').
Over those reaches the script computes
  Q_gw, T_gw     total groundwater (bucket) inflow [m3/s] and its flow-weighted
                 temperature (GWOUT outflow x CHRTOUT T_gwInflow)
  Q_sfc, T_sfc   total surface (non-bucket) lateral inflow, q_lateral minus the
                 bucket outflow, and its flow-weighted temperature (T_lateral)
  T_src          flow-weighted temperature of all water entering the network
                 upstream of the gauge (groundwater + surface)
  T_eq, T_air    length-weighted mean equilibrium temperature and smoothed air
                 temperature of the upstream reaches
and, at the gauge reach, the simulated riverTemp and streamflow.  Results are
averaged to days (CHRTOUT times are end-of-hour and shifted back 1 h, as in
eval_usgs_seasonal.py), merged with USGS daily mean temperature, and
summarised by gauge and month.

Reading the breakdown: in this model a reach temperature is its inflow
temperature relaxed toward T_eq, so the gauge temperature lies between the
source temperature (T_src, upstream mixing) and the basin T_eq.  If T_src and
T_eq are both above the observation the bias is not from the inputs; if T_src
is far below and the gauge sits near T_src, the sources are too cold or the
exchange with the atmosphere is too weak (Kexch); if T_eq itself is below the
observation, the equilibrium (air temperature, shortwave/Ksw) is too cold.

Needs a run with water_temperature = 2 and output_gw = 1 (CHRTOUT T_* and
GWOUT outflow / T_bucket).

usage: python temp_source_breakdown.py OUTDIR OBSDIR START END OUT.csv [GAUGES]
  GAUGES  comma-separated USGS site numbers (default: the unregulated gauges in
          HUC8 02040101/02040102 that are not downstream of Cannonsville or
          Pepacton)
  OUT.csv daily table (one row per gauge-day); the monthly summary is printed
"""
import sys, os, glob
from collections import defaultdict
import numpy as np, pandas as pd
import netCDF4 as nc

outdir, obsdir, start, end, outcsv = sys.argv[1:6]
gauges_arg = sys.argv[6] if len(sys.argv) > 6 else None
route_link = os.path.join(os.path.dirname(os.path.abspath(outdir)), 'DOMAIN', 'Route_Link.nc')
LAKES = {2613174: 'Cannonsville', 1748473: 'Pepacton'}
EVAL_HUC8 = {'02040101', '02040102'}

# --- network ---
rl = nc.Dataset(route_link)
link = rl['link'][:].astype(np.int64); to = rl['to'][:].astype(np.int64)
length = rl['Length'][:].astype(float)
wbody = rl['NHDWaterbodyComID'][:].astype(np.int64) if 'NHDWaterbodyComID' in rl.variables else np.zeros_like(link)
pos = {c: i for i, c in enumerate(link)}
ups = defaultdict(list)
for i, t in enumerate(to):
    if t in pos:
        ups[pos[t]].append(i)

def upstream(i):
    seen, stack = {i}, [i]
    while stack:
        k = stack.pop()
        for u in ups[k]:
            if u not in seen:
                seen.add(u); stack.append(u)
    return np.array(sorted(seen))

def downstream_of_lake(i):
    """True if reach i is downstream of one of the evaluated lakes."""
    up = upstream(i)
    return any(wbody[u] in LAKES for u in up)

sites = pd.read_csv(f'{obsdir}/usgs_sites.csv', dtype={'site_no': str, 'huc8': str})
sites = sites[sites.huc8.isin(EVAL_HUC8) & sites.comid.isin(pos)]
if gauges_arg:
    sites = sites[sites.site_no.isin(gauges_arg.split(','))]
else:
    sites = sites[[not downstream_of_lake(pos[c]) for c in sites.comid]]
obs = pd.read_csv(f'{obsdir}/usgs_daily_temperature.csv', dtype={'site_no': str}, parse_dates=['date'])
obs = obs[(obs.date >= start) & (obs.date <= end) & obs.site_no.isin(sites.site_no)]
sites = sites[sites.site_no.isin(obs.site_no.unique())]
print(f'{len(sites)} gauges:', ', '.join(f'{s} ({n[:30]})' for s, n in zip(sites.site_no, sites.station_nm)))
upidx = {s: upstream(pos[c]) for s, c in zip(sites.site_no, sites.comid)}
gidx = {s: pos[c] for s, c in zip(sites.site_no, sites.comid)}

# --- hourly loop ---
t0 = pd.Timestamp(start); t1 = pd.Timestamp(end) + pd.Timedelta(days=1)
files = []
for f in sorted(glob.glob(os.path.join(outdir, '*.CHRTOUT_DOMAIN1'))):
    t = pd.to_datetime(os.path.basename(f)[:12], format='%Y%m%d%H%M') - pd.Timedelta(hours=1)
    if t0 <= t < t1:
        files.append((t, f))
print(f'{len(files)} hourly files')
rows = []
order = None
for n, (t, f) in enumerate(files):
    c = nc.Dataset(f)
    fid = c['feature_id'][:].astype(np.int64)
    if order is None:
        o = np.array([np.where(fid == L)[0][0] for L in link])   # Route_Link order
        order = o
    v = {k: np.ma.filled(c[k][:].astype(float), np.nan)[order]
         for k in ('streamflow', 'q_lateral', 'riverTemp', 'T_eq', 'T_airSmooth', 'T_lateral', 'T_gwInflow')}
    g = nc.Dataset(f.replace('CHRTOUT', 'GWOUT'))
    gq = dict(zip(g['feature_id'][:].astype(np.int64), np.ma.filled(g['outflow'][:].astype(float), 0.0)))
    qgw = np.array([gq.get(L, 0.0) for L in link])
    qsfc = np.maximum(0.0, v['q_lateral'] - qgw)
    for k in ('T_eq', 'T_airSmooth', 'T_lateral', 'T_gwInflow'):
        v[k][v[k] < -9000] = np.nan
    for s in sites.site_no:
        u = upidx[s]
        def wmean(T, w):
            m = np.isfinite(T[u]) & (w[u] > 0)
            return np.sum(T[u][m] * w[u][m]) / np.sum(w[u][m]) if m.any() else np.nan
        Qgw, Qs = qgw[u].sum(), qsfc[u].sum()
        Tgw, Ts = wmean(v['T_gwInflow'], qgw), wmean(v['T_lateral'], qsfc)
        Tsrc = (np.nansum([Qgw * Tgw if Qgw > 0 else 0, Qs * Ts if Qs > 0 else 0]) / (Qgw + Qs)) if Qgw + Qs > 0 else np.nan
        rows.append(dict(site_no=s, time=t, Q_gauge=v['streamflow'][gidx[s]], T_model=v['riverTemp'][gidx[s]],
                         Q_gw=Qgw, T_gw=Tgw, Q_sfc=Qs, T_sfc=Ts, T_src=Tsrc,
                         T_eq=wmean(v['T_eq'], length), T_air=wmean(v['T_airSmooth'], length)))
    if n % 500 == 0:
        print(f'  {n}/{len(files)} {t}', flush=True)

h = pd.DataFrame(rows)
h['date'] = h.time.dt.floor('D')
d = h.groupby(['site_no', 'date']).mean(numeric_only=True).reset_index()
d = d.merge(obs[['site_no', 'date', 'tmean']].rename(columns={'tmean': 'T_obs'}), on=['site_no', 'date'], how='left')
d['gw_frac'] = d.Q_gw / (d.Q_gw + d.Q_sfc)
d.to_csv(outcsv, index=False)

d['month'] = d.date.dt.month
cols = ['T_obs', 'T_model', 'T_src', 'T_gw', 'T_sfc', 'gw_frac', 'T_eq', 'T_air', 'Q_gauge']
m = d.dropna(subset=['T_obs'])
pd.set_option('display.width', 200); pd.set_option('display.max_columns', 20)
print('\nAll gauges, days with observations (means):')
print(m.groupby('month')[cols].mean().round(2))
print('\nBy gauge and month:')
names = dict(zip(sites.site_no, sites.station_nm.str[:28]))
m = m.assign(gauge=m.site_no + ' ' + m.site_no.map(names))
print(m.groupby(['gauge', 'month'])[cols].mean().round(2))
