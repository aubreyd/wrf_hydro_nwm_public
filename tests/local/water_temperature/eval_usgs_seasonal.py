"""Seasonal evaluation of simulated stream temperature against USGS daily observations.

Compares daily-mean CHRTOUT riverTemp with USGS daily mean water temperature
(parameter 00010) at gauge reaches, and reports bias / RMSE / r overall, by
season (DJF/MAM/JJA/SON), by gauge group, and per gauge.

Gauge selection and grouping:
  - Only gauges in HUC8 02040101 and 02040102 are used (the HUCs fully
    contained in the DRB model domain); requires a 'huc8' column in
    usgs_sites.csv.
  - Each gauge is labelled 'below reservoir' if its reach is downstream of the
    Cannonsville (2613174) or Pepacton (1748473) lake reaches in Route_Link
    (NHDWaterbodyComID, followed via 'to'), else 'above reservoir'
    (unregulated, baseflow-controlled).

Inputs:
  OUTDIR      model output directory with hourly *.CHRTOUT_DOMAIN1 files
              (not read if CACHE.csv already exists; pass any string then)
  OBSDIR      directory with usgs_sites.csv (site_no, comid, huc8, ...) and
              usgs_daily_temperature.csv (site_no, date, tmean, comid, ...),
              as written by DRB_yates/fetch_usgs_temperature.py plus huc8
  START, END  evaluation period (YYYY-MM-DD, inclusive)
  CACHE.csv   riverTemp at gauge reaches, hourly; written on first use
              (extraction of ~14k hourly files takes ~6 min), read thereafter
  ROUTE_LINK  optional; default OUTDIR/../DOMAIN/Route_Link.nc

Hourly CHRTOUT times are end-of-hour; they are shifted back 1 h before daily
averaging so 01..24 UTC fall on the same day.

usage: python eval_usgs_seasonal.py OUTDIR OBSDIR START END CACHE.csv [ROUTE_LINK]
example (DRB spun-up run, evaluation year after spinup):
  python eval_usgs_seasonal.py run_dirs/run.drb_spinup_surft/OUTDIR \
      NASA_WaterTemperature/USGS_OBS_2017-2019 2018-10-01 2019-10-10 sim_surft.csv
"""
import sys, os, glob
import numpy as np, pandas as pd
import netCDF4 as nc

outdir, obsdir, start, end, cache = sys.argv[1:6]
route_link = sys.argv[6] if len(sys.argv) > 6 else os.path.join(os.path.dirname(os.path.abspath(outdir)), 'DOMAIN', 'Route_Link.nc')
RESERVOIRS = {2613174: 'Cannonsville', 1748473: 'Pepacton'}   # the reservoirs set up in the model
EVAL_HUC8 = {'02040101', '02040102'}   # HUCs fully contained in the DRB model domain

sites = pd.read_csv(f'{obsdir}/usgs_sites.csv', dtype={'site_no': str, 'huc8': str})
sites = sites[sites.huc8.isin(EVAL_HUC8)]
obs = pd.read_csv(f'{obsdir}/usgs_daily_temperature.csv', dtype={'site_no': str}, parse_dates=['date'])
obs = obs[obs.site_no.isin(sites.site_no) & obs.tmean.notna()]

if os.path.exists(cache):
    sim = pd.read_csv(cache, parse_dates=['time'])
else:
    files = sorted(glob.glob(f'{outdir}/*.CHRTOUT_DOMAIN1'))
    with nc.Dataset(files[0]) as d:
        fid = d['feature_id'][:]
    idx = {int(c): int(np.where(fid == c)[0][0]) for c in sites.comid if (fid == c).any()}
    cols = list(idx)
    pos = [idx[c] for c in cols]
    rows, times = [], []
    for f in files:
        with nc.Dataset(f) as d:
            rows.append(d['riverTemp'][:][pos])
        b = os.path.basename(f)[:12]
        times.append(pd.Timestamp(f'{b[:4]}-{b[4:6]}-{b[6:8]} {b[8:10]}:{b[10:12]}'))
    sim = pd.DataFrame(np.array(rows), columns=[str(c) for c in cols])
    sim.insert(0, 'time', times)
    sim.to_csv(cache, index=False)

sim = sim.set_index('time')
# hourly stamps are end-of-hour; shift 1 h so 01..24 UTC fall on the same day
daily = sim.shift(-1, freq='h').resample('D').mean()
daily = daily.loc[start:end]
long = daily.reset_index().melt(id_vars='time', var_name='comid', value_name='sim')
long['comid'] = long.comid.astype(int)
long = long.rename(columns={'time': 'date'})
m = obs.merge(long, on=['comid', 'date'], how='inner').dropna(subset=['sim'])
m['err'] = m.sim - m.tmean
m['season'] = m.date.dt.month.map({12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                                    6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'})

# --- classify gauges: below a reservoir if downstream of its lake reaches in Route_Link ---
with nc.Dataset(route_link) as rl:
    link = rl['link'][:].astype('int64'); to = rl['to'][:].astype('int64'); wb = rl['NHDWaterbodyComID'][:].astype('int64')
down = dict(zip(link.tolist(), to.tolist()))
below = {}
for lake, name in RESERVOIRS.items():
    for l in link[wb == lake].tolist():
        x = down.get(l, 0)
        while x and x in down:
            if wb[link == x][0] != lake:
                below.setdefault(x, name)
            x = down[x]
m['group'] = np.where(m.comid.map(lambda c: c in below), 'below reservoir', 'above reservoir')
m['reservoir'] = m.comid.map(lambda c: below.get(c, '-'))

def stats(g):
    return pd.Series({'n': len(g), 'obs': g.tmean.mean(), 'sim': g.sim.mean(), 'bias': g.err.mean(),
                      'rmse': np.sqrt((g.err ** 2).mean()),
                      'r': np.corrcoef(g.sim, g.tmean)[0, 1] if len(g) > 2 else np.nan})

print(f'{start} .. {end}: {m.site_no.nunique()} sites, {len(m)} site-days')
print('\nALL:', stats(m).round(2).to_dict())
print('\nBy season:')
print(m.groupby('season').apply(stats).loc[['DJF', 'MAM', 'JJA', 'SON']].round(2).to_string())
print('\nBy group:')
print(m.groupby('group').apply(stats).round(2).to_string())
print('\nBy group and season:')
print(m.groupby(['group', 'season']).apply(stats).reindex(['DJF', 'MAM', 'JJA', 'SON'], level=1).round(2).to_string())
print('\nBy site (bias by season, overall RMSE):')
piv = m.pivot_table(index=['group', 'reservoir', 'site_no', 'station_nm'], columns='season', values='err', aggfunc='mean')
piv = piv.reindex(columns=['DJF', 'MAM', 'JJA', 'SON'])
piv['RMSE'] = m.groupby(['group', 'reservoir', 'site_no', 'station_nm']).err.apply(lambda e: np.sqrt((e ** 2).mean()))
piv['n'] = m.groupby(['group', 'reservoir', 'site_no', 'station_nm']).size()
print(piv.round(2).to_string())
