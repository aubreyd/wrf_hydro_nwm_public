"""Seasonal evaluation of simulated stream temperature against USGS daily observations.

Compares daily-mean CHRTOUT riverTemp with USGS daily mean water temperature
(parameter 00010) at gauge reaches, and reports bias / RMSE / r overall, by
season (DJF/MAM/JJA/SON), by gauge group, and per gauge.

Gauge selection and grouping:
  - Only gauges in the HUC8s given by --huc8 are used (default 02040101 and
    02040102, the HUCs fully contained in the DRB model domain); requires a
    'huc8' column in usgs_sites.csv.
  - Each gauge is labelled 'below reservoir' if its reach is downstream of the
    lake reaches of a --lakes reservoir in Route_Link (NHDWaterbodyComID,
    followed via 'to'), else 'above reservoir' (unregulated, baseflow-
    controlled).  Default: DRB Cannonsville (2613174) and Pepacton (1748473).

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

If usgs_daily_temperature.csv has tmax / tmin, a last table compares daily
maximum, minimum and diurnal range (model days in local standard time,
LOCAL_UTC_OFFSET_H, as USGS daily values are) by gauge group and season.

usage: python eval_usgs_seasonal.py OUTDIR OBSDIR START END CACHE.csv [ROUTE_LINK]
                                    [--huc8 H1,H2] [--lakes ID:name,...]
example (DRB spun-up run, evaluation year after spinup):
  python eval_usgs_seasonal.py run_dirs/run.drb_spinup_surft/OUTDIR \
      NASA_WaterTemperature/USGS_OBS_2017-2019 2018-10-01 2019-10-10 sim_surft.csv
example (Lake Lanier domain, Upper Chattahoochee):
  python eval_usgs_seasonal.py run.lanier_spinup/OUTDIR USGS_OBS_LANIER_2017-2019 \
      2018-10-01 2019-10-10 sim_lanier.csv --huc8 03130001 --lakes 6495106:Lanier
"""
import os, glob, argparse
import numpy as np, pandas as pd
import netCDF4 as nc

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for a in ('outdir', 'obsdir', 'start', 'end', 'cache'):
    ap.add_argument(a)
ap.add_argument('route_link', nargs='?')
ap.add_argument('--huc8', default='02040101,02040102',
                help='comma-separated HUC8s to evaluate (default: DRB 02040101,02040102)')
ap.add_argument('--lakes', default='2613174:Cannonsville,1748473:Pepacton',
                help='comma-separated LAKE_ID:name of the reservoirs set up in the model '
                     '(default: DRB Cannonsville, Pepacton)')
args = ap.parse_args()
outdir, obsdir, start, end, cache = args.outdir, args.obsdir, args.start, args.end, args.cache
route_link = args.route_link or os.path.join(os.path.dirname(os.path.abspath(outdir)), 'DOMAIN', 'Route_Link.nc')
RESERVOIRS = {int(k): v for k, v in (x.split(':') for x in args.lakes.split(','))}
EVAL_HUC8 = set(args.huc8.split(','))

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

# --- daily maximum / minimum and diurnal range ---
# USGS daily values are for local standard time days, so the model's daily
# max/min use local-standard-time days (UTC offset below; -5 h = EST for the
# DRB).  The daily-mean statistics above keep UTC days.
LOCAL_UTC_OFFSET_H = -5
if {'tmax', 'tmin'}.issubset(obs.columns):
    loc = sim.shift(-1, freq='h').shift(LOCAL_UTC_OFFSET_H, freq='h')
    dx = loc.resample('D').max().loc[start:end]
    dn = loc.resample('D').min().loc[start:end]
    mm = []
    for name, tab in (('simmax', dx), ('simmin', dn)):
        t = tab.reset_index().melt(id_vars='time', var_name='comid', value_name=name)
        t['comid'] = t.comid.astype(int)
        mm.append(t.rename(columns={'time': 'date'}).set_index(['comid', 'date']))
    x = obs.merge(mm[0].join(mm[1]).reset_index(), on=['comid', 'date'], how='inner')
    x = x.dropna(subset=['tmax', 'tmin', 'simmax', 'simmin'])
    x['season'] = x.date.dt.month.map({12: 'DJF', 1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                                       6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON', 11: 'SON'})
    x['group'] = np.where(x.comid.map(lambda c: c in below), 'below reservoir', 'above reservoir')

    def mxstats(g):
        emax, emin = g.simmax - g.tmax, g.simmin - g.tmin
        return pd.Series({'n': len(g), 'obs_max': g.tmax.mean(), 'bias_max': emax.mean(),
                          'rmse_max': np.sqrt((emax ** 2).mean()), 'bias_min': emin.mean(),
                          'obs_range': (g.tmax - g.tmin).mean(), 'sim_range': (g.simmax - g.simmin).mean()})

    print(f'\nDaily maximum / minimum (local standard time days, UTC{LOCAL_UTC_OFFSET_H:+d} h):')
    print(x.groupby(['group', 'season'])[['tmax', 'tmin', 'simmax', 'simmin']].apply(mxstats)
          .reindex(['DJF', 'MAM', 'JJA', 'SON'], level=1).round(2).to_string())
