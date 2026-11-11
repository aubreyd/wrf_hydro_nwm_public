"""Hydrograph + stream temperature plots at USGS gauges: observed vs simulated daily means.

Inputs as for eval_usgs_seasonal.py (OBSDIR with usgs_sites.csv incl. huc8,
usgs_daily_temperature.csv and usgs_daily_discharge.csv).

usage: python plot_gauges.py OUTDIR OBSDIR HUC8S CACHE.csv FIGDIR TITLE [START]
  HUC8S       comma-separated HUC8s of the gauges to plot
  CACHE.csv   hourly streamflow + riverTemp at gauge reaches (written on first use)
  START       first day to plot and evaluate (YYYY-MM-DD), e.g. the end of spin-up;
              default: the start of the model output
One PNG per gauge (two panels, shared time axis) and gauge_stats.csv.
"""
import sys, os, glob
import numpy as np, pandas as pd, netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

outdir, obsdir, hucs, cache, figdir, title = sys.argv[1:7]
start = pd.Timestamp(sys.argv[7]) if len(sys.argv) > 7 else None
os.makedirs(figdir, exist_ok=True)

sites = pd.read_csv(f'{obsdir}/usgs_sites.csv', dtype={'site_no': str, 'huc8': str})
sites = sites[sites.huc8.isin(hucs.split(','))].sort_values('site_no')
tobs = pd.read_csv(f'{obsdir}/usgs_daily_temperature.csv', dtype={'site_no': str}, parse_dates=['date'])
qobs = pd.read_csv(f'{obsdir}/usgs_daily_discharge.csv', dtype={'site_no': str}, parse_dates=['date'])

comids = sorted(set(int(c) for c in sites.comid))
if os.path.exists(cache):
    h = pd.read_csv(cache, parse_dates=['time'])
else:
    files = sorted(glob.glob(f'{outdir}/*.CHRTOUT_DOMAIN1'))
    fid = nc.Dataset(files[0])['feature_id'][:].astype(int)
    pos = [int(np.where(fid == c)[0][0]) for c in comids]
    rows = []
    for f in files:
        with nc.Dataset(f) as d:
            q = np.ma.filled(d['streamflow'][:].astype(float), np.nan)[pos]
            t = np.ma.filled(d['riverTemp'][:].astype(float), np.nan)[pos]
        b = os.path.basename(f)[:12]
        rows.append([pd.Timestamp(f'{b[:4]}-{b[4:6]}-{b[6:8]} {b[8:10]}:00')] + list(q) + list(t))
    h = pd.DataFrame(rows, columns=['time'] + [f'Q_{c}' for c in comids] + [f'T_{c}' for c in comids])
    h.to_csv(cache, index=False)
# hourly stamps are end-of-hour; shift 1 h so 01..24 UTC fall on the same day (as eval_usgs_seasonal.py)
daily = h.set_index('time').shift(-1, freq='h').resample('D').mean()
if start is not None:
    daily = daily[daily.index >= start]
t0, t1 = daily.index.min(), daily.index.max()

OBS, SIM = '#3d3d3d', '#2a78d6'          # observed: neutral ink; simulated: one blue (entity colour, fixed)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6})


def panel(ax, o, s, ylabel, zero_floor=False):
    if o is not None and len(o):
        ax.plot(o.index, o.values, color=OBS, lw=1.2, label='Observed (USGS)', zorder=3)
    ax.plot(s.index, s.values, color=SIM, lw=1.2, label='WRF-Hydro', zorder=2)
    ax.set_ylabel(ylabel)
    if zero_floor:
        ax.set_ylim(bottom=0)
    ax.set_xlim(t0, t1)


def stats(o, s):
    j = pd.concat([o, s], axis=1, join='inner').dropna()
    if len(j) < 10:
        return ''
    e = j.iloc[:, 1] - j.iloc[:, 0]
    return f'bias {e.mean():+.2f}, RMSE {np.sqrt((e**2).mean()):.2f}, r {j.corr().iloc[0, 1]:.2f} (n={len(j)})'


summary = []
for _, r in sites.iterrows():
    s_no, c = r.site_no, int(r.comid)
    to = tobs[tobs.site_no == s_no].set_index('date').tmean.dropna()
    qo = qobs[qobs.site_no == s_no].set_index('date').q_cms.dropna()
    to = to[(to.index >= t0) & (to.index <= t1)]; qo = qo[(qo.index >= t0) & (qo.index <= t1)]
    if len(to) == 0 and len(qo) == 0:
        continue
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 5.6), sharex=True, gridspec_kw={'hspace': 0.2})
    panel(a1, qo, daily[f'Q_{c}'], 'Streamflow (m³/s)', zero_floor=True)
    panel(a2, to, daily[f'T_{c}'], 'Water temperature (°C)')
    qs = stats(qo, daily[f'Q_{c}'])
    ts = stats(to, daily[f'T_{c}'])
    fig.suptitle(f'{s_no}  {r.station_nm.title()}  (reach {c})', x=0.07, y=0.995, ha='left',
                 fontsize=10, fontweight='bold')
    a1.set_title(('Q (m³/s): ' + qs) if qs else 'no observed streamflow', loc='right', fontsize=8, color='#444444')
    a2.set_title(('T (°C): ' + ts) if ts else 'no observed temperature', loc='right', fontsize=8, color='#444444')
    h_, l_ = a2.get_legend_handles_labels()
    fig.legend(h_, l_, loc='lower left', bbox_to_anchor=(0.08, -0.005), frameon=False, fontsize=8, ncol=2)
    a2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    a2.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    fig.text(0.99, 0.005, title + f'; daily means {t0:%Y-%m-%d} to {t1:%Y-%m-%d}',
             ha='right', va='bottom', fontsize=7, color='#777777')
    fn = os.path.join(figdir, f'gauge_{s_no}.png')
    fig.savefig(fn, dpi=150, bbox_inches='tight'); plt.close(fig)
    summary.append((s_no, r.station_nm, fn, qs, ts))
    print(fn)

pd.DataFrame(summary, columns=['site_no', 'station_nm', 'figure', 'Q_stats', 'T_stats']).to_csv(
    os.path.join(figdir, 'gauge_stats.csv'), index=False)
