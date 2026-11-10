"""Add stream/bucket temperature parameters to copies of Route_Link and GWBUCKPARM.

Writes the optional per-reach Route_Link variables read by
init_stream_temp_params (module_stream_temperature_eq.F)
    Kexch    bulk surface heat exchange coefficient [W m-2 C-1]  (model default 50)
    Ksw      shortwave equilibrium coefficient      [W m-2 C-1]  (model default 400)
    TauTair  air temperature smoothing e-folding time [s]        (model default 86400)
    Shade    shade fraction of incoming shortwave   [0-1]        (model default 0)
    Albedo   water surface albedo                   [0-1]        (model default 0)
and the optional per-bucket GWBUCKPARM variable read by readBucket_nhd
    Dtherm   bucket thermal damping depth [m water equivalent]   (model default 0)

Any parameter not given is not written (the model then uses its default).

usage:
  python add_temp_params.py ROUTE_LINK_IN ROUTE_LINK_OUT GWBUCKPARM_IN GWBUCKPARM_OUT \
         [--Kexch V] [--Ksw V] [--TauTair V] [--Shade V] [--Albedo V] [--Dtherm V] [--test-pattern]

  V is a uniform value for every reach (bucket).  --test-pattern writes a
  spatially varying pattern by stream order for all six parameters, with one
  invalid value in each of Kexch, Ksw, Shade and Dtherm (to exercise the
  model's replacement of invalid values by the default); it is meant for code
  testing, not for physics.
"""
import argparse, shutil
import numpy as np
import netCDF4 as nc

UNITS = {'Kexch': 'W m-2 C-1', 'Ksw': 'W m-2 C-1', 'TauTair': 's', 'Shade': '1', 'Albedo': '1',
         'Dtherm': 'm'}
LONG = {'Kexch': 'stream temperature bulk surface heat exchange coefficient',
        'Ksw': 'stream temperature shortwave equilibrium coefficient',
        'TauTair': 'e-folding time of the air temperature smoothing for stream temperature',
        'Shade': 'shade fraction of incoming shortwave on the reach',
        'Albedo': 'water surface albedo of the reach',
        'Dtherm': 'groundwater bucket thermal damping depth (water equivalent)'}


def add_var(path, name, dim, values):
    with nc.Dataset(path, 'a') as d:
        if name in d.variables:
            d[name][:] = values
        else:
            v = d.createVariable(name, 'f4', (dim,))
            v[:] = values
            v.units = UNITS[name]
            v.long_name = LONG[name]
    print(f'  {path}: {name} min {np.nanmin(values):.4g} max {np.nanmax(values):.4g}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rl_in'); ap.add_argument('rl_out'); ap.add_argument('gw_in'); ap.add_argument('gw_out')
    for p in ('Kexch', 'Ksw', 'TauTair', 'Shade', 'Albedo', 'Dtherm'):
        ap.add_argument('--' + p, type=float)
    ap.add_argument('--test-pattern', action='store_true')
    a = ap.parse_args()

    shutil.copyfile(a.rl_in, a.rl_out)
    shutil.copyfile(a.gw_in, a.gw_out)
    with nc.Dataset(a.rl_out) as d:
        rl_dim = d['link'].dimensions[0]
        order = d['order'][:].astype(int)
    with nc.Dataset(a.gw_out) as d:
        gw_dim = d['ComID'].dimensions[0]
        nb = len(d['ComID'][:])
    nr = len(order)

    vals = {}
    if a.test_pattern:
        vals['Kexch'] = np.where(order >= 4, 75.0, 50.0)
        vals['Ksw'] = np.where(order == 1, 300.0, 400.0)
        vals['TauTair'] = np.where(order >= 5, 172800.0, 86400.0)
        vals['Shade'] = np.where(order <= 2, 0.5, 0.1)
        vals['Albedo'] = np.full(nr, 0.08)
        vals['Dtherm'] = np.where(np.arange(nb) % 2 == 0, 1.0, 0.0)
        vals['Kexch'][0] = -1.0        # invalid: <= 0
        vals['Ksw'][1] = np.nan        # invalid: NaN
        vals['Shade'][2] = 1.5         # invalid: > 1
        vals['Dtherm'][3] = -5.0       # invalid: < 0
    for p in ('Kexch', 'Ksw', 'TauTair', 'Shade', 'Albedo'):
        if getattr(a, p) is not None:
            vals[p] = np.full(nr, getattr(a, p))
    if a.Dtherm is not None:
        vals['Dtherm'] = np.full(nb, a.Dtherm)

    for p, v in vals.items():
        if p == 'Dtherm':
            add_var(a.gw_out, p, gw_dim, v)
        else:
            add_var(a.rl_out, p, rl_dim, v)


if __name__ == '__main__':
    main()
