"""Compare WRF-Hydro runs that differ only in MPI rank count.

For each run against the reference run, reports per variable how many files
differ and the largest absolute difference, for the main outputs and for every
numeric variable in the RESTART / HYDRO_RST files both runs wrote.  Results
that depend on the domain decomposition show up as differences; float
rounding from reductions summed in a different order shows up as differences
of a few ulps that may then grow.

usage: python check_mpi_ranks.py REF_RUN RUN [RUN ...]
   RUN: run directories with OUTDIR/
example (DRB):
   python check_mpi_ranks.py run.lkh_new run.mpi_n1 run.mpi_n4 run.mpi_n8 run.mpi_n32
"""
import sys, os, glob
import numpy as np
import netCDF4 as nc

OUT_VARS = [('*.CHRTOUT_DOMAIN1', ['streamflow', 'riverTemp']),
            ('*.LAKEOUT_DOMAIN1', ['water_sfc_elev', 'outflow']),
            ('*.GWOUT_DOMAIN1', ['depth', 'outflow']),
            ('*.LDASOUT_DOMAIN1', ['TG', 'SOIL_T', 'SNEQV', 'LH', 'LAKE_T3D', 'LAKE_DEPTH2D'])]


def arr(v):
    return np.ma.filled(np.asarray(v[:], dtype=float), np.nan)


def maxdiff(x, y):
    if x.shape != y.shape:
        return np.inf
    k = np.isfinite(x) & np.isfinite(y)
    nan_mismatch = np.any(np.isfinite(x) != np.isfinite(y))
    d = float(np.max(np.abs(x[k] - y[k]))) if k.any() else 0.0
    return np.inf if nan_mismatch else d


def compare(ref, run):
    print(f'== {run} vs {ref}')
    ok = True
    for pat, vs in OUT_VARS:
        stats = {v: [0, 0, 0.0, None] for v in vs}
        for g in sorted(glob.glob(os.path.join(run, 'OUTDIR', pat))):
            f = os.path.join(ref, 'OUTDIR', os.path.basename(g))
            if not os.path.exists(f):
                continue
            a, b = nc.Dataset(f), nc.Dataset(g)
            for v in vs:
                if v in a.variables and v in b.variables:
                    s = stats[v]
                    s[0] += 1
                    d = maxdiff(arr(a[v]), arr(b[v]))
                    if d > 0:
                        s[1] += 1
                        if s[3] is None:
                            s[3] = os.path.basename(g)[:12]
                        s[2] = max(s[2], d)
        for v, (n, nd, d, first) in stats.items():
            if n:
                tag = 'identical' if nd == 0 else f'{nd} differ, max |diff| {d:.3g}, first {first}'
                print(f'  {v:15s} {n:4d} files: {tag}')
                ok &= nd == 0
    for pat in ('RESTART.*_DOMAIN1', 'HYDRO_RST.*_DOMAIN1'):
        for g in sorted(glob.glob(os.path.join(run, 'OUTDIR', pat))):
            f = os.path.join(ref, 'OUTDIR', os.path.basename(g))
            if not os.path.exists(f):
                continue
            a, b = nc.Dataset(f), nc.Dataset(g)
            diffs = []
            for v in a.variables:
                if v in b.variables and a[v].dtype.kind in 'fi':
                    d = maxdiff(arr(a[v]), arr(b[v]))
                    if d > 0:
                        diffs.append(f'{v}({d:.2g})')
            print(f'  {os.path.basename(g)}: ' + (', '.join(diffs) if diffs else 'identical'))
            ok &= not diffs
    print('  RANK INDEPENDENCE:', 'bit-identical' if ok else 'DIFFERENCES')
    return ok


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ref = sys.argv[1]
    ok = True
    for run in sys.argv[2:]:
        ok &= compare(ref, run)
    sys.exit(0 if ok else 1)
