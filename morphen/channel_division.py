"""
Instrument-aware sub-band division for (combined) WSClean imaging.

WSClean's default `-channels-out` division splits the *sorted set of unique
channel frequencies* into equal-count chunks, and `-gap-channel-division`
splits at the largest spacing between neighbouring channels. When visibilities
from instruments with different channel widths are concatenated (e.g.
e-MERLIN 1 MHz + VLA 2 MHz), the channels interleave, and the largest
"gaps" after the true inter-band gap are the single-instrument band edges,
so gap division produces 1-channel output images.

This module reads the MS and:

1. builds `-channel-division-frequencies` so that the band is split at the
   real gaps and, inside each contiguous block, evenly over the region where
   all instruments overlap (single-instrument edges are merged into the outer
   sub-bands). Exactly nc - 1 splits are returned, so WSClean assigns one
   output channel per part and does no further division;
2. reports, for each predicted sub-band, how the (unflagged) visibility
   weight is shared between instruments. A large spread of this fraction
   means each sub-band has a different effective uv-coverage/PSF, which biases
   the sub-band fluxes and the fitted spectral index;
3. optionally writes a copy of the MS where WEIGHT_SPECTRUM is rescaled on a
   fine frequency grid so every instrument contributes the same fraction of
   weight at every frequency (the total weight per bin is preserved).

Only casatools and numpy are required.

Command line:
    python channel_division.py MS --nc 8
    python channel_division.py MS --nc 8 --mode weight --json division.json
    python channel_division.py MS --balance-weights OUT.ms --bin-width 32e6
"""
import argparse
import json
import os
import shutil

import numpy as np

try:
    import casatools
except ImportError:  # pragma: no cover
    casatools = None


ROW_CHUNK = 50000


def _table():
    if casatools is None:
        raise ImportError('casatools is required by channel_division.')
    return casatools.table()


def _parallel_hands(npol):
    """Correlation indices of the parallel hands (XX/YY or RR/LL)."""
    if npol == 4:
        return [0, 3]
    if npol == 2:
        return [0, 1]
    return [0]


def _query_string(ddid, obs_ids=None, field=None):
    q = f'DATA_DESC_ID=={ddid} && ANTENNA1!=ANTENNA2'
    if obs_ids is not None:
        q += ' && OBSERVATION_ID IN [' + ','.join(str(o) for o in obs_ids) + ']'
    if field is not None:
        q += ' && FIELD_ID IN [' + ','.join(str(f) for f in np.atleast_1d(field)) + ']'
    return q


def read_channel_table(ms, spws=None, field=None, with_weights=True,
                       group_by='telescope', verbose=True):
    """
    Read one record per (spw, instrument, channel) of a measurement set.

    Parameters
    ----------
    ms : str
        Path to the measurement set.
    spws : list of int, optional
        Restrict to these SPW ids (same meaning as WSClean's -spws).
    field : int or list of int, optional
        Restrict the weight sums to these FIELD_IDs.
    with_weights : bool
        If True, sum the unflagged parallel-hand weights per channel. This
        requires one pass through WEIGHT_SPECTRUM/FLAG.
    group_by : 'telescope' or 'observation'
        How instruments are identified: by OBSERVATION.TELESCOPE_NAME, or by
        OBSERVATION_ID.

    Returns
    -------
    dict with numpy arrays 'freq', 'width', 'spw', 'ddid', 'inst' (int),
    'weight' (NaN when with_weights is False), and the list 'instruments'
    (labels indexed by 'inst') plus 'inst_obs' (observation ids per label).
    """
    ms = ms.rstrip('/')
    tb = _table()

    tb.open(ms + '/SPECTRAL_WINDOW', nomodify=True)
    chan_freq = [tb.getcell('CHAN_FREQ', i) for i in range(tb.nrows())]
    chan_width = [np.abs(tb.getcell('CHAN_WIDTH', i)) for i in range(tb.nrows())]
    tb.close()

    tb.open(ms + '/DATA_DESCRIPTION', nomodify=True)
    dd_spw = tb.getcol('SPECTRAL_WINDOW_ID')
    tb.close()

    tb.open(ms + '/OBSERVATION', nomodify=True)
    tel_names = list(tb.getcol('TELESCOPE_NAME'))
    tb.close()

    tb.open(ms, nomodify=True)
    colnames = tb.colnames()
    ddid_col = tb.getcol('DATA_DESC_ID')
    obs_col = tb.getcol('OBSERVATION_ID')
    has_wspec = ('WEIGHT_SPECTRUM' in colnames and
                 tb.iscelldefined('WEIGHT_SPECTRUM', 0))

    # which observations (and so instruments) actually have rows in each ddid
    pairs = np.unique(np.stack([ddid_col, obs_col]), axis=1).T
    del ddid_col, obs_col

    if group_by == 'observation':
        label_of_obs = {o: f'OBS{o}:{tel_names[o]}' for o in range(len(tel_names))}
    else:
        label_of_obs = {o: tel_names[o] for o in range(len(tel_names))}
    instruments = []
    for o in range(len(tel_names)):
        if label_of_obs[o] not in instruments:
            instruments.append(label_of_obs[o])
    inst_obs = {lab: [o for o in range(len(tel_names)) if label_of_obs[o] == lab]
                for lab in instruments}

    out = {k: [] for k in ('freq', 'width', 'spw', 'ddid', 'inst', 'weight')}
    for ddid in np.unique(pairs[:, 0]):
        spw = int(dd_spw[ddid])
        if spws is not None and spw not in spws:
            continue
        labs = sorted({label_of_obs[o] for d, o in pairs if d == ddid},
                      key=instruments.index)
        for lab in labs:
            nchan = len(chan_freq[spw])
            wsum = np.full(nchan, np.nan)
            if with_weights:
                wsum = np.zeros(nchan)
                sub = tb.query(_query_string(ddid, inst_obs[lab], field),
                               columns='WEIGHT_SPECTRUM,WEIGHT,FLAG')
                nrow = sub.nrows()
                for start in range(0, nrow, ROW_CHUNK):
                    n = min(ROW_CHUNK, nrow - start)
                    flag = sub.getcol('FLAG', start, n)
                    if has_wspec:
                        w = sub.getcol('WEIGHT_SPECTRUM', start, n)
                    else:
                        w = np.broadcast_to(sub.getcol('WEIGHT', start, n)[:, None, :],
                                            flag.shape)
                    hands = _parallel_hands(flag.shape[0])
                    wsum += np.where(flag[hands], 0.0, w[hands]).sum(axis=(0, 2))
                sub.close()
            out['freq'].append(chan_freq[spw])
            out['width'].append(chan_width[spw])
            out['spw'].append(np.full(nchan, spw))
            out['ddid'].append(np.full(nchan, ddid))
            out['inst'].append(np.full(nchan, instruments.index(lab)))
            out['weight'].append(wsum)
    tb.close()

    table = {k: np.concatenate(v) for k, v in out.items()}
    used = np.unique(table['inst'])
    # keep only instruments with data, re-indexed
    remap = {old: new for new, old in enumerate(used)}
    table['inst'] = np.array([remap[i] for i in table['inst']])
    table['instruments'] = [instruments[i] for i in used]
    table['inst_obs'] = {instruments[i]: inst_obs[instruments[i]] for i in used}
    table['has_weight_spectrum'] = has_wspec
    order = np.argsort(table['freq'], kind='stable')
    for k in ('freq', 'width', 'spw', 'ddid', 'inst', 'weight'):
        table[k] = table[k][order]
    if verbose:
        print(f' >> channel_division: {len(table["freq"])} channels '
              f'({len(np.unique(table["freq"]))} unique) from '
              f'{len(table["instruments"])} instrument(s): '
              f'{", ".join(table["instruments"])}')
    return table


def find_blocks(table, gap_factor=3.0):
    """
    Contiguous frequency blocks of the combined coverage.

    A gap is a frequency interval covered by no channel of any instrument and
    wider than gap_factor times the largest channel width. The spacing between
    interleaved channels of different instruments is therefore never a gap.

    Returns a list of dicts with 'lo', 'hi' (channel centres, Hz), 'mask'
    (channel mask), 'common' (lo, hi) where all instruments present overlap,
    and 'insts' (instrument indices present).
    """
    f = table['freq']
    w = table['width']
    lo_edge = f - w / 2
    hi_edge = f + w / 2
    tol = gap_factor * np.max(w)

    # f is sorted; walk the channels and break where the next lower edge
    # starts beyond the running maximum upper edge by more than tol
    starts = [0]
    run_hi = hi_edge[0]
    for i in range(1, len(f)):
        if lo_edge[i] - run_hi > tol:
            starts.append(i)
        run_hi = max(run_hi, hi_edge[i])
    stops = starts[1:] + [len(f)]

    blocks = []
    for a, b in zip(starts, stops):
        mask = np.zeros(len(f), bool)
        mask[a:b] = True
        insts = np.unique(table['inst'][mask])
        c_lo = max(f[mask & (table['inst'] == i)].min() for i in insts)
        c_hi = min(f[mask & (table['inst'] == i)].max() for i in insts)
        if c_hi <= c_lo:  # instruments do not overlap inside this block
            c_lo, c_hi = f[a], f[b - 1]
        blocks.append({'lo': f[a], 'hi': f[b - 1], 'mask': mask,
                       'common': (c_lo, c_hi), 'insts': insts})
    return blocks


def _allocate(nc, sizes):
    """Largest-remainder allocation of nc parts, at least one per block."""
    sizes = np.asarray(sizes, float)
    nb = len(sizes)
    alloc = np.ones(nb, int)
    rest = nc - nb
    if rest > 0:
        share = sizes / sizes.sum() * nc
        extra = np.maximum(np.floor(share) - 1, 0).astype(int)
        # do not over-allocate
        while extra.sum() > rest:
            extra[np.argmax(extra)] -= 1
        alloc += extra
        remainder = share - alloc
        for i in np.argsort(-remainder)[:nc - alloc.sum()]:
            alloc[i] += 1
    return alloc


def _snap(split, ufreq):
    """Move a split to the midpoint between the adjacent unique channels."""
    idx = np.searchsorted(ufreq, split)
    idx = int(np.clip(idx, 1, len(ufreq) - 1))
    return 0.5 * (ufreq[idx - 1] + ufreq[idx])


def compute_channel_division(table, nc, mode='bandwidth', gap_factor=3.0,
                             min_chan=2, verbose=True):
    """
    Compute split frequencies for WSClean's -channel-division-frequencies.

    Parameters
    ----------
    table : dict or str
        Output of read_channel_table, or a path to an MS.
    nc : int
        Number of output channels (-channels-out).
    mode : 'bandwidth' or 'weight'
        Inside each block, split the common region into equal bandwidth or
        into equal total (unflagged) weight.
    gap_factor : float
        See find_blocks.
    min_chan : int
        Warn when a sub-band has fewer channels than this from an instrument
        that is present in its block.

    Returns
    -------
    splits : list of float
        nc - 1 split frequencies in Hz (sorted).
    subbands : list of dict
        Predicted imaging table (see predict_subbands).
    """
    if isinstance(table, str):
        table = read_channel_table(table, with_weights=True, verbose=verbose)
    nc = int(nc)
    ufreq = np.unique(table['freq'])
    if nc <= 1:
        return [], predict_subbands(table, [])

    blocks = find_blocks(table, gap_factor=gap_factor)
    gap_splits = [0.5 * (b0['hi'] + b1['lo']) for b0, b1 in zip(blocks[:-1], blocks[1:])]

    if nc < len(blocks):
        # not enough outputs: split only at the widest gaps
        widths = [b1['lo'] - b0['hi'] for b0, b1 in zip(blocks[:-1], blocks[1:])]
        keep = sorted(np.argsort(widths)[::-1][:nc - 1])
        splits = sorted(gap_splits[i] for i in keep)
        print(f' !! channel_division: nc={nc} is smaller than the number of '
              f'frequency blocks ({len(blocks)}); some sub-bands will span a gap.')
        return splits, predict_subbands(table, splits, verbose=verbose,
                                        min_chan=min_chan)

    sizes = [b['common'][1] - b['common'][0] for b in blocks]
    alloc = _allocate(nc, sizes)

    splits = list(gap_splits)
    for b, k in zip(blocks, alloc):
        if k <= 1:
            continue
        c_lo, c_hi = b['common']
        if mode == 'weight' and np.all(np.isfinite(table['weight'])):
            m = b['mask'] & (table['freq'] >= c_lo) & (table['freq'] <= c_hi)
            fq = table['freq'][m]
            cw = np.cumsum(table['weight'][m])
            if cw[-1] > 0:
                targets = cw[-1] * np.arange(1, k) / k
                inner = [fq[np.searchsorted(cw, t)] for t in targets]
            else:
                inner = list(c_lo + (c_hi - c_lo) * np.arange(1, k) / k)
        else:
            if mode == 'weight':
                print(' !! channel_division: no weights available, '
                      'using equal bandwidth.')
            inner = list(c_lo + (c_hi - c_lo) * np.arange(1, k) / k)
        splits += [_snap(s, ufreq) for s in inner]

    splits = sorted(set(splits))
    if len(splits) != nc - 1:
        print(f' !! channel_division: produced {len(splits)} splits for '
              f'nc={nc} (expected {nc - 1}).')
    subbands = predict_subbands(table, splits, verbose=verbose, min_chan=min_chan,
                                blocks=blocks)
    return splits, subbands


def predict_subbands(table, splits, verbose=True, min_chan=2, blocks=None):
    """
    Predicted sub-bands for a list of split frequencies: frequency range,
    number of unique channels (as in WSClean's imaging table), channels and
    weight fraction per instrument.
    """
    f = table['freq']
    edges = [-np.inf] + list(splits) + [np.inf]
    insts = table['instruments']
    has_w = np.all(np.isfinite(table['weight']))
    if blocks is None:
        blocks = find_blocks(table)
    subbands = []
    for i in range(len(edges) - 1):
        m = (f > edges[i]) & (f < edges[i + 1])
        sb = {'index': i,
              'fmin': float(f[m].min()) if m.any() else np.nan,
              'fmax': float(f[m].max()) if m.any() else np.nan,
              'n_unique': int(len(np.unique(f[m])))}
        wtot = table['weight'][m].sum() if has_w else np.nan
        for j, lab in enumerate(insts):
            mj = m & (table['inst'] == j)
            sb[f'nchan[{lab}]'] = int(mj.sum())
            if has_w:
                sb[f'wfrac[{lab}]'] = (float(table['weight'][mj].sum() / wtot)
                                      if wtot > 0 else np.nan)
        subbands.append(sb)
        # warn about sub-bands missing an instrument that exists in the block
        for b in blocks:
            if m.any() and b['lo'] <= sb['fmin'] <= b['hi']:
                for j in b['insts']:
                    if sb[f'nchan[{insts[j]}]'] < min_chan:
                        print(f' !! channel_division: sub-band {i} '
                              f'({sb["fmin"]/1e6:.1f}-{sb["fmax"]/1e6:.1f} MHz) has '
                              f'only {sb[f"nchan[{insts[j]}]"]} channel(s) of '
                              f'{insts[j]}.')
    if verbose:
        report_division(subbands, insts, has_w)
    return subbands


def report_division(subbands, instruments, has_weights=True, spread_warn=0.1):
    """Print the predicted imaging table and warn on unbalanced weights."""
    head = f'   {"#":>2} {"Freq (MHz)":>17} {"uniq":>5}'
    for lab in instruments:
        head += f' {"n[" + lab[:8] + "]":>12}'
    if has_weights and len(instruments) > 1:
        for lab in instruments:
            head += f' {"w[" + lab[:8] + "]":>12}'
    print(' >> channel_division: predicted imaging table')
    print(head)
    for sb in subbands:
        line = (f'   {sb["index"]:>2} {sb["fmin"]/1e6:8.1f}-{sb["fmax"]/1e6:8.1f} '
                f'{sb["n_unique"]:>5}')
        for lab in instruments:
            line += f' {sb[f"nchan[{lab}]"]:>12}'
        if has_weights and len(instruments) > 1:
            for lab in instruments:
                line += f' {sb[f"wfrac[{lab}]"]:>12.3f}'
        print(line)
    if has_weights and len(instruments) > 1 and len(subbands) > 1:
        for lab in instruments:
            fr = np.array([sb[f'wfrac[{lab}]'] for sb in subbands], float)
            spread = np.nanmax(fr) - np.nanmin(fr)
            if spread > spread_warn:
                print(f' !! channel_division: the {lab} weight fraction varies by '
                      f'{spread:.2f} across sub-bands ({np.nanmin(fr):.2f}-'
                      f'{np.nanmax(fr):.2f}). Sub-band PSFs differ, which can bias '
                      f'sub-band fluxes/spectral index. Consider '
                      f'balance_instrument_weights().')
                break


def format_division_frequencies(splits):
    """Comma-separated list in Hz, for -channel-division-frequencies."""
    return ','.join(f'{s:.6e}' for s in splits)


def auto_division_args(ms, nc, mode='bandwidth', spws=None, field=None,
                       with_weights=True, csv_file=None, verbose=True):
    """
    WSClean arguments for an instrument-aware division of ms into nc outputs.

    Returns the argument string (' -channel-division-frequencies ... ') or an
    empty string when nc <= 1 or no split is needed.
    """
    if int(nc) <= 1:
        return ''
    table = read_channel_table(ms, spws=spws, field=field,
                               with_weights=with_weights or mode == 'weight',
                               verbose=verbose)
    splits, subbands = compute_channel_division(table, nc, mode=mode,
                                                verbose=verbose)
    if csv_file is not None:
        _write_csv(subbands, csv_file)
    if not splits:
        return ''
    return ' -channel-division-frequencies ' + format_division_frequencies(splits) + ' '


def _write_csv(subbands, csv_file):
    keys = list(subbands[0].keys())
    with open(csv_file, 'w') as fh:
        fh.write(','.join(keys) + '\n')
        for sb in subbands:
            fh.write(','.join(str(sb[k]) for k in keys) + '\n')


def _balance_bins(table, blocks, bin_width):
    """Frequency bin edges, bins covering each block fully."""
    edges = []
    for b in blocks:
        c_lo, c_hi = b['common']
        nb = max(1, int(round((c_hi - c_lo) / bin_width)))
        e = list(c_lo + (c_hi - c_lo) * np.arange(nb + 1) / nb)
        # outer bins absorb the single-instrument edges of the block
        e[0] = b['lo'] - 1.0
        e[-1] = b['hi'] + 1.0
        edges.append(np.array(e))
    return edges


def compute_balance_factors(table, bin_width=32e6, target=None, gap_factor=3.0):
    """
    Per-channel weight scale factors that make every instrument contribute
    the same fraction of weight in every frequency bin.

    target : dict {instrument label: fraction}, optional
        Target weight fraction per instrument. Default: the global fraction.

    Returns an array of factors aligned with table['freq'] and the target.
    """
    if not np.all(np.isfinite(table['weight'])):
        raise ValueError('The channel table has no weights '
                         '(use read_channel_table(..., with_weights=True)).')
    insts = table['instruments']
    w = table['weight']
    tot = w.sum()
    if target is None:
        target = {lab: float(w[table['inst'] == j].sum() / tot)
                  for j, lab in enumerate(insts)}
    t = np.array([target[lab] for lab in insts], float)

    factors = np.ones_like(w)
    blocks = find_blocks(table, gap_factor=gap_factor)
    for edges in _balance_bins(table, blocks, bin_width):
        for a, b in zip(edges[:-1], edges[1:]):
            m = (table['freq'] >= a) & (table['freq'] < b)
            wi = np.array([w[m & (table['inst'] == j)].sum() for j in range(len(insts))])
            present = wi > 0
            if present.sum() < 2:
                continue
            ti = np.where(present, t, 0.0)
            ti = ti / ti.sum()
            wbin = wi.sum()
            for j in np.where(present)[0]:
                factors[m & (table['inst'] == j)] = ti[j] * wbin / wi[j]
    return factors, target


def balance_instrument_weights(ms, out_ms, bin_width=32e6, target=None,
                               spws=None, field=None, overwrite=False,
                               verbose=True):
    """
    Write a copy of ms whose WEIGHT_SPECTRUM is rescaled per instrument so that
    every instrument contributes the same fraction of weight in every
    frequency bin of width bin_width (Hz). The input MS is not modified.

    The total weight in each bin is preserved, so the relative weighting of
    WSClean's sub-bands in the MFS image is unchanged. WEIGHT is scaled by the
    mean factor of each SPW, and SIGMA/SIGMA_SPECTRUM by 1/sqrt(factor).

    Returns the path of the new MS.
    """
    ms = ms.rstrip('/')
    out_ms = out_ms.rstrip('/')
    if os.path.abspath(ms) == os.path.abspath(out_ms):
        raise ValueError('out_ms must differ from ms.')
    if os.path.exists(out_ms):
        if not overwrite:
            raise FileExistsError(f'{out_ms} exists (use overwrite=True).')
        shutil.rmtree(out_ms)

    table = read_channel_table(ms, spws=spws, field=field, with_weights=True,
                               verbose=verbose)
    if not table['has_weight_spectrum']:
        raise ValueError('The MS has no WEIGHT_SPECTRUM; create it first, e.g. '
                         "initweights(vis, wtmode='weight', dowtsp=True).")
    factors, target = compute_balance_factors(table, bin_width=bin_width,
                                              target=target)
    if verbose:
        print(' >> channel_division: target weight fractions: ' +
              ', '.join(f'{k}={v:.3f}' for k, v in target.items()))
        print(f' >> channel_division: scale factors range '
              f'{factors.min():.3f}-{factors.max():.3f}')

    print(f' >> channel_division: copying {ms} -> {out_ms}')
    shutil.copytree(ms, out_ms, symlinks=True)

    tb = _table()
    tb.open(out_ms, nomodify=False)
    colnames = tb.colnames()
    has_sspec = ('SIGMA_SPECTRUM' in colnames and
                 tb.iscelldefined('SIGMA_SPECTRUM', 0))
    insts = table['instruments']
    for ddid in np.unique(table['ddid']):
        for j, lab in enumerate(insts):
            m = (table['ddid'] == ddid) & (table['inst'] == j)
            if not m.any():
                continue
            # channel order inside the SPW is the original CHAN_FREQ order
            spw = int(table['spw'][m][0])
            fac = _factors_in_spw_order(ms, spw, table['freq'][m], factors[m])
            if np.allclose(fac, 1.0):
                continue
            # the query includes auto-correlations and every field: the
            # scaling is a property of the instrument/frequency, not of the
            # row selection used for the weight sums
            sub = tb.query(f'DATA_DESC_ID=={ddid} && OBSERVATION_ID IN ['
                           + ','.join(str(o) for o in table['inst_obs'][lab]) + ']')
            nrow = sub.nrows()
            mean_fac = float(fac.mean())
            for start in range(0, nrow, ROW_CHUNK):
                n = min(ROW_CHUNK, nrow - start)
                ws = sub.getcol('WEIGHT_SPECTRUM', start, n)
                sub.putcol('WEIGHT_SPECTRUM', ws * fac[None, :, None], start, n)
                wt = sub.getcol('WEIGHT', start, n)
                sub.putcol('WEIGHT', wt * mean_fac, start, n)
                sg = sub.getcol('SIGMA', start, n)
                sub.putcol('SIGMA', sg / np.sqrt(mean_fac), start, n)
                if has_sspec:
                    ss = sub.getcol('SIGMA_SPECTRUM', start, n)
                    sub.putcol('SIGMA_SPECTRUM', ss / np.sqrt(fac)[None, :, None],
                               start, n)
            sub.close()
            if verbose:
                print(f'    ddid {ddid:2d} spw {spw:2d} {lab:>10}: {nrow} rows, '
                      f'factor {fac.min():.3f}-{fac.max():.3f}')
    tb.flush()
    tb.close()

    with open(out_ms + '_balance_factors.json', 'w') as fh:
        json.dump({'source_ms': ms, 'bin_width': bin_width, 'target': target,
                   'freq': table['freq'].tolist(),
                   'inst': [insts[i] for i in table['inst']],
                   'factor': factors.tolist()}, fh)
    return out_ms


def _factors_in_spw_order(ms, spw, freqs_sorted, fac_sorted):
    """Map factors from frequency-sorted order to the SPW's channel order."""
    tb = _table()
    tb.open(ms + '/SPECTRAL_WINDOW', nomodify=True)
    cf = tb.getcell('CHAN_FREQ', spw)
    tb.close()
    lookup = dict(zip(np.round(freqs_sorted, 3), fac_sorted))
    return np.array([lookup[v] for v in np.round(cf, 3)])


def parse_spws(opt_args):
    """SPW list from a '-spws 0,1,2' WSClean argument string, or None."""
    toks = opt_args.split()
    if '-spws' in toks:
        return [int(s) for s in toks[toks.index('-spws') + 1].split(',')]
    return None


def parse_field(opt_args):
    """Field id(s) from a '-field 0' WSClean argument string, or None."""
    toks = opt_args.split()
    if '-field' in toks:
        val = toks[toks.index('-field') + 1]
        if val.lower() == 'all':
            return None
        return [int(s) for s in val.split(',')]
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Instrument-aware -channel-division-frequencies for WSClean, '
                    'and optional per-instrument weight balancing.')
    parser.add_argument('ms', help='Measurement set.')
    parser.add_argument('--nc', type=int, default=4, help='-channels-out.')
    parser.add_argument('--mode', default='bandwidth', choices=['bandwidth', 'weight'])
    parser.add_argument('--spws', default=None, help='Comma-separated SPW ids.')
    parser.add_argument('--field', default=None, help='Comma-separated field ids.')
    parser.add_argument('--no-weights', action='store_true',
                        help='Do not read weights (faster, no balance report).')
    parser.add_argument('--group-by', default='telescope',
                        choices=['telescope', 'observation'])
    parser.add_argument('--json', default=None, help='Write splits and table to JSON.')
    parser.add_argument('--balance-weights', default=None, metavar='OUT_MS',
                        help='Write a weight-balanced copy of the MS.')
    parser.add_argument('--bin-width', type=float, default=32e6,
                        help='Balancing bin width in Hz (default 32e6).')
    parser.add_argument('--overwrite', action='store_true')
    a = parser.parse_args()

    spws = [int(s) for s in a.spws.split(',')] if a.spws else None
    field = [int(s) for s in a.field.split(',')] if a.field else None

    if a.balance_weights:
        out = balance_instrument_weights(a.ms, a.balance_weights,
                                         bin_width=a.bin_width, spws=spws,
                                         field=field, overwrite=a.overwrite)
        print(' >> Balanced copy; division of the new MS:')
        a.ms = out

    tab = read_channel_table(a.ms, spws=spws, field=field,
                             with_weights=not a.no_weights, group_by=a.group_by)
    splits, subbands = compute_channel_division(tab, a.nc, mode=a.mode)
    arg = format_division_frequencies(splits)
    print(f' >> -channels-out {a.nc} -channel-division-frequencies {arg}')
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump({'ms': a.ms, 'nc': a.nc, 'mode': a.mode, 'splits': splits,
                       'subbands': subbands}, fh, indent=1, default=float)
