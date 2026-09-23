"""
 ____  _     _
|  _ \(_)___| |_ __ _ _ __   ___ ___  ___
| | | | / __| __/ _` | '_ \ / __/ _ \/ __|
| |_| | \__ \ || (_| | | | | (_|  __/\__ \
|____/|_|___/\__\__,_|_| |_|\___\___||___/
  ____                                     _                         
 / ___|   ___    ___   _ __ ___     ___   | |   ___     __ _   _   _ 
| |      / _ \  / __| | '_ ` _ \   / _ \  | |  / _ \   / _` | | | | |
| |___  | (_) | \__ \ | | | | | | | (_) | | | | (_) | | (_| | | |_| |
 \____|  \___/  |___/ |_| |_| |_|  \___/  |_|  \___/   \__, |  \__, |
                                                       |___/   |___/ 

"""



def luminosity_distance_cosmo(z,Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    D_L = cosmo.luminosity_distance(z).value
    # print('D_l = ', D_L)  # 946.9318492873492 Mpc
    return(D_L)

def comoving_distance_cosmo(z,Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    D_C = cosmo.comoving_distance(z).value
    # print('D_l = ', D_L)  # 946.9318492873492 Mpc
    return(D_C)

def angular_distance_cosmo(z, Om0=0.308):
    h = 67.8# * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    return(d_A)

def arcsec_to_pc(z, cell_size=1, Om0=0.308):
    h = 67.8# * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    theta = 1 * u.arcsec
    distance_pc = (theta * d_A).to(u.pc, u.dimensionless_angles())
    # unit is Mpc only now
    # print('Linear Distance = ', distance_pc)  # 3.384745689510495 Mpc
    return (distance_pc)

    
    

def pc_to_arcsec(parsecs, redshift, h=67.8, Om0=0.308):
    """
    Convert a distance in parsecs to an angular size in arcseconds,
    based on the angular diameter distance to an object at redshift z.
    
    Parameters:
    parsecs (float): Distance in parsecs
    redshift (float): Redshift of the object
    h (float): Hubble constant in km/s/Mpc
    Om0 (float): Matter density parameter
    
    Returns:
    float: Angular size in arcseconds
    """
    # Create the cosmology object
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    
    # Calculate the angular diameter distance in Mpc
    d_A = cosmo.angular_diameter_distance(redshift)
    
    # Convert parsecs to the same units as d_A (Mpc)
    distance_mpc = parsecs * u.pc.to(u.Mpc)
    
    # Calculate the angle in radians: theta = physical_size / d_A
    theta_rad = distance_mpc / d_A.value
    
    # Convert radians to arcseconds
    theta_arcsec = theta_rad * u.rad.to(u.arcsec)
    
    return theta_arcsec

def pixsize_to_pc(z, cell_size, Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    theta = cell_size * u.arcsec
    distance_pc = (theta * d_A).to(u.pc, u.dimensionless_angles())  # unit is Mpc only now

    # print('Linear Distance = ', distance_pc)  # 3.384745689510495 Mpc
    return (distance_pc.value)



def cosmo_stats(imagename,z,results=None):
    """
    Get beam shape info in physical units.
    """
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(imagename)
    pc_scale, bmaj, bmin, BA_pc = beam_physical_area(imagename, z=z)
    results['arcsec_to_pc'] = pc_scale.value
    results['bmaj_pc'] = bmaj.value
    results['bmin_pc'] = bmin.value
    results['BA_pc'] = BA_pc.value
    return(results)

def find_z_NED_old(source_name, return_luminosity_distance=False):
    """
    Find the redshift of a source (by name) from NED.

    Parameters
    ----------
    source_name : str
        Source name.
            Example: 'VV705'
    return_luminosity_distance : bool, optional
        If True, also calculate and return the luminosity distance.
        Default: False

    Returns
    -------
    redshift_NED : float, None
        The redshift value from NED
    luminosity_distance : astropy.units.Quantity, optional
        The luminosity distance (only returned if return_luminosity_distance=True)
    """
    from astroquery.ipac.ned import Ned
    
    result_table = Ned.query_object(source_name)
    redshift_NED = result_table['Redshift'].data.data
    
    if redshift_NED.shape[0] == 0:
        return None
    else:
        z = redshift_NED[0]
        
        if return_luminosity_distance:
            from astropy.cosmology import Planck18
            
            # Calculate luminosity distance
            lum_dist = Planck18.luminosity_distance(z)
            
            # Print the luminosity distance
            print(f"Luminosity distance for {source_name} (z={z:.4f}): {lum_dist:.2f}")
            
            return z, lum_dist
        else:
            return z


# Redshift lookup
# ---------------
# `find_z` is the entry point; `find_z_NED` is kept as a thin backwards-compatible
# alias. The lookup is a chain of independent services, tried in order and each
# with its own short timeout, so a single degraded service cannot stall a batch:
#
#   0. local cache / manual override        instant, and fully reproducible
#   1. NED ObjectLookup  (REST/JSON)        ~0.5 s   name -> z, dz, refcode, position
#   2. SIMBAD TAP        (ADQL on `ident`)  ~0.1 s   name -> z
#   3. Sesame (CDS) name resolver           ~1 s     name -> position
#      then a cone search around that position:
#        NED TAP  (NEDTAP.objdir)           ~1 s
#        SIMBAD TAP (basic)                 ~0.1 s
#   4. astroquery `Ned.query_object`        15-60 s  legacy CGI, circuit-broken
#
# Benchmarked on 28 GOALS/local-merger names (2026-09-08, see
# test_redshift_services.py and redshift_benchmark.log):
#   NED ObjectLookup  27/28  avg 0.90 s  <- the workhorse; NED's own name resolver
#                                           (NED-NNS) takes unspaced names such as
#                                           'CGCG436-030' directly, and returns dz
#                                           and the bibcode as well
#   SIMBAD TAP ident  28/28  avg 0.22 s  <- fastest; needs the right catalogue
#                                           prefix, which is why _normalize_source_name
#                                           maps CGCG -> 'Z' (SIMBAD files
#                                           'CGCG 436-030' only as 'Z 436-30')
#   Sesame + cone     28/28  avg 0.3-1.6 s  <- name-format independent; resolves
#                                           what no name lookup does, e.g.
#                                           'IRAS F17132+5313'
#   astroquery NED    28/28  avg 3.18 s  <- the old sole path, and the whole
#                                           problem: it scrapes the legacy
#                                           objsearch CGI, so it is fine in one
#                                           run and times out at 20-60 s in the
#                                           next (6/15 on 2026-05-12). Kept last.
#
# Why the values can differ slightly between services: NED and SIMBAD quote
# redshifts from different references (e.g. Arp 220: NED 0.018398, SIMBAD 0.018100
# -> ~90 km/s, a ~1.6% difference in luminosity distance). `prefer` therefore
# fixes a deterministic service order rather than racing them, and every resolved
# value is cached with the service and refcode it came from, so a re-run of an
# analysis reproduces the number that went into it.
_NED_AVAILABLE = True
_NED_TIMEOUT_STRIKES = 0
_NED_STRIKE_LIMIT = 3  # cumulative timeouts across all sources before disabling NED

_NED_OBJLOOKUP_URL = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
_NED_TAP_URL = "https://ned.ipac.caltech.edu/tap/sync"
_SESAME_URL = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV"

_Z_CACHE = None  # lazily loaded dict, see _load_z_cache()


def z_cache_path():
    """
    Path of the on-disk redshift cache (JSON).

    Override with the environment variable ``MORPHEN_Z_CACHE``; defaults to
    ``~/.morphen/redshift_cache.json``.
    """
    import os
    return os.environ.get(
        'MORPHEN_Z_CACHE',
        os.path.join(os.path.expanduser('~'), '.morphen', 'redshift_cache.json'))


def _cache_key(name):
    """Cache keys ignore case and all whitespace, so 'Mrk 331' == 'MRK331'."""
    return ''.join(str(name).split()).upper()


def _load_z_cache():
    global _Z_CACHE
    if _Z_CACHE is None:
        import json
        import os
        path = z_cache_path()
        if os.path.isfile(path):
            try:
                with open(path) as f:
                    _Z_CACHE = json.load(f)
            except Exception:
                _Z_CACHE = {}
        else:
            _Z_CACHE = {}
    return _Z_CACHE


def _save_z_cache():
    import json
    import os
    path = z_cache_path()
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(_load_z_cache(), f, indent=1, sort_keys=True)
    except Exception as e:
        print(f"[cache]     could not write {path} ({type(e).__name__}: {e})")


def set_z(source_name, z, z_err=None, note=None, resolved_name=None):
    """
    Pin a redshift for a source by hand.

    The value is written to the cache with ``service='manual'`` and is used in
    preference to any service afterwards (and is never overwritten by
    ``refresh=True``). Use it for sources no service resolves, or to force the
    exact value quoted in a paper.

    Parameters
    ----------
    source_name : str
    z : float
    z_err : float, optional
    note : str, optional
        Free text kept with the entry, e.g. a bibcode.
    resolved_name : str, optional
    """
    import datetime
    cache = _load_z_cache()
    cache[_cache_key(source_name)] = {
        'z': float(z),
        'z_err': None if z_err is None else float(z_err),
        'service': 'manual',
        'resolved_name': resolved_name or source_name,
        'query_name': source_name,
        'ra': None, 'dec': None, 'separation_arcsec': None,
        'refcode': note,
        'date': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    _save_z_cache()
    print(f"[manual]    z={float(z):.6f} pinned for {source_name}")
    return cache[_cache_key(source_name)]


def forget_z(source_name):
    """Drop a source from the redshift cache (including manual entries)."""
    cache = _load_z_cache()
    if cache.pop(_cache_key(source_name), None) is not None:
        _save_z_cache()
        print(f"[cache]     dropped {source_name}")
        return True
    return False


def _normalize_source_name(name):
    """
    Generate candidate name variants for a source, handling common
    catalog formatting issues that cause NED/SIMBAD lookup failures.

    Returns a list of candidates in priority order (original first).

    Catalogs handled
    ----------------
    IRAS/IRASF  : 'IRASF17132+5313' -> 'IRAS F17132+5313', 'IRAS 17132+5313'
    Zwicky      : 'IIIZw035'        -> 'III Zw 035', 'III Zw 35'
    CGCG        : 'CGCG436-030'     -> 'CGCG 436-030', 'Z 436-30'
                  (SIMBAD files CGCG objects under the prefix 'Z')
    VV          : 'VV250'           -> 'VV 250', 'VV 250a', 'VV 250b'
    MCG         : 'MCG+02-01-051'   -> 'MCG +02-01-051', 'MCG 02-01-051'
    ESO         : 'ESO148-IG002'    -> 'ESO 148-IG 002'
    Generic     : 'NGC5256'         -> 'NGC 5256'  (alpha + digits)
    """
    import re
    candidates = [name]

    def add(c):
        if c not in candidates:
            candidates.append(c)

    # -- IRAS Faint Source Catalog: IRASF[hhmm+ddmm] --------------------------
    m = re.match(r'^IRASF[\s]?(\d{4,5}[+-]\d{4})$', name, re.IGNORECASE)
    if m:
        coords = m.group(1)
        add(f"IRAS F{coords}")
        add(f"IRAS f{coords}")
        add(f"IRAS {coords}")
        return candidates

    # -- Plain IRAS: IRAS[hhmm+ddmm] ------------------------------------------
    m = re.match(r'^IRAS[\s]?(\d{4,5}[+-]\d{4})$', name, re.IGNORECASE)
    if m:
        coords = m.group(1)
        add(f"IRAS {coords}")
        add(f"IRAS F{coords}")
        return candidates

    # -- Zwicky: [I|II|III|IV]Zw[NNN] -----------------------------------------
    m = re.match(r'^(I{1,3}V?|IV)[\s]?(Zw)[\s]?(\d+)$', name, re.IGNORECASE)
    if m:
        roman, zw, num = m.group(1).upper(), 'Zw', m.group(3)
        num_stripped = str(int(num))
        num_padded   = num.zfill(3)
        add(f"{roman} {zw} {num_padded}")
        add(f"{roman} {zw} {num_stripped}")
        return candidates

    # -- CGCG: CGCG[nnn]-[nnn]; SIMBAD knows these as 'Z nnn-nn' ---------------
    m = re.match(r'^(?:CGCG|Z)[\s]?(\d+)-(\d+)$', name, re.IGNORECASE)
    if m:
        field, num = m.group(1), m.group(2)
        add(f"CGCG {field}-{num.zfill(3)}")
        add(f"Z {field}-{str(int(num))}")
        add(f"Z {field}-{num.zfill(3)}")
        return candidates

    # -- VV catalog: VV[NNN] ---------------------------------------------------
    m = re.match(r'^VV[\s]?(\d+)([ab]?)$', name, re.IGNORECASE)
    if m:
        num, suffix = m.group(1), m.group(2).lower()
        base = f"VV {num}"
        add(base)
        if not suffix:
            # Try component suffixes - NED/SIMBAD often require them
            add(f"VV {num}a")
            add(f"VV {num}b")
        else:
            add(f"VV {num}{suffix}")
        return candidates

    # -- MCG: MCG[+/-][ll]-[gg]-[nnn] -----------------------------------------
    m = re.match(r'^MCG[\s]?([+-]?\d{2}-\d{2}-\d{3})$', name, re.IGNORECASE)
    if m:
        coords = m.group(1)
        sign = '+' if not coords.startswith('-') else ''
        add(f"MCG {sign}{coords}")
        add(f"MCG {coords.lstrip('+-')}")
        return candidates

    # -- ESO: ESO[nnn]-[IG|G]?[nnn] -------------------------------------------
    m = re.match(r'^ESO[\s]?(\d+)-?([A-Z]*)[\s]?(\d+)$', name, re.IGNORECASE)
    if m:
        field, kind, num = m.group(1), m.group(2).upper(), m.group(3)
        if kind:
            add(f"ESO {field}-{kind} {num}")
        add(f"ESO {field}-{num}")
        return candidates

    # -- Generic: insert space between alpha prefix and digits -----------------
    # e.g. 'UGC12150' -> 'UGC 12150', 'Mrk331' -> 'Mrk 331'
    m = re.match(r'^([A-Za-z]+)(\d+.*)$', name)
    if m:
        add(f"{m.group(1)} {m.group(2)}")

    return candidates


def _clean(value):
    """Return a plain float, or None for masked/NaN/empty table values."""
    import numpy as np
    if value is None:
        return None
    if hasattr(value, 'mask') and bool(np.all(value.mask)):
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(v) else v


# -- service 1: NED ObjectLookup (REST/JSON) ---------------------------------
def _z_ned_objectlookup(name, timeout=15):
    """
    NED's REST name resolver. Returns a record dict or None.

    This is NED's own name-normalising service (NED-NNS) and it is ~30x faster
    than `Ned.query_object`, which goes through the legacy objsearch CGI. It
    returns the redshift together with its uncertainty, bibcode and quality
    flag, plus NED's preferred name and position.
    """
    import requests
    r = requests.get(_NED_OBJLOOKUP_URL, params={'name': name}, timeout=timeout)
    r.raise_for_status()
    payload = r.json()
    preferred = payload.get('Preferred') or {}
    zblock = preferred.get('Redshift') or {}
    z = _clean(zblock.get('Value'))
    if z is None:
        return None
    position = preferred.get('Position') or {}
    return {
        'z': z,
        'z_err': _clean(zblock.get('Uncertainty')),
        'service': 'NED/ObjectLookup',
        'resolved_name': preferred.get('Name'),
        'ra': _clean(position.get('RA')),
        'dec': _clean(position.get('Dec')),
        'separation_arcsec': 0.0,
        'refcode': zblock.get('RefCode'),
    }


# -- service 2: SIMBAD TAP, exact identifier ---------------------------------
def _z_simbad_ident(name, timeout=15):
    """SIMBAD TAP lookup on the `ident` table (exact identifier match)."""
    from astroquery.simbad import Simbad
    Simbad.TIMEOUT = timeout
    escaped = name.replace("'", "''")
    table = Simbad.query_tap(
        f"SELECT TOP 1 main_id, ra, dec, rvz_redshift "
        f"FROM basic JOIN ident ON ident.oidref = basic.oid "
        f"WHERE ident.id = '{escaped}'")
    if table is None or len(table) == 0:
        return None
    z = _clean(table['rvz_redshift'][0])
    if z is None:
        return None
    return {
        'z': z,
        'z_err': None,
        'service': 'SIMBAD/ident',
        'resolved_name': str(table['main_id'][0]),
        'ra': _clean(table['ra'][0]),
        'dec': _clean(table['dec'][0]),
        'separation_arcsec': 0.0,
        'refcode': None,
    }


# -- service 3a: Sesame, name -> position ------------------------------------
def _coords_sesame(name, timeout=15):
    """
    CDS Sesame name resolver: name -> (ra, dec) in degrees, or None.

    Sesame queries the SIMBAD, NED and VizieR name dictionaries in turn, so it
    resolves names that no single catalogue service matches by string.
    """
    import re
    import requests
    r = requests.get(f"{_SESAME_URL}?{requests.utils.quote(name)}", timeout=timeout)
    r.raise_for_status()
    m = re.search(r'^%J\s+([-\d.]+)\s+([-+\d.]+)', r.text, re.MULTILINE)
    if m is None:
        return None
    return float(m.group(1)), float(m.group(2))


# -- service 3b: cone searches ------------------------------------------------
def _z_ned_cone(ra, dec, radius_arcsec=30.0, timeout=30):
    """Nearest NED object with a redshift within `radius_arcsec` of (ra, dec)."""
    import io
    import requests
    from astropy.io.votable import parse_single_table

    radius_deg = radius_arcsec / 3600.0
    # `z IS NOT NULL` in the query, not in a loop over the nearest rows: in a
    # crowded field (M82, say) the 20 nearest NED entries can all be knots and
    # X-ray sources with no redshift at all.
    query = (f"SELECT TOP 5 prefname, z, "
             f"DISTANCE(POINT('J2000', ra, dec), POINT('J2000', {ra}, {dec})) AS d "
             f"FROM NEDTAP.objdir "
             f"WHERE CONTAINS(POINT('J2000', ra, dec), "
             f"CIRCLE('J2000', {ra}, {dec}, {radius_deg})) = 1 "
             f"AND z IS NOT NULL "
             f"ORDER BY d ASC")
    r = requests.post(_NED_TAP_URL,
                      data={'REQUEST': 'doQuery', 'LANG': 'ADQL',
                            'FORMAT': 'votable', 'QUERY': query},
                      timeout=timeout)
    r.raise_for_status()
    table = parse_single_table(io.BytesIO(r.content)).to_table()
    for row in table:
        z = _clean(row['z'])
        if z is None:
            continue
        return {
            'z': z,
            'z_err': None,
            'service': 'NED/TAP-cone',
            'resolved_name': str(row['prefname']).strip() or None,
            'ra': ra, 'dec': dec,
            'separation_arcsec': float(row['d']) * 3600.0,
            'refcode': None,
        }
    return None


def _z_simbad_cone(ra, dec, radius_arcsec=30.0, timeout=15):
    """Nearest SIMBAD object with a redshift within `radius_arcsec` of (ra, dec)."""
    from astroquery.simbad import Simbad
    Simbad.TIMEOUT = timeout
    radius_deg = radius_arcsec / 3600.0
    table = Simbad.query_tap(
        f"SELECT TOP 5 main_id, otype, rvz_redshift, "
        f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) AS d "
        f"FROM basic "
        f"WHERE CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius_deg})) = 1 "
        f"AND rvz_redshift IS NOT NULL "
        f"ORDER BY d ASC")
    if table is None or len(table) == 0:
        return None
    for row in table:
        z = _clean(row['rvz_redshift'])
        if z is None:
            continue
        return {
            'z': z,
            'z_err': None,
            'service': 'SIMBAD/cone',
            'resolved_name': str(row['main_id']),
            'ra': ra, 'dec': dec,
            'separation_arcsec': float(row['d']) * 3600.0,
            'refcode': None,
        }
    return None


# -- service 4: legacy astroquery NED (slow, circuit-broken) ------------------
def _z_ned_astroquery(name, timeout=20, verbose=True):
    """
    `Ned.query_object` - the original path, kept as a last resort.

    It scrapes NED's legacy objsearch CGI and routinely takes 15-60 s or times
    out, so it is circuit-broken: after `_NED_STRIKE_LIMIT` cumulative timeouts
    it is skipped for the rest of the session.
    """
    global _NED_AVAILABLE, _NED_TIMEOUT_STRIKES
    from requests.exceptions import ReadTimeout, ConnectionError as ReqConnectionError
    from astroquery.ipac.ned import Ned

    if not _NED_AVAILABLE:
        return None
    Ned.TIMEOUT = timeout
    try:
        table = Ned.query_object(name)
    except ReqConnectionError:
        _NED_AVAILABLE = False
        if verbose:
            print("[NED]       connection error - disabling for this session.")
        return None
    except (ReadTimeout, TimeoutError):
        _NED_TIMEOUT_STRIKES += 1
        if verbose:
            print(f"[NED]       timeout on '{name}' "
                  f"(strike {_NED_TIMEOUT_STRIKES}/{_NED_STRIKE_LIMIT}).")
        if _NED_TIMEOUT_STRIKES >= _NED_STRIKE_LIMIT:
            _NED_AVAILABLE = False
            if verbose:
                print("[NED]       too many timeouts - disabling for this session.")
        return None
    data = table['Redshift'].data.data
    if data.shape[0] == 0:
        return None
    z = _clean(data[0])
    if z is None:
        return None
    _NED_TIMEOUT_STRIKES = 0
    return {
        'z': z,
        'z_err': None,
        'service': 'NED/astroquery',
        'resolved_name': str(table['Object Name'][0]),
        'ra': _clean(table['RA'][0]),
        'dec': _clean(table['DEC'][0]),
        'separation_arcsec': 0.0,
        'refcode': None,
    }


def _image_center_coords(imagename):
    """Sky coordinates of the centre of a FITS image, in degrees."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(imagename) as hdul:
        header = hdul[0].header
        shape = hdul[0].data.shape
    wcs = WCS(header).celestial
    ny, nx = shape[-2], shape[-1]
    ra, dec = wcs.all_pix2world([[(nx - 1) / 2.0, (ny - 1) / 2.0]], 0)[0]
    return float(ra) % 360.0, float(dec)


def resolve_source_coordinates(source_name, timeout=15, verbose=True):
    """
    Resolve a source name to sky coordinates, robustly.

    Same service chain as `find_z` (NED ObjectLookup -> SIMBAD TAP -> Sesame),
    so it resolves the same awkward names and answers in ~0.5 s instead of the
    15-60 s of `Ned.query_object`. Positions already stored in the redshift
    cache are reused.

    Parameters
    ----------
    source_name : str
    timeout : int, optional
        Per-service timeout in seconds. Default: 15.
    verbose : bool, optional

    Returns
    -------
    astropy.coordinates.SkyCoord or None
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    cached = _load_z_cache().get(_cache_key(source_name))
    if cached is not None and cached.get('ra') is not None:
        return SkyCoord(cached['ra'], cached['dec'], unit=(u.deg, u.deg))

    for name in _normalize_source_name(source_name):
        for query in (_z_ned_objectlookup, _z_simbad_ident):
            try:
                record = query(name, timeout=timeout)
            except Exception:
                continue
            if record is not None and record.get('ra') is not None:
                if verbose:
                    print(f"[{record['service']}] {source_name} -> "
                          f"{record['resolved_name']} "
                          f"RA={record['ra']:.6f} Dec={record['dec']:+.6f}")
                return SkyCoord(record['ra'], record['dec'], unit=(u.deg, u.deg))
        try:
            position = _coords_sesame(name, timeout=timeout)
        except Exception:
            position = None
        if position is not None:
            if verbose:
                print(f"[Sesame]    {source_name} -> "
                      f"RA={position[0]:.6f} Dec={position[1]:+.6f}")
            return SkyCoord(position[0], position[1], unit=(u.deg, u.deg))

    if verbose:
        print(f"[!] Could not resolve coordinates for {source_name}")
    return None


def find_z(source_name=None, coordinates=None, imagename=None,
           prefer='ned', radius_arcsec=30.0, use_cache=True, refresh=False,
           full=False, return_luminosity_distance=False,
           timeout=15, verbose=True):
    """
    Find the redshift of a source, robustly.

    Tries several independent services in a fixed order (see the module note
    above), caching every success so a repeated call - or a re-run of a whole
    analysis - never depends on a service being up, and always returns the same
    number that went into the previous run.

    A source can be identified by name, by coordinates, or by a FITS image
    whose centre lands on it; the coordinate paths bypass name resolution
    entirely and are the reliable option for objects with awkward names.

    Parameters
    ----------
    source_name : str, optional
        Source name, e.g. 'CGCG436-030', 'VV705', 'IRASF17132+5313'.
    coordinates : tuple or SkyCoord, optional
        (ra, dec) in degrees, or a SkyCoord. Skips name resolution and goes
        straight to the cone searches; `source_name`, if also given, is then
        used only as the label and cache key.
    imagename : str, optional
        FITS image; its centre is used as the search position. Used only if
        `coordinates` is not given and name lookup fails.
    prefer : {'ned', 'simbad'}, optional
        Which service is consulted first. NED and SIMBAD quote redshifts from
        different references, so this fixes which value you get (differences are
        typically ~100 km/s, i.e. ~1% in luminosity distance at z~0.03).
        Default: 'ned'.
    radius_arcsec : float, optional
        Cone-search radius for the positional fallbacks. Default: 30.
    use_cache : bool, optional
        Read from / write to the on-disk cache (see `z_cache_path`). Default: True.
    refresh : bool, optional
        Ignore a cached value and query again (manual entries set with `set_z`
        are still honoured). Default: False.
    full : bool, optional
        Return the full record dict (z, z_err, service, resolved_name, ra, dec,
        separation_arcsec, refcode, date) instead of a bare float. Default: False.
    return_luminosity_distance : bool, optional
        Also return the Planck18 luminosity distance. Default: False.
    timeout : int, optional
        Per-service timeout in seconds. Default: 15.
    verbose : bool, optional
        Print which service answered. Default: True.

    Returns
    -------
    z : float, dict or None
        The redshift (or the full record if `full=True`), None if nothing
        resolved.
    luminosity_distance : astropy.units.Quantity, optional
        Only if `return_luminosity_distance=True`.

    Examples
    --------
    >>> find_z('CGCG436-030')
    >>> find_z('CGCG436-030', full=True)['service']
    >>> find_z(coordinates=(20.010964, 14.361750))
    >>> find_z('IRASF17132+5313', imagename='IRASF17132+5313_C_band.fits')
    >>> set_z('MyOddSource', 0.0421, note='2019ApJ...123..456X')  # pin by hand
    """
    import datetime

    if source_name is None and coordinates is None and imagename is None:
        raise ValueError("give at least one of source_name, coordinates, imagename")

    label = source_name if source_name is not None else (imagename or 'position')
    record = None

    # -- 0. cache / manual override -------------------------------------------
    if use_cache and source_name is not None:
        cached = _load_z_cache().get(_cache_key(source_name))
        if cached is not None and (not refresh or cached.get('service') == 'manual'):
            if verbose:
                print(f"[cache]     z={cached['z']:.6f} for {source_name} "
                      f"(via {cached['service']}, {cached.get('date', '?')})")
            record = cached

    # -- name-based services ---------------------------------------------------
    # Explicit coordinates override name resolution: passing them is how you say
    # "this position, whatever the name resolves to".
    if record is None and source_name is not None and coordinates is None:
        name_services = [('NED/ObjectLookup', _z_ned_objectlookup),
                         ('SIMBAD/ident', _z_simbad_ident)]
        if str(prefer).lower().startswith('simbad'):
            name_services.reverse()

        candidates = _normalize_source_name(source_name)
        for service_name, query in name_services:
            for name in candidates:
                try:
                    record = query(name, timeout=timeout)
                except Exception as e:
                    if verbose:
                        print(f"[{service_name}] failed on '{name}' "
                              f"({type(e).__name__}: {e})")
                    break  # service is down; move to the next one
                if record is not None:
                    record['query_name'] = name
                    break
            if record is not None:
                break

    # -- positional fallback ---------------------------------------------------
    if record is None:
        position = None
        if coordinates is not None:
            if hasattr(coordinates, 'ra'):  # SkyCoord
                position = (float(coordinates.ra.deg), float(coordinates.dec.deg))
            else:
                position = (float(coordinates[0]), float(coordinates[1]))
        elif source_name is not None:
            for name in _normalize_source_name(source_name):
                try:
                    position = _coords_sesame(name, timeout=timeout)
                except Exception as e:
                    if verbose:
                        print(f"[Sesame]    failed ({type(e).__name__}: {e})")
                    break
                if position is not None:
                    if verbose:
                        print(f"[Sesame]    {source_name} -> "
                              f"RA={position[0]:.6f} Dec={position[1]:+.6f} "
                              + (f"(as '{name}')" if name != source_name else ""))
                    break
        if position is None and imagename is not None:
            try:
                position = _image_center_coords(imagename)
                if verbose:
                    print(f"[image]     centre of {imagename} -> "
                          f"RA={position[0]:.6f} Dec={position[1]:+.6f}")
            except Exception as e:
                if verbose:
                    print(f"[image]     could not read WCS from {imagename} "
                          f"({type(e).__name__}: {e})")

        if position is not None:
            cone_services = [_z_ned_cone, _z_simbad_cone]
            if str(prefer).lower().startswith('simbad'):
                cone_services.reverse()
            for cone in cone_services:
                try:
                    record = cone(position[0], position[1],
                                  radius_arcsec=radius_arcsec,
                                  timeout=max(timeout, 30))
                except Exception as e:
                    if verbose:
                        print(f"[cone]      {cone.__name__} failed "
                              f"({type(e).__name__}: {e})")
                    continue
                if record is not None:
                    record['query_name'] = source_name
                    break

    # -- last resort: the legacy astroquery NED path ---------------------------
    # One attempt on the name as given: NED's CGI runs its own name interpreter,
    # and each attempt costs 20-60 s when the service is degraded.
    if record is None and source_name is not None:
        try:
            record = _z_ned_astroquery(source_name, timeout=max(timeout, 20),
                                       verbose=verbose)
        except Exception as e:
            if verbose:
                print(f"[NED]       query_object failed "
                      f"({type(e).__name__}: {e})")
            record = None
        if record is not None:
            record['query_name'] = source_name

    if record is None:
        if verbose:
            print(f"[!] No redshift found for {label}. Tried "
                  f"{_normalize_source_name(source_name) if source_name else 'position'}. "
                  f"Pin it by hand with set_z('{label}', z).")
        return None

    # -- report / cache --------------------------------------------------------
    if 'date' not in record:
        record['date'] = datetime.datetime.now().isoformat(timespec='seconds')
        if verbose:
            resolved = record.get('resolved_name')
            note = f" (as '{resolved}')" if resolved and resolved != source_name else ""
            separation = record.get('separation_arcsec') or 0.0
            if separation > 1.0:
                note += f" [{separation:.1f}\" away]"
            print(f"[{record['service']}] z={record['z']:.6f} for {label}{note}")
        if use_cache and source_name is not None:
            _load_z_cache()[_cache_key(source_name)] = record
            _save_z_cache()

    separation = record.get('separation_arcsec') or 0.0
    if verbose and separation > radius_arcsec / 2:
        print(f"[!]         nearest match is {separation:.1f}\" from the search "
              f"position - check it is the right object.")

    result = record if full else record['z']
    if return_luminosity_distance:
        from astropy.cosmology import Planck18
        lum_dist = Planck18.luminosity_distance(record['z'])
        print(f"Luminosity distance for {label} "
              f"(z={record['z']:.4f}): {lum_dist:.2f}")
        return result, lum_dist
    return result


def find_z_NED(source_name, return_luminosity_distance=False,
               ned_timeout=20, verbose=True, **kwargs):
    """
    Backwards-compatible alias for `find_z` (see it for the full interface).

    Kept so existing notebooks keep working; new code should call `find_z`,
    which also accepts coordinates, a FITS image, and `full=True`.
    """
    return find_z(source_name,
                  return_luminosity_distance=return_luminosity_distance,
                  timeout=ned_timeout, verbose=verbose, **kwargs)



# def find_z_NED(source_name, return_luminosity_distance=False,
#                max_retries=3, timeout=60, backoff=5):
#     """
#     Find the redshift of a source (by name) from NED.

#     Parameters
#     ----------
#     source_name : str
#         Source name. Example: 'VV705'
#     return_luminosity_distance : bool, optional
#         If True, also calculate and return the luminosity distance.
#         Default: False
#     max_retries : int, optional
#         Number of attempts before giving up. Default: 3
#     timeout : int, optional
#         Per-attempt read timeout in seconds. Default: 60
#     backoff : int, optional
#         Seconds to wait between retries (doubles each attempt). Default: 5

#     Returns
#     -------
#     redshift_NED : float or None
#         The redshift value from NED, or None if not found / all retries failed.
#     luminosity_distance : astropy.units.Quantity, optional
#         The luminosity distance (only returned if return_luminosity_distance=True)
#     """
#     import time
#     from astroquery.ipac.ned import Ned
#     from astroquery.exceptions import RemoteServiceError
#     from requests.exceptions import ReadTimeout, ConnectionError as ReqConnectionError

#     last_exc = None
#     wait = backoff

#     for attempt in range(1, max_retries + 1):
#         try:
#             Ned.TIMEOUT = timeout
#             result_table = Ned.query_object(source_name)
#             redshift_NED = result_table['Redshift'].data.data

#             if redshift_NED.shape[0] == 0:
#                 return None

#             z = redshift_NED[0]

#             if return_luminosity_distance:
#                 from astropy.cosmology import Planck18
#                 lum_dist = Planck18.luminosity_distance(z)
#                 print(f"Luminosity distance for {source_name} (z={z:.4f}): {lum_dist:.2f}")
#                 return z, lum_dist
#             else:
#                 return z

#         except (ReadTimeout, ReqConnectionError, TimeoutError) as e:
#             last_exc = e
#             print(f"  [{source_name}] NED timeout on attempt {attempt}/{max_retries}, "
#                   f"retrying in {wait}s...")
#             time.sleep(wait)
#             wait *= 2  # exponential backoff

#         except RemoteServiceError as e:
#             # NED returned an error (e.g. source not found in NED)
#             print(f"  [{source_name}] NED service error: {e}")
#             return None

#         except Exception as e:
#             print(f"  [{source_name}] Unexpected error querying NED: {e}")
#             return None

#     print(f"  [{source_name}] All {max_retries} NED attempts failed: {last_exc}")
#     return None
