#!/usr/bin/env bash
# CASA steps shared by the conda and uv recipes (run from %post as root):
#
#   bash /opt/morphen/bin/casa-post.sh /path/to/env/bin/python
#
#  1. bake casarundata + measures into /opt/casa/data (read-only at runtime;
#     /opt/casa/casasiteconfig.py points casaconfig at it, see that file)
#  2. extract the casaplotms/casaviewer AppImages (no FUSE in Singularity)
#     and replace each .AppImage with a shim that runs the extracted AppRun
#  3. check that casatools imports against the baked data
set -euo pipefail
PY=${1:?usage: casa-post.sh /path/to/python}
export MORPHEN_PY=$PY

# 1. CASA data
mkdir -p /opt/casa/data
"$PY" -m casaconfig --measurespath /opt/casa/data --update-all
"$PY" -m casaconfig --measurespath /opt/casa/data --current-data
chmod -R a+rX /opt/casa

# 2. AppImages
# --- casaplotms / casaviewer without FUSE ---------------------------------
# The pip wheels ship type-2 AppImages that need FUSE (not available in
# Singularity), and their old runtime has no APPIMAGE_EXTRACT_AND_RUN.
# Extract each one once at build time and replace the .AppImage file (its
# path is hard-coded in casaplotms/private/plotmstool.py and
# casaviewer/private/config.py) with a small sh shim that execs the
# extracted AppRun. Path-agnostic: the env's python locates the packages
# (set MORPHEN_PY to that python if it is not first on PATH). Idempotent.
PY="${MORPHEN_PY:-python}"
for _pkg in casaplotms casaviewer; do
    _dir=$("$PY" -c "import importlib.util as u; s=u.find_spec('$_pkg'); print(s.submodule_search_locations[0] if s else '')")
    [ -n "$_dir" ] || { echo "skip $_pkg: not installed"; continue; }
    _img="$_dir/__bin__/$_pkg-x86_64.AppImage"
    _app="$_dir/__bin__/$_pkg.AppDir"
    if [ -f "$_img" ] && ! grep -q 'morphen-appimage-shim' "$_img" 2>/dev/null; then
        rm -rf "$_app" "$_dir/__bin__/squashfs-root"
        ( cd "$_dir/__bin__" && "./$_pkg-x86_64.AppImage" --appimage-extract >/dev/null 2>&1 ) \
            || { echo "ERROR: extracting $_img failed" >&2; exit 1; }
        mv "$_dir/__bin__/squashfs-root" "$_app"
        chmod -R a+rX "$_app"          # extracted dirs are created 0700
        rm -f "$_img"
        cat > "$_img" <<'SHIM'
#!/bin/sh
# morphen-appimage-shim: the AppImage was extracted at image build time
# (no FUSE in Singularity); run the extracted AppRun instead.
# The app links its own casacore, which ignores casaconfig: point its
# measures.directory at the same data the site config uses (CASARCFILES).
_d="${TMPDIR:-/tmp}/${USER:-$(id -un 2>/dev/null || id -u)}-morphen"
if [ -z "${CASARCFILES:-}" ]; then
    _m=/opt/casa/data
    if [ -n "${MORPHEN_CASADATA:-}" ] && [ -f "$MORPHEN_CASADATA/geodetic/readme.txt" ]; then
        _m=$MORPHEN_CASADATA
    fi
    if [ -d "$_m/geodetic" ] && mkdir -p "$_d" 2>/dev/null &&
       printf 'measures.directory: %s\n' "$_m" > "$_d/casarc.$$" &&
       mv -f "$_d/casarc.$$" "$_d/casarc"; then
        export CASARCFILES="$_d/casarc"
    fi
fi
# Its casacore also reads/creates ~/.casa/rc and fails to load measures
# (e.g. plotms "Error during cache loading") when $HOME is read-only.
if ! [ -w "${HOME:-/nonexistent}" ]; then
    mkdir -p "$_d/apphome" 2>/dev/null && export HOME="$_d/apphome"
fi
# No usable X display (none set, or a local :N socket not visible, as under
# --containall): render offscreen so plotms/imview exports still work.
if [ -z "${QT_QPA_PLATFORM:-}" ]; then
    case "${DISPLAY:-}" in
        '') export QT_QPA_PLATFORM=offscreen ;;
        :*) _n=${DISPLAY#:}; _n=${_n%%.*}
            [ -S "/tmp/.X11-unix/X$_n" ] || export QT_QPA_PLATFORM=offscreen ;;
    esac
fi
exec "$(dirname "$0")/@PKG@.AppDir/AppRun" "$@"
SHIM
        sed -i "s/@PKG@/$_pkg/" "$_img"
        chmod 755 "$_img"
    fi
    [ -x "$_app/AppRun" ] || { echo "ERROR: $_app/AppRun missing" >&2; exit 1; }
    echo "$_pkg: $(du -sh "$_app" | cut -f1) at $_app"
done
unset _pkg _dir _img _app

# 3. sanity: site config picked up, baked data used, nothing to update
cd "${TMPDIR:-/tmp}"
"$PY" - <<'PYEOF'
import casaconfig.config as c
assert '/opt/casa/casasiteconfig.py' in c.load_success(), c.load_success()
assert c.measurespath == '/opt/casa/data' and not c.measures_auto_update, c.measurespath
import casatools
me = casatools.measures()
me.doframe(me.observatory('VLA'))
print('casa-post: casatools %s OK, measurespath %s' % (casatools.version_string(), c.measurespath))
PYEOF
