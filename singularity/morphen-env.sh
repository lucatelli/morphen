# Sourced from the image's %environment (POSIX sh).
#
# When $HOME is missing or read-only (--containall with a read-only -H,
# --no-home, some HPC setups), point every tool that writes under $HOME at
# ${TMPDIR:-/tmp}/<user>-morphen instead. Without this, jupyter lab fails to
# start (Permission denied: ~/.local), astropy cannot cache IERS downloads and
# matplotlib rebuilds its font cache on every start. With a writable $HOME
# this does nothing. Explicitly set variables are never overridden.
#
# Note: under --containall the container's /tmp is a 16 MB tmpfs (Singularity
# "sessiondir max size"); on HPC set SINGULARITYENV_TMPDIR to a scratch bind.
if [ -z "${HOME:-}" ] || ! [ -d "$HOME" ] || ! [ -w "$HOME" ]; then
    _mu="${USER:-$(id -un 2>/dev/null || id -u)}"
    _mb="${TMPDIR:-/tmp}/${_mu}-morphen"
    : "${XDG_CACHE_HOME:=$_mb/cache}"
    : "${XDG_CONFIG_HOME:=$_mb/config}"
    : "${XDG_DATA_HOME:=$_mb/share}"
    : "${XDG_STATE_HOME:=$_mb/state}"
    : "${MPLCONFIGDIR:=$_mb/config/matplotlib}"
    : "${IPYTHONDIR:=$_mb/ipython}"
    : "${JUPYTER_CONFIG_DIR:=$_mb/jupyter}"
    : "${JUPYTER_DATA_DIR:=$_mb/share/jupyter}"
    : "${JUPYTER_RUNTIME_DIR:=$_mb/share/jupyter/runtime}"
    : "${JUPYTERLAB_SETTINGS_DIR:=$_mb/jupyter/lab/user-settings}"
    : "${JUPYTERLAB_WORKSPACES_DIR:=$_mb/jupyter/lab/workspaces}"
    export XDG_CACHE_HOME XDG_CONFIG_HOME XDG_DATA_HOME XDG_STATE_HOME MPLCONFIGDIR \
           IPYTHONDIR JUPYTER_CONFIG_DIR JUPYTER_DATA_DIR JUPYTER_RUNTIME_DIR \
           JUPYTERLAB_SETTINGS_DIR JUPYTERLAB_WORKSPACES_DIR
    mkdir -p "$XDG_CACHE_HOME/astropy" "$XDG_CONFIG_HOME/astropy" "$XDG_DATA_HOME" \
             "$XDG_STATE_HOME" "$MPLCONFIGDIR" "$IPYTHONDIR" "$JUPYTER_CONFIG_DIR" \
             "$JUPYTER_RUNTIME_DIR" 2>/dev/null
    chmod 700 "$_mb" "$JUPYTER_RUNTIME_DIR" 2>/dev/null
    unset _mu _mb
fi
