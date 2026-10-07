#!/usr/bin/env python
"""Fix petrofit 0.5.x's `_astropy_init.py` for astropy >= 5.

Run once, inside the activated environment, after creating it from
`environment_py310.yml` (or any env on python < 3.12):

    conda env create -f environment_py310.yml
    conda activate morphen310_v1_0
    python patch_petrofit.py

Why: on python < 3.12 pip can only install petrofit 0.5.0, whose
`petrofit/_astropy_init.py` imports `update_default_config` from
`astropy.config.configuration`. astropy 5 removed it, so `import petrofit`
raises ImportError. petrofit 0.6.0 (python >= 3.12, see `environment.yml`)
no longer has this file, and nothing is patched there.

The file is located with `importlib.util.find_spec`, which does not execute
petrofit's `__init__` (it would fail), so the miniconda path and env name are
never hard-coded. The original is kept as `_astropy_init.py.bak`. Running the
script again is a no-op.
"""
import importlib.metadata
import importlib.util
import os
import shutil
import sys

PATCHED = """\
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['__version__']

# this indicates whether or not we are in the package's setup.py
try:
    _ASTROPY_SETUP_
except NameError:
    import builtins
    builtins._ASTROPY_SETUP_ = False

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''


if not _ASTROPY_SETUP_:  # noqa
    import os
    from warnings import warn

    # Create the test function for self test
    from astropy.tests.runner import TestRunner
    test = TestRunner.make_test_runner_in(os.path.dirname(__file__))
    test.__test__ = False
    __all__ += ['test']

    # Configuration update code removed for astropy 5.0+ compatibility
    # The update_default_config function was removed in astropy 5.0
    # This configuration management is not critical for petrofit functionality
    # (patched by morphen/patch_petrofit.py)
"""


def main():
    spec = importlib.util.find_spec("petrofit")
    if spec is None or not spec.submodule_search_locations:
        sys.exit(f"[patch_petrofit] petrofit is not installed in {sys.prefix}")

    version = importlib.metadata.version("petrofit")
    target = os.path.join(spec.submodule_search_locations[0], "_astropy_init.py")
    print(f"[patch_petrofit] petrofit {version} in {sys.prefix}")

    if not os.path.isfile(target):
        print("[patch_petrofit] no _astropy_init.py (petrofit >= 0.6): nothing to do")
        return

    with open(target) as f:
        source = f.read()
    if "from astropy.config.configuration import" not in source:
        print(f"[patch_petrofit] already patched: {target}")
        return
    if not version.startswith("0.5"):
        sys.exit(f"[patch_petrofit] unexpected petrofit {version} with the old "
                 f"astropy config import; refusing to patch {target}")

    backup = target + ".bak"
    if not os.path.exists(backup):
        shutil.copy2(target, backup)
    with open(target, "w") as f:
        f.write(PATCHED)
    print(f"[patch_petrofit] patched {target} (original kept as {backup})")


if __name__ == "__main__":
    main()
