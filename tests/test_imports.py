"""Every module must at least import.

The suite drives the viewer through offscreen widgets but never imports most
of `gui`, so a bad import there used to survive a green run: a linter's
autofix once removed a re-exported constant and nothing noticed. Importing is
safe — it is constructing a QApplication that the convention forbids.
"""

import importlib
import os
import pkgutil
import subprocess
import sys
from pathlib import Path

import pytest

import microscope

MODULES = sorted(
    info.name
    for info in pkgutil.walk_packages(microscope.__path__, "microscope.")
    if not info.name.endswith("__main__")      # that one runs the program
)


def test_the_package_has_modules_to_check():
    assert len(MODULES) > 25


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    importlib.import_module(name)


def test_the_public_api_is_what_it_says():
    for name in microscope.__all__:
        assert hasattr(microscope, name), name


def test_the_viewer_starts_without_matplotlib():
    """matplotlib was half the viewer's start-up time, paid for a spectrum most
    files do not have; it loads with the first spectrum tab now. Checked in a
    fresh interpreter, since this one has long since imported everything."""
    script = ("import sys, microscope.gui.app; "
              "print(any(m.split('.')[0] == 'matplotlib' for m in sys.modules))")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
    proc = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                          env=env, cwd=str(Path(__file__).parent.parent))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False"
