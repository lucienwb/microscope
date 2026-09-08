"""Every module must at least import.

The suite drives the viewer through offscreen widgets but never imports most
of `gui`, so a bad import there used to survive a green run: a linter's
autofix once removed a re-exported constant and nothing noticed. Importing is
safe — it is constructing a QApplication that the convention forbids.
"""

import importlib
import pkgutil

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
