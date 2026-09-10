# Development

```bash
git clone https://github.com/lucienwb/microscope
cd microscope
pip install -e ".[dev]"
```

## Checks

All three must pass, and CI runs them on every push:

```bash
pytest                                  # the test suite, headless
ruff check microscope tests scripts     # lint
mypy                                    # types (the Qt widgets are excluded)
```

The tests need no display and no OpenGL — no test may construct a Qt
application — so they run the same on a laptop and a CI runner. Parser tests
use real program outputs from the cclib test suite, kept in `tests/data/`.

Ruff's name checks earn their keep in a codebase that moves code between
modules: a name used but never imported is exactly what slips through a
refactor, and it has caught several.

## Checking pictures

A test can say a function returned the right numbers; it cannot say a figure
looks right. Rendering work is checked by rendering to a file and looking at it,
and the command line makes that one line:

```bash
scope -s tests/data/trp.log -o /tmp/look.png
```

## The animations

Every animation in this documentation is recorded by one script, rendering
offscreen through the same code the viewer uses, so a clip cannot drift from
what the program draws:

```bash
python scripts/record_demo.py              # all of them
python scripts/record_demo.py style orbit  # just these
```

The 3-D clips need a working OpenGL driver; the Lewis ones do not.

## This site

```bash
pip install -e ".[docs]"
mkdocs serve          # live preview at http://127.0.0.1:8000
mkdocs build --strict # what CI runs: any broken link fails the build
```

The command-line and style-file references are generated from the code at build
time, by `scripts/mkdocs_hooks.py`, so they cannot fall out of date. Pushing to
`main` rebuilds and publishes the site through `.github/workflows/docs.yml`.

MkDocs is pinned below 2.0: that release removes the plugin system the theme and
the API reference are built on, with no migration path.

## Ground rules

- **Pure visualization.** It draws what programs produced and does not compute
  chemistry. Evaluating orbitals from a basis set was considered and rejected.
- **Few dependencies.** `numpy`, `scipy`, `matplotlib`, `pandas`, `PySide6`,
  `PyOpenGL`, and nothing else at run time. cclib is an optional fallback, never
  imported unless needed.
- **Parsers in-house.** Every format is read by code in this repository.
