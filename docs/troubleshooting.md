# Troubleshooting

Most problems come down to one of two things: the machine cannot give
microscope the OpenGL it draws with, or a file says something the parsers do
not expect. Each has its own section below.

## Does this machine have the OpenGL microscope needs?

Every 3-D picture, in the window or from `scope -s`, is drawn with **OpenGL 3.3
or later** (core profile). A desktop or laptop from the last ten years has it.
A virtual machine, a remote session or a cluster node might not. The quickest
test is to make a figure:

```bash
scope -s tests/data/trp.log -o check.png
```

If that writes a picture of tryptophan, the graphics are fine. For more
detail, this script asks for the same kind of context microscope does, and
says which driver answered:

```python
from OpenGL import GL
from PySide6.QtGui import (QGuiApplication, QOffscreenSurface, QOpenGLContext,
                           QSurfaceFormat)

app = QGuiApplication([])
fmt = QSurfaceFormat()
fmt.setVersion(3, 3)
fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
context = QOpenGLContext()
context.setFormat(fmt)
surface = QOffscreenSurface()
surface.setFormat(fmt)
surface.create()
if not (context.create() and context.makeCurrent(surface)):
    raise SystemExit("no OpenGL 3.3 context: microscope cannot draw here")
print("OpenGL", GL.glGetString(GL.GL_VERSION).decode())
print("renderer", GL.glGetString(GL.GL_RENDERER).decode())
```

The version it prints should be 3.3 or higher. A renderer named `llvmpipe` is
Mesa drawing in software: slower than a graphics card, but correct.

## Linux

### The window will not open

Qt prints which platform plugin failed. The most common one is on Qt 6.5
and later, which needs `libxcb-cursor0` for its X11 plugin and says so:

```bash
sudo apt install libxcb-cursor0          # Debian, Ubuntu
sudo dnf install xcb-util-cursor         # Fedora
```

A minimal system, such as a container or a fresh server image, may be missing
more of Qt's runtime libraries. On Debian and Ubuntu these are the ones the
project's own CI installs:

```bash
sudo apt install libegl1 libgl1 libx11-6 libxkbcommon0 libfontconfig1 \
                 libfreetype6 libdbus-1-3
```

If it still fails, `QT_DEBUG_PLUGINS=1 scope` makes Qt name the exact library
it could not load.

### A cluster node, or any machine with no display

`scope -s` needs no window. The reliable way to run it where there is no
display at all is under a virtual one:

```bash
xvfb-run -a scope -s opt.log -o figure.png
```

`xvfb-run` comes with the `xvfb` package. Under it, Mesa's software renderer
does the drawing. Writing an orbital to a cube file (`--mo homo -o homo.cube`),
a Lewis structure to `.cdxml` or `.mol`, or a spectrum plot draws nothing in
3-D, so those run on a node as they are.

### A black window, or an error about the OpenGL context

This usually means the driver in use offers only an old OpenGL. That is common
in virtual machines and under WSL. Mesa's software renderer does support
OpenGL 3.3, and this makes Mesa use it:

```bash
LIBGL_ALWAYS_SOFTWARE=1 scope mycalc.log
```

## Working on a remote machine

**Over `ssh -X`**, the window is drawn by your own computer through the
forwarded connection. That connection usually offers only a much older OpenGL
than 3.3, so the viewer tends not to start. Two ways round it:

- run `scope -s` on the remote machine and copy the image back, or
- use a remote desktop (VNC, X2Go, NoMachine) that draws on the remote machine
  and sends you the picture.

**Over Windows Remote Desktop**, a session may be given a basic OpenGL 1.1
driver instead of the graphics card, and the window is then blank. If it
works when you sit at the machine but not over Remote Desktop, that is the
reason. Your administrator can let sessions use the graphics card; a figure
from `scope -s`, run at the machine or from a scheduled task, avoids the
question.

## Files

### "could not recognize this output file"

For `.log` and `.out` files the program is recognized from the contents, not
the name. The [file formats](reference/formats.md) page lists what is read
natively. If [cclib](https://cclib.github.io) is installed, microscope tries
it for anything the native parsers do not handle:

```bash
pip install cclib
```

If a Gaussian, ORCA or Q-Chem file is not recognized, that is a bug. Please
[open an issue](https://github.com/lucienwb/microscope/issues) and attach the
file, or its first hundred lines if it is large.

### "this file's basis-set conventions were not recognized"

This is almost always a Molden file. Programs follow the Molden format
loosely: each orders, normalizes and signs its basis functions a little
differently. microscope tries each known convention, and keeps the one under
which the orbitals come out orthonormal, as a correct reading must. This
warning means none of them did, so the orbitals may be drawn wrong. The
structure and the rest of the file are unaffected.

To check, make a cube file of the same orbital with the program's own tool
(`cubegen`, `orca_plot`) and open that. It is drawn from the grid as it is,
with no basis set involved. Then please
[open an issue](https://github.com/lucienwb/microscope/issues) with the Molden
file and the name of the program that wrote it.

### The Lewis structure's charges are not the file's

The Lewis view writes a note at the bottom of the page when the formal charges
in the drawing do not add up to the charge the calculation states. Bond orders
and charges are perceived from the geometry alone, so a strained or unusual
structure can be drawn as a different resonance form, or with the charge on a
different atom. Around a metal, the ligands are counted as ions, so the charge
drawn on the metal is its oxidation state, and the page says that too. A file
with no hydrogens (most crystal structures) gets no formal charges at all,
since the geometry cannot settle them. [How it works](about/how-it-works.md#lewis-perception)
explains the perception.

## Still stuck?

[Open an issue](https://github.com/lucienwb/microscope/issues) with:

- `scope --version`
- your operating system
- what the OpenGL check above prints
- the file, if the problem is with one file
