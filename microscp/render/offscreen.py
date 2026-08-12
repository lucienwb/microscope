"""Offscreen high-resolution rendering (publication PNG export)."""

from __future__ import annotations

import numpy as np
from OpenGL import GL
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QOffscreenSurface, QOpenGLContext, QSurfaceFormat

from ..core.molecule import Molecule
from .camera import OrthoCamera
from .glrenderer import MoleculeRenderer, create_fbo, delete_fbo
from .scene import build_scene
from .styles import Style

MAX_RENDER_DIM = 8192


def render_molecule_image(molecule: Molecule, style: Style, camera: OrthoCamera,
                          width: int, height: int, supersample: int = 3,
                          transparent: bool = False) -> QImage:
    """Render to a QImage of (width, height); supersampled for anti-aliasing."""
    ss = max(1, int(supersample))
    while ss > 1 and max(width, height) * ss > MAX_RENDER_DIM:
        ss -= 1
    w, h = int(width) * ss, int(height) * ss

    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)

    ctx = QOpenGLContext()
    ctx.setFormat(fmt)
    if not ctx.create():
        raise RuntimeError("could not create an OpenGL context")
    surface = QOffscreenSurface()
    surface.setFormat(ctx.format())
    surface.create()
    if not ctx.makeCurrent(surface):
        raise RuntimeError("could not make the OpenGL context current")
    try:
        renderer = MoleculeRenderer()
        renderer.initialize()
        if molecule.bonds is None:
            molecule.perceive_bonds()
        renderer.set_scene(build_scene(molecule, style))
        fbo, color, depth = create_fbo(w, h)
        try:
            bg = (*style.background, 0.0 if transparent else 1.0)
            renderer.draw(camera.view_matrix(), camera.proj_matrix(w / h), w, h,
                          background=bg, sample_shading=False)
            GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
            raw = GL.glReadPixels(0, 0, w, h, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        finally:
            delete_fbo(fbo, color, depth)
        arr = np.frombuffer(raw, dtype=np.uint8).reshape(h, w, 4)[::-1]
        img = QImage(arr.tobytes(), w, h, w * 4, QImage.Format.Format_RGBA8888).copy()
        if ss > 1:
            img = img.scaled(int(width), int(height),
                             Qt.AspectRatioMode.IgnoreAspectRatio,
                             Qt.TransformationMode.SmoothTransformation)
        return img
    finally:
        ctx.doneCurrent()
