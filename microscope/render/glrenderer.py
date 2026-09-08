"""OpenGL renderer: instanced ray-cast impostor spheres and cylinders."""

from __future__ import annotations

import ctypes

import numpy as np
from OpenGL import GL

from .scene import SceneBuffers
from .shaders import CYLINDER_FRAG, CYLINDER_VERT, MESH_FRAG, MESH_VERT, SPHERE_FRAG, SPHERE_VERT


class RendererError(RuntimeError):
    pass


def _compile_shader(source: str, kind) -> int:
    shader = GL.glCreateShader(kind)
    GL.glShaderSource(shader, source)
    GL.glCompileShader(shader)
    if not GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS):
        log = GL.glGetShaderInfoLog(shader).decode(errors="replace")
        raise RendererError(f"shader compile failed:\n{log}")
    return shader


def _link_program(vert_src: str, frag_src: str) -> int:
    vs = _compile_shader(vert_src, GL.GL_VERTEX_SHADER)
    fs = _compile_shader(frag_src, GL.GL_FRAGMENT_SHADER)
    prog = GL.glCreateProgram()
    GL.glAttachShader(prog, vs)
    GL.glAttachShader(prog, fs)
    GL.glLinkProgram(prog)
    GL.glDeleteShader(vs)
    GL.glDeleteShader(fs)
    if not GL.glGetProgramiv(prog, GL.GL_LINK_STATUS):
        log = GL.glGetProgramInfoLog(prog).decode(errors="replace")
        raise RendererError(f"program link failed:\n{log}")
    return prog


def create_fbo(width: int, height: int):
    """Create and bind an RGBA8 + depth24 framebuffer; returns (fbo, color, depth)."""
    fbo = GL.glGenFramebuffers(1)
    color = GL.glGenRenderbuffers(1)
    depth = GL.glGenRenderbuffers(1)
    GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, color)
    GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_RGBA8, width, height)
    GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, depth)
    GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT24, width, height)
    GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo)
    GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                 GL.GL_RENDERBUFFER, color)
    GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                 GL.GL_RENDERBUFFER, depth)
    if GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER) != GL.GL_FRAMEBUFFER_COMPLETE:
        raise RendererError("framebuffer incomplete")
    return fbo, color, depth


def delete_fbo(fbo, color, depth) -> None:
    GL.glDeleteFramebuffers(1, [fbo])
    GL.glDeleteRenderbuffers(1, [color])
    GL.glDeleteRenderbuffers(1, [depth])


class MoleculeRenderer:
    SPHERE_STRIDE = 8 * 4
    CYLINDER_STRIDE = 13 * 4

    def __init__(self):
        self._ready = False
        self._nspheres = 0
        self._ncylinders = 0
        self._supports_sample_shading = False
        self._meshes: list[tuple[int, int, int, tuple]] = []  # vao, vbo, nverts, rgba
        self._quad_color = (0.0, 0.0, 0.0)
        self._quad_width = 0.0

    def set_style_params(self, quad_color: tuple | None = None,
                         quad_width: float = 0.0) -> None:
        """Houkmol seam-line ("quadrant") parameters from the Style."""
        self._quad_color = tuple(quad_color) if quad_color else (0.0, 0.0, 0.0)
        self._quad_width = float(quad_width) if quad_color else 0.0

    def initialize(self) -> None:
        self._sphere_prog = _link_program(SPHERE_VERT, SPHERE_FRAG)
        self._cyl_prog = _link_program(CYLINDER_VERT, CYLINDER_FRAG)
        self._mesh_prog = _link_program(MESH_VERT, MESH_FRAG)

        quad = np.array([-1, -1, 1, -1, -1, 1, 1, 1], dtype=np.float32)
        self._quad_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._quad_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, quad.nbytes, quad, GL.GL_STATIC_DRAW)
        self._sphere_vbo = GL.glGenBuffers(1)
        self._cyl_vbo = GL.glGenBuffers(1)

        self._sphere_vao = self._make_vao(
            self._sphere_vbo, self.SPHERE_STRIDE,
            ((1, 3, 0), (2, 1, 12), (3, 3, 16), (4, 1, 28)),
        )
        self._cyl_vao = self._make_vao(
            self._cyl_vbo, self.CYLINDER_STRIDE,
            ((1, 3, 0), (2, 3, 12), (3, 1, 24), (4, 3, 28), (5, 3, 40)),
        )

        self._supports_sample_shading = GL.glGetIntegerv(GL.GL_MAJOR_VERSION) >= 4
        self._ready = True

    def _make_vao(self, instance_vbo, stride, attribs) -> int:
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._quad_vbo)
        GL.glEnableVertexAttribArray(0)
        GL.glVertexAttribPointer(0, 2, GL.GL_FLOAT, GL.GL_FALSE, 8, ctypes.c_void_p(0))
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, instance_vbo)
        for loc, size, offset in attribs:
            GL.glEnableVertexAttribArray(loc)
            GL.glVertexAttribPointer(loc, size, GL.GL_FLOAT, GL.GL_FALSE,
                                     stride, ctypes.c_void_p(offset))
            GL.glVertexAttribDivisor(loc, 1)
        GL.glBindVertexArray(0)
        return vao

    def set_scene(self, scene: SceneBuffers) -> None:
        self._nspheres = len(scene.spheres)
        self._ncylinders = len(scene.cylinders)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._sphere_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, scene.spheres.nbytes,
                        scene.spheres, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._cyl_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, scene.cylinders.nbytes,
                        scene.cylinders, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

    def set_meshes(self, meshes) -> None:
        """Upload translucent triangle meshes: (vertices, normals, rgba) tuples.

        Vertices/normals are (N, 3) float32 with consecutive vertex triples
        forming triangles. Requires a current GL context (like set_scene).
        """
        for vao, vbo, _count, _color in self._meshes:
            GL.glDeleteVertexArrays(1, [vao])
            GL.glDeleteBuffers(1, [vbo])
        self._meshes = []
        for verts, normals, rgba in meshes:
            if len(verts) == 0:
                continue
            data = np.hstack([verts, normals]).astype(np.float32)
            vbo = GL.glGenBuffers(1)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, data.nbytes, data, GL.GL_STATIC_DRAW)
            vao = GL.glGenVertexArrays(1)
            GL.glBindVertexArray(vao)
            for loc, offset in ((0, 0), (1, 12)):
                GL.glEnableVertexAttribArray(loc)
                GL.glVertexAttribPointer(loc, 3, GL.GL_FLOAT, GL.GL_FALSE,
                                         24, ctypes.c_void_p(offset))
            GL.glBindVertexArray(0)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
            self._meshes.append((vao, vbo, len(verts), tuple(rgba)))

    def draw(self, view: np.ndarray, proj: np.ndarray, width: int, height: int,
             background=(1.0, 1.0, 1.0, 1.0), pick: bool = False,
             sample_shading: bool = True) -> None:
        GL.glViewport(0, 0, int(width), int(height))
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthFunc(GL.GL_LESS)
        GL.glDisable(GL.GL_BLEND)
        GL.glClearColor(*background)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        use_ss = sample_shading and not pick and self._supports_sample_shading
        if use_ss:
            GL.glEnable(GL.GL_SAMPLE_SHADING)
            GL.glMinSampleShading(1.0)

        light = np.array([-0.35, 0.55, 0.80], dtype=np.float32)
        light /= np.linalg.norm(light)
        vm = np.ascontiguousarray(view, dtype=np.float32)
        pm = np.ascontiguousarray(proj, dtype=np.float32)

        # Houkmol seam planes: world x/y axes in view space, so the seam
        # lines rotate with the molecule
        quad_a = np.ascontiguousarray(vm[:3, 0])
        quad_b = np.ascontiguousarray(vm[:3, 1])

        batches = ((self._sphere_prog, self._sphere_vao, self._nspheres),
                   (self._cyl_prog, self._cyl_vao, self._ncylinders))
        for prog, vao, count in batches:
            if count == 0:
                continue
            GL.glUseProgram(prog)
            GL.glUniformMatrix4fv(GL.glGetUniformLocation(prog, "uView"), 1, GL.GL_TRUE, vm)
            GL.glUniformMatrix4fv(GL.glGetUniformLocation(prog, "uProj"), 1, GL.GL_TRUE, pm)
            GL.glUniform3fv(GL.glGetUniformLocation(prog, "uLightDir"), 1, light)
            GL.glUniform1i(GL.glGetUniformLocation(prog, "uPick"), 1 if pick else 0)
            loc = GL.glGetUniformLocation(prog, "uQuadColor")
            if loc != -1:
                GL.glUniform3f(loc, *self._quad_color)
                GL.glUniform1f(GL.glGetUniformLocation(prog, "uQuadWidth"),
                               self._quad_width)
                GL.glUniform3fv(GL.glGetUniformLocation(prog, "uQuadA"), 1, quad_a)
                GL.glUniform3fv(GL.glGetUniformLocation(prog, "uQuadB"), 1, quad_b)
            GL.glBindVertexArray(vao)
            GL.glDrawArraysInstanced(GL.GL_TRIANGLE_STRIP, 0, 4, count)

        if self._meshes and not pick:
            GL.glUseProgram(self._mesh_prog)
            GL.glUniformMatrix4fv(GL.glGetUniformLocation(self._mesh_prog, "uView"),
                                  1, GL.GL_TRUE, vm)
            GL.glUniformMatrix4fv(GL.glGetUniformLocation(self._mesh_prog, "uProj"),
                                  1, GL.GL_TRUE, pm)
            GL.glUniform3fv(GL.glGetUniformLocation(self._mesh_prog, "uLightDir"),
                            1, light)
            color_loc = GL.glGetUniformLocation(self._mesh_prog, "uColor")
            # translucent surfaces, two passes: prime the depth buffer with the
            # nearest surface layer, then blend exactly that layer once —
            # interior/back faces never bleed through as facet noise
            GL.glColorMask(False, False, False, False)
            for vao, _vbo, count, _rgba in self._meshes:
                GL.glBindVertexArray(vao)
                GL.glDrawArrays(GL.GL_TRIANGLES, 0, count)
            GL.glColorMask(True, True, True, True)
            GL.glDepthFunc(GL.GL_LEQUAL)
            GL.glDepthMask(GL.GL_FALSE)
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFuncSeparate(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA,
                                   GL.GL_ONE, GL.GL_ONE_MINUS_SRC_ALPHA)
            for vao, _vbo, count, rgba in self._meshes:
                GL.glUniform4f(color_loc, *rgba)
                GL.glBindVertexArray(vao)
                GL.glDrawArrays(GL.GL_TRIANGLES, 0, count)
            GL.glDepthMask(GL.GL_TRUE)
            GL.glDepthFunc(GL.GL_LESS)
            GL.glDisable(GL.GL_BLEND)

        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        if use_ss:
            GL.glDisable(GL.GL_SAMPLE_SHADING)

    def pick_atom(self, x: float, y: float, view, proj, width: int, height: int) -> int:
        """Pick at GL pixel coords (origin bottom-left). Returns atom index or -1."""
        if not self._ready or self._nspheres == 0:
            return -1
        prev_fbo = GL.glGetIntegerv(GL.GL_FRAMEBUFFER_BINDING)
        fbo, color, depth = create_fbo(int(width), int(height))
        try:
            self.draw(view, proj, width, height, background=(0, 0, 0, 0),
                      pick=True, sample_shading=False)
            GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
            data = GL.glReadPixels(int(x), int(y), 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
            ident = int(data[0]) | (int(data[1]) << 8) | (int(data[2]) << 16)
            return ident - 1
        finally:
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, int(prev_fbo))
            delete_fbo(fbo, color, depth)
