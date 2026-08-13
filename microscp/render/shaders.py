"""GLSL sources: ray-cast impostor spheres and cylinders, orthographic camera.

Each primitive is drawn as an instanced screen-aligned quad; the fragment
shader computes the exact surface point, normal, and depth per pixel, which
gives perfectly smooth geometry at any zoom (the CYLview look).
"""

_LIGHTING = """
uniform vec3 uLightDir;   // view space, normalized

vec3 shade(vec3 n, vec3 base) {
    vec3 l = normalize(uLightDir);
    float wrap = clamp((dot(n, l) + 0.5) / 1.5, 0.0, 1.0);
    vec3 col = base * (0.38 + 0.62 * wrap);
    vec3 h = normalize(l + vec3(0.0, 0.0, 1.0));
    float spec = pow(max(dot(n, h), 0.0), 90.0) * 0.28;
    col += vec3(spec);
    col *= 0.86 + 0.14 * n.z;   // soft silhouette falloff
    return col;
}
"""

SPHERE_VERT = """#version 330 core
layout(location = 0) in vec2 corner;
layout(location = 1) in vec3 iCenter;
layout(location = 2) in float iRadius;
layout(location = 3) in vec3 iColor;
layout(location = 4) in float iQuad;
uniform mat4 uView;
uniform mat4 uProj;
out vec2 vCorner;
flat out vec3 vColor;
flat out vec3 vCenter;
flat out float vRadius;
flat out float vQuad;
flat out int vId;

void main() {
    vec3 c = (uView * vec4(iCenter, 1.0)).xyz;
    vCorner = corner;
    vColor = iColor;
    vCenter = c;
    vRadius = iRadius;
    vQuad = iQuad;
    vId = gl_InstanceID;
    vec3 pos = c + vec3(corner * iRadius, 0.0);
    gl_Position = uProj * vec4(pos, 1.0);
}
"""

SPHERE_FRAG = """#version 330 core
in vec2 vCorner;
flat in vec3 vColor;
flat in vec3 vCenter;
flat in float vRadius;
flat in float vQuad;
flat in int vId;
uniform mat4 uProj;
uniform int uPick;
uniform vec3 uQuadColor;
uniform float uQuadWidth; // 0 = no seam lines
uniform vec3 uQuadA;      // view-space seam plane normals (world axes,
uniform vec3 uQuadB;      // so the seams rotate with the molecule)
out vec4 fragColor;
""" + _LIGHTING + """
void main() {
    float r2 = dot(vCorner, vCorner);
    if (r2 > 1.0) discard;
    float nz = sqrt(1.0 - r2);
    vec3 n = vec3(vCorner, nz);
    float zView = vCenter.z + vRadius * nz;
    float ndcZ = uProj[2][2] * zView + uProj[3][2];
    gl_FragDepth = ndcZ * 0.5 + 0.5;
    if (uPick == 1) {
        int id = vId + 1;
        fragColor = vec4(float(id & 0xFF) / 255.0,
                         float((id >> 8) & 0xFF) / 255.0,
                         float((id >> 16) & 0xFF) / 255.0, 1.0);
        return;
    }
    vec3 base = vColor;
    // Houkmol "quadrants": two orthogonal great-circle seam lines
    if (vQuad > 0.5 && uQuadWidth > 0.0 &&
        (abs(dot(n, uQuadA)) < uQuadWidth || abs(dot(n, uQuadB)) < uQuadWidth))
        base = uQuadColor;
    fragColor = vec4(shade(n, base), 1.0);
}
"""

CYLINDER_VERT = """#version 330 core
layout(location = 0) in vec2 corner;
layout(location = 1) in vec3 iA;
layout(location = 2) in vec3 iB;
layout(location = 3) in float iRadius;
layout(location = 4) in vec3 iColorA;
layout(location = 5) in vec3 iColorB;
uniform mat4 uView;
uniform mat4 uProj;
out vec2 vPos;
flat out vec3 vA;
flat out vec3 vB;
flat out float vRadius;
flat out vec3 vColorA;
flat out vec3 vColorB;

void main() {
    vec3 a = (uView * vec4(iA, 1.0)).xyz;
    vec3 b = (uView * vec4(iB, 1.0)).xyz;
    vA = a;
    vB = b;
    vRadius = iRadius;
    vColorA = iColorA;
    vColorB = iColorB;
    vec2 axis = b.xy - a.xy;
    float len = length(axis);
    vec2 p = (len > 1e-6) ? axis / len : vec2(1.0, 0.0);
    vec2 q = vec2(-p.y, p.x);
    float s = corner.x * 0.5 + 0.5;
    vec2 base = mix(a.xy, b.xy, s) + p * (iRadius * corner.x) + q * (iRadius * corner.y);
    float z = max(a.z, b.z) + iRadius;
    vPos = base;
    gl_Position = uProj * vec4(base, z, 1.0);
}
"""

CYLINDER_FRAG = """#version 330 core
in vec2 vPos;
flat in vec3 vA;
flat in vec3 vB;
flat in float vRadius;
flat in vec3 vColorA;
flat in vec3 vColorB;
uniform mat4 uProj;
uniform int uPick;
out vec4 fragColor;
""" + _LIGHTING + """
void main() {
    vec3 axisVec = vB - vA;
    float L = length(axisVec);
    if (L < 1e-6) discard;
    vec3 d = axisVec / L;
    vec3 F = vec3(vPos, 0.0);       // fragment ray origin (orthographic)
    vec3 e = vec3(0.0, 0.0, -1.0);  // ray direction
    vec3 u = F - vA;
    float ed = dot(e, d);
    float A = 1.0 - ed * ed;
    if (A < 1e-5) discard;          // end-on: joint spheres cover this case
    float ud = dot(u, d);
    float B = 2.0 * (dot(u, e) - ud * ed);
    float C = dot(u, u) - ud * ud - vRadius * vRadius;
    float disc = B * B - 4.0 * A * C;
    if (disc < 0.0) discard;
    float s = (-B - sqrt(disc)) / (2.0 * A);   // nearest intersection
    vec3 H = F + s * e;
    float t = dot(H - vA, d);
    if (t < 0.0 || t > L) discard;
    vec3 n = (H - (vA + t * d)) / vRadius;
    float ndcZ = uProj[2][2] * H.z + uProj[3][2];
    gl_FragDepth = ndcZ * 0.5 + 0.5;
    if (uPick == 1) {
        fragColor = vec4(0.0, 0.0, 0.0, 1.0);  // bonds occlude but are not pickable
        return;
    }
    vec3 base = (t < 0.5 * L) ? vColorA : vColorB;
    fragColor = vec4(shade(n, base), 1.0);
}
"""

MESH_VERT = """#version 330 core
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
uniform mat4 uView;
uniform mat4 uProj;
out vec3 vNormal;

void main() {
    vNormal = mat3(uView) * normal;
    gl_Position = uProj * (uView * vec4(position, 1.0));
}
"""

MESH_FRAG = """#version 330 core
in vec3 vNormal;
uniform vec4 uColor;    // rgb + opacity
out vec4 fragColor;
""" + _LIGHTING + """
void main() {
    vec3 n = normalize(vNormal);
    if (n.z < 0.0) n = -n;   // two-sided: light whichever face shows
    fragColor = vec4(shade(n, uColor.rgb), uColor.a);
}
"""
