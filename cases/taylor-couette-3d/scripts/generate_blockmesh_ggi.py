#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [ "foamlib" ]
# ///
# DISCLAIMER: script was generated mostly by AI agents

import math, sys, os
from datetime import datetime
from typing import Tuple

sys.path.insert(0, os.path.dirname(__file__))
from case_config import *

n_radial_inner = int(10*resolution) # Radial cells in inner zone
n_radial_outer = int(10*resolution)   # Radial cells in outer zone
n_circumferential = int(32*resolution)  # Circumferential cells per sector
n_axial = int(20*resolution)          # Axial cells
n_sectors = 8         # Number of sectors (O-grid)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def cyl_to_cart(r: float, theta_deg: float, z: float) -> Tuple[float, float, float]:
    """Convert cylindrical to Cartesian coordinates."""
    theta_rad = math.radians(theta_deg)
    x = r * math.cos(theta_rad)
    y = r * math.sin(theta_rad)
    return (x, y, z)

def format_vertex(x: float, y: float, z: float) -> str:
    """Format vertex as blockMeshDict entry."""
    return f"({x:.8f} {y:.8f} {z:.8f})"

def arc_point(r: float, theta_deg: float, z: float) -> Tuple[float, float, float]:
    """Calculate arc midpoint at 22.5° offset from theta."""
    theta_rad = math.radians(theta_deg + 22.5)  # Midpoint between sectors
    x = r * math.cos(theta_rad)
    y = r * math.sin(theta_rad)
    return (x, y, z)

# ============================================================================
# VERTEX GENERATION
# ============================================================================

# Sector angles
angles = [i * 360 / n_sectors for i in range(n_sectors)]  # [0, 45, 90, 135, 180, 225, 270, 315]

vertices = []
vertex_comments = []

# INNER RING VERTICES (0-31)
# Inner surface (r_i): 16 vertices (8 sectors × 2 axial levels)
for z in [0.0, H]:
    for angle in angles:
        x, y, z_coord = cyl_to_cart(r_i, angle, z)
        vertices.append((x, y, z_coord))
        vertex_comments.append(f"Inner ring, r={r_i}, θ={angle}°, z={z}")

# Outer surface of inner ring (r_interface): 16 vertices (8 sectors × 2 axial levels)
for z in [0.0, H]:
    for angle in angles:
        x, y, z_coord = cyl_to_cart(r_interface, angle, z)
        vertices.append((x, y, z_coord))
        vertex_comments.append(f"Inner ring, r={r_interface}, θ={angle}°, z={z}")

# OUTER RING VERTICES (32-63)
# Inner surface of outer ring (r_interface): 16 vertices (8 sectors × 2 axial levels)
# These are DUPLICATES of inner ring outer surface (disconnected topology!)
for z in [0.0, H]:
    for angle in angles:
        x, y, z_coord = cyl_to_cart(r_interface, angle, z)
        vertices.append((x, y, z_coord))
        vertex_comments.append(f"Outer ring, r={r_interface}, θ={angle}°, z={z} [DUP]")

# Outer surface (r_o): 16 vertices (8 sectors × 2 axial levels)
for z in [0.0, H]:
    for angle in angles:
        x, y, z_coord = cyl_to_cart(r_o, angle, z)
        vertices.append((x, y, z_coord))
        vertex_comments.append(f"Outer ring, r={r_o}, θ={angle}°, z={z}")

print(f"\nGenerated {len(vertices)} vertices")

# ============================================================================
# BLOCK GENERATION
# ============================================================================

blocks = []

# Inner ring blocks (0-7)
for i in range(n_sectors):
    i_next = (i + 1) % n_sectors
    # Correct hex ordering: bottom face goes inner→outer→outer→inner
    # Then top face has same pattern
    v0 = i           # inner, bottom, angle i
    v1 = 16 + i      # outer, bottom, angle i
    v2 = 16 + i_next # outer, bottom, angle i+1
    v3 = i_next      # inner, bottom, angle i+1
    v4 = 8 + i       # inner, top, angle i
    v5 = 24 + i      # outer, top, angle i
    v6 = 24 + i_next # outer, top, angle i+1
    v7 = 8 + i_next  # inner, top, angle i+1

    blocks.append({
        'vertices': (v0, v1, v2, v3, v4, v5, v6, v7),
        'cells': (n_radial_inner, n_circumferential // n_sectors, n_axial),
        'grading': '(1 1 1)',
        'comment': f'Inner ring sector {i}'
    })

# Outer ring blocks (8-15)
for i in range(n_sectors):
    i_next = (i + 1) % n_sectors
    # Correct hex ordering: bottom face goes inner→outer→outer→inner
    v0 = 32 + i           # inner, bottom, angle i
    v1 = 48 + i           # outer, bottom, angle i
    v2 = 48 + i_next      # outer, bottom, angle i+1
    v3 = 32 + i_next      # inner, bottom, angle i+1
    v4 = 40 + i           # inner, top, angle i
    v5 = 56 + i           # outer, top, angle i
    v6 = 56 + i_next      # outer, top, angle i+1
    v7 = 40 + i_next      # inner, top, angle i+1

    blocks.append({
        'vertices': (v0, v1, v2, v3, v4, v5, v6, v7),
        'cells': (n_radial_outer, n_circumferential // n_sectors, n_axial),
        'grading': '(1 1 1)',
        'comment': f'Outer ring sector {i}'
    })

print(f"Generated {len(blocks)} blocks")

# ============================================================================
# EDGE GENERATION (ARCS)
# ============================================================================

edges = []

# Inner ring arcs
for i in range(n_sectors):
    i_next = (i + 1) % n_sectors

    # Bottom level arcs
    for r, v_base in [(r_i, 0), (r_interface, 16)]:
        v1 = v_base + i
        v2 = v_base + i_next
        arc_mid = arc_point(r, angles[i], 0.0)
        edges.append(f"arc {v1} {v2} {format_vertex(*arc_mid)}")

    # Top level arcs
    for r, v_base in [(r_i, 8), (r_interface, 24)]:
        v1 = v_base + i
        v2 = v_base + i_next
        arc_mid = arc_point(r, angles[i], H)
        edges.append(f"arc {v1} {v2} {format_vertex(*arc_mid)}")

# Outer ring arcs
for i in range(n_sectors):
    i_next = (i + 1) % n_sectors

    # Bottom level arcs
    for r, v_base in [(r_interface, 32), (r_o, 48)]:
        v1 = v_base + i
        v2 = v_base + i_next
        arc_mid = arc_point(r, angles[i], 0.0)
        edges.append(f"arc {v1} {v2} {format_vertex(*arc_mid)}")

    # Top level arcs
    for r, v_base in [(r_interface, 40), (r_o, 56)]:
        v1 = v_base + i
        v2 = v_base + i_next
        arc_mid = arc_point(r, angles[i], H)
        edges.append(f"arc {v1} {v2} {format_vertex(*arc_mid)}")

print(f"Generated {len(edges)} edge arcs")

# ============================================================================
# BOUNDARY PATCHES
# ============================================================================

# Helper to generate face list for a boundary patch
def get_bottom_faces():
    """Bottom faces (z=0) of all blocks."""
    faces = []
    # Inner ring: bottom face is (v0, v1, v2, v3) = (inner_i, outer_i, outer_i+1, inner_i+1)
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({i} {16+i} {16+i_next} {i_next})")
    # Outer ring
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({32+i} {48+i} {48+i_next} {32+i_next})")
    return faces

def get_top_faces():
    """Top faces (z=H) of all blocks."""
    faces = []
    # Inner ring: top face is (v4, v5, v6, v7) = (inner_i, outer_i, outer_i+1, inner_i+1) at top
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({8+i} {24+i} {24+i_next} {8+i_next})")
    # Outer ring
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({40+i} {56+i} {56+i_next} {40+i_next})")
    return faces

def get_inner_wall_faces():
    """Inner cylinder wall (r=r_i)."""
    faces = []
    # Connects inner vertices: (bottom_i, bottom_i+1, top_i+1, top_i)
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({i} {i_next} {8+i_next} {8+i})")
    return faces

def get_outer_wall_faces():
    """Outer cylinder wall (r=r_o)."""
    faces = []
    # Connects outer vertices: (bottom_i, bottom_i+1, top_i+1, top_i)
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({48+i} {48+i_next} {56+i_next} {56+i})")
    return faces

def get_ggi_inner_faces():
    """GGI inner interface (outer surface of inner ring at r=r_interface)."""
    faces = []
    # Outer surface of inner ring: (outer_bottom_i, outer_bottom_i+1, outer_top_i+1, outer_top_i)
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({16+i} {16+i_next} {24+i_next} {24+i})")
    return faces

def get_ggi_outer_faces():
    """GGI outer interface (inner surface of outer ring at r=r_interface)."""
    faces = []
    # Inner surface of outer ring: (inner_bottom_i, inner_bottom_i+1, inner_top_i+1, inner_top_i)
    for i in range(n_sectors):
        i_next = (i + 1) % n_sectors
        faces.append(f"({32+i} {32+i_next} {40+i_next} {40+i})")
    return faces

# ============================================================================
# WRITE BLOCKMESHDICT
# ============================================================================

output_path = "system/blockMeshDict"

print(f"\nWriting blockMeshDict to: {output_path}")

with open(output_path, 'w') as f:
    # Header
    f.write(f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\\    /   O peration     | Version:     5.0                                |
|   \\\\  /    A nd           | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}            |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Taylor-Couette 3D mesh with GGI interfaces
// Generated by: generate_blockmesh_ggi.py
//
// Geometry:
//   Inner radius (r_i):      {r_i} m
//   Interface radius:        {r_interface} m
//   Outer radius (r_o):      {r_o} m
//   Height (H):              {H} m
//
// Mesh parameters:
//   Radial (inner):          {n_radial_inner} cells
//   Radial (outer):          {n_radial_outer} cells
//   Circumferential:         {n_circumferential} cells (per sector: {n_circumferential // n_sectors})
//   Axial:                   {n_axial} cells
//   Sectors:                 {n_sectors}
//
// Topology:
//   - Inner ring: vertices 0-31, blocks 0-7
//   - Outer ring: vertices 32-63, blocks 8-15
//   - Interface vertices are DUPLICATED (disconnected topology for GGI)

convertToMeters 1;

vertices
(
""")

    # Write vertices
    for i, (v, comment) in enumerate(zip(vertices, vertex_comments)):
        f.write(f"    {format_vertex(*v)} // {i} - {comment}\n")

    f.write(");\n\n")

    # Write blocks
    f.write("blocks\n(\n")
    for i, block in enumerate(blocks):
        verts = ' '.join(map(str, block['vertices']))
        cells = ' '.join(map(str, block['cells']))
        f.write(f"    hex ({verts}) ({cells}) simpleGrading {block['grading']} // {i} - {block['comment']}\n")
    f.write(");\n\n")

    # Write edges
    f.write("edges\n(\n")
    for edge in edges:
        f.write(f"    {edge}\n")
    f.write(");\n\n")

    # Write boundary
    f.write("boundary\n(\n")

    # Bottom patch
    f.write("    bottom\n    {\n")
    f.write("        type patch;\n")
    f.write("        faces\n        (\n")
    for face in get_bottom_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n\n")

    # Top patch
    f.write("    top\n    {\n")
    f.write("        type patch;\n")
    f.write("        faces\n        (\n")
    for face in get_top_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n\n")

    # Inner wall
    f.write("    innerWall\n    {\n")
    f.write("        type wall;\n")
    f.write("        faces\n        (\n")
    for face in get_inner_wall_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n\n")

    # Outer wall
    f.write("    outerWall\n    {\n")
    f.write("        type wall;\n")
    f.write("        faces\n        (\n")
    for face in get_outer_wall_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n\n")

    # GGI inner patch
    f.write("    GGI_inner\n    {\n")
    f.write("        type            ggi;\n")
    f.write("        shadowPatch     GGI_outer;\n")
    f.write("        zone            GGI_inner_zone;\n")
    f.write("        bridgeOverlap   false;\n")
    f.write("        rotationAxis    (0 0 1);\n")
    f.write("        rotationAngle   0;\n")
    f.write("        separationOffset (0 0 0);\n")
    f.write("        faces\n        (\n")
    for face in get_ggi_inner_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n\n")

    # GGI outer patch
    f.write("    GGI_outer\n    {\n")
    f.write("        type            ggi;\n")
    f.write("        shadowPatch     GGI_inner;\n")
    f.write("        zone            GGI_outer_zone;\n")
    f.write("        bridgeOverlap   false;\n")
    f.write("        rotationAxis    (0 0 1);\n")
    f.write("        rotationAngle   0;\n")
    f.write("        separationOffset (0 0 0);\n")
    f.write("        faces\n        (\n")
    for face in get_ggi_outer_faces():
        f.write(f"            {face}\n")
    f.write("        );\n    }\n")

    f.write(");\n\n")

    # Empty sections
    f.write("mergePatchPairs\n(\n);\n\n")

    f.write("// ************************************************************************* //\n")

print(f"✓ BlockMeshDict written to {output_path}")
