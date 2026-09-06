from __future__ import annotations

import json
import struct
from pathlib import Path
from typing import Any, Dict, Union

import aerosandbox as asb
import numpy as np

PathLike = Union[str, Path]


def export_airplane_glb(airplane: asb.Airplane, out_path: PathLike) -> Path:
    """Write the design airplane to a binary glTF (.glb) + a .meta.json sidecar.

    This is the AUTHORITATIVE CAD/geometry artifact: it meshes the real solved
    airplane (twin fins, twin pods, twin booms), so any consumer (the sibling sim's
    Commandant viewer) renders exactly the run's airframe instead of rebuilding a
    surrogate. The GLB writer is stdlib-only (struct+json), so no extra dependency.
    It emits POSITION + NORMAL + a PBR material: without normals/material a renderer
    like Cesium can't shade the mesh and draws NOTHING (only the orientation axes) --
    which is exactly the "no aircraft in Commandant" symptom.

    Frame: emit the aircraft in the plain AeroSandbox BODY frame -- X=aft, Y=starboard,
    Z=up. Commandant applies a fixed HeadingPitchRoll(0,90,180) stand-up correction that
    is authored for exactly this frame, so the plane renders level (top up, nose/wings
    horizontal). (An earlier Cesium-style remap made the correction stand it on its tail.)
    """
    points, faces = airplane.mesh_body(method="tri")
    points = np.ascontiguousarray(np.asarray(points, dtype=np.float32))
    faces = np.ascontiguousarray(np.asarray(faces, dtype=np.uint32).reshape(-1, 3))

    # Smooth per-vertex normals: sum each triangle's face normal onto its 3 vertices,
    # then normalize. Required for the renderer to shade the surface.
    normals = np.zeros(points.shape, dtype=np.float64)
    tris = points[faces]
    fn = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    for k in range(3):
        np.add.at(normals, faces[:, k], fn)
    ln = np.linalg.norm(normals, axis=1, keepdims=True)
    ln[ln == 0] = 1.0
    normals = np.ascontiguousarray(normals / ln, dtype=np.float32)
    idx = np.ascontiguousarray(faces.reshape(-1), dtype=np.uint32)

    pos_bytes, nrm_bytes, idx_bytes = points.tobytes(), normals.tobytes(), idx.tobytes()
    nrm_off, idx_off = len(pos_bytes), len(pos_bytes) + len(nrm_bytes)
    pad = lambda buf, fill: buf + fill * ((-len(buf)) % 4)   # 4-byte chunk align
    bin_blob = pad(pos_bytes + nrm_bytes + idx_bytes, b"\x00")
    mins, maxs = points.min(axis=0).tolist(), points.max(axis=0).tolist()

    gltf = {
        "asset": {"version": "2.0", "generator": "solar-airplane/geometry"},
        "scene": 0, "scenes": [{"nodes": [0]}],
        "nodes": [{"mesh": 0, "name": str(airplane.name)}],
        "meshes": [{"name": str(airplane.name), "primitives": [
            {"attributes": {"POSITION": 0, "NORMAL": 1}, "indices": 2, "material": 0, "mode": 4}]}],
        "materials": [{
            "name": "airframe",
            # Dark solar-slate, not white: a light base reads as a flat white blob
            # against terrain. Some metallic + mid roughness gives visible 3-D shading.
            "pbrMetallicRoughness": {"baseColorFactor": [0.20, 0.24, 0.32, 1.0],
                                     "metallicFactor": 0.35, "roughnessFactor": 0.5},
            "doubleSided": True,
        }],
        "buffers": [{"byteLength": len(bin_blob)}],
        "bufferViews": [
            {"buffer": 0, "byteOffset": 0, "byteLength": len(pos_bytes), "target": 34962},
            {"buffer": 0, "byteOffset": nrm_off, "byteLength": len(nrm_bytes), "target": 34962},
            {"buffer": 0, "byteOffset": idx_off, "byteLength": len(idx_bytes), "target": 34963},
        ],
        "accessors": [
            {"bufferView": 0, "componentType": 5126, "count": int(points.shape[0]),
             "type": "VEC3", "min": mins, "max": maxs},
            {"bufferView": 1, "componentType": 5126, "count": int(normals.shape[0]), "type": "VEC3"},
            {"bufferView": 2, "componentType": 5125, "count": int(idx.shape[0]), "type": "SCALAR"},
        ],
    }
    json_blob = pad(json.dumps(gltf, separators=(",", ":")).encode("utf-8"), b" ")

    out = Path(out_path).with_suffix(".glb")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        total = 12 + 8 + len(json_blob) + 8 + len(bin_blob)
        f.write(struct.pack("<III", 0x46546C67, 2, total))       # header
        f.write(struct.pack("<II", len(json_blob), 0x4E4F534A))  # JSON chunk
        f.write(json_blob)
        f.write(struct.pack("<II", len(bin_blob), 0x004E4942))   # BIN chunk
        f.write(bin_blob)

    out.with_suffix(".meta.json").write_text(json.dumps({
        "name": str(airplane.name),
        "frame": {"up": "+Z", "forward": "-X", "right": "+Y",
                  "note": "AeroSandbox body frame (x-aft, y-starboard, z-up); Commandant's "
                          "HeadingPitchRoll(0,90,180) stands it up level."},
        "reference": {"span_m": float(points[:, 1].max() - points[:, 1].min()),   # Y = starboard
                      "length_m": float(points[:, 0].max() - points[:, 0].min()),  # X = fore/aft
                      "vertices": int(points.shape[0]), "triangles": int(faces.shape[0])},
    }, indent=2))
    return out


def export_xflr5_xml_from_soln(
    *,
    airplane_sol: asb.Airplane,
    soln: Dict[str, Any],
    out_path: PathLike,
    include_fuselages: bool = False,
) -> None:
    """Exports an XFLR5 XML from an already-solved airplane, using parameters stored in `soln`.

    AeroSandbox's XFLR5 exporter supports exactly 3 lifting surfaces: main wing, elevator, and fin.
    This project has a twin-fin configuration, so we synthesize a single centerline fin using
    values from the grouped output sections (e.g. `soln['Geometry']`, `soln['V Stab']`, `soln['HStab']`).
    """
    out_path = Path(out_path)

    geom = soln["Geometry"]
    vstab = soln["V Stab"]
    hstab = soln["HStab"]

    boom_length = float(geom["boom_length"])
    vstab_root_chord = float(vstab["vstab_root_chord"])
    vstab_span = float(vstab["vstab_span"])
    vstab_tip_chord = float(hstab["hstab_chordlen"])

    tail_airfoil_name = soln.get("Airfoils", {}).get("tail_airfoil", "naca0010")
    tail_airfoil = asb.Airfoil(tail_airfoil_name)

    if len(getattr(airplane_sol, "wings", [])) < 2:
        raise ValueError("airplane_sol must have at least 2 wings (main wing, elevator).")

    mainwing_xflr = airplane_sol.wings[0]
    elevator_xflr = airplane_sol.wings[1]

    fin_xflr = asb.Wing(
        name="Vertical Stabilizer (XFLR5)",
        symmetric=False,
        xsecs=[
            asb.WingXSec(
                xyz_le=[0.0, 0.0, 0.0],
                chord=vstab_root_chord,
                twist=0.0,
                airfoil=tail_airfoil,
            ),
            asb.WingXSec(
                xyz_le=[vstab_root_chord / 4, 0.0, vstab_span],
                chord=vstab_tip_chord,
                twist=0.0,
                airfoil=tail_airfoil,
            ),
        ],
    ).translate([boom_length, 0.0, 0.0])

    airplane_sol.export_XFLR5_xml(
        str(out_path),
        include_fuselages=include_fuselages,
        mainwing=mainwing_xflr,
        elevator=elevator_xflr,
        fin=fin_xflr,
    )


def export_cadquery_step(
    *,
    airplane_sol: asb.Airplane,
    out_path: PathLike,
    minimum_airfoil_TE_thickness: float = 0.001,
) -> None:
    """Exports a CadQuery STEP file from an already-solved airplane."""
    out_path = Path(out_path)
    airplane_sol.export_cadquery_geometry(
        str(out_path),
        minimum_airfoil_TE_thickness=float(minimum_airfoil_TE_thickness),
    )

