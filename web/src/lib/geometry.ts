// Reconstruct renderable geometry for each component from the solution + mass properties.
// Ported from the original Three.js viewer (reconstructGeometry). Dimensions come from the soln;
// centers come from each component's mass-properties CG (so dragged parts move with their CG).

import type { MassProperties } from "./types";

export type Geom =
  | { type: "box"; center: [number, number, number]; Lx: number; Ly: number; Lz: number; mass: number }
  | { type: "cylinder"; center: [number, number, number]; radius: number; length: number; axis: "x"; mass: number }
  | { type: "point"; center: [number, number, number]; radius: number; mass: number }
  // Tapered flat-LE planform (root chord at center, tip chord at the ends). Used for the
  // main wing and the vertical/horizontal stabilizers so taper renders in 3D.
  | {
      type: "wing";
      center: [number, number, number];
      rootChord: number;
      tipChord: number;
      span: number;
      thickness: number;
      axis: "y" | "z"; // span direction: "y" = horizontal surface, "z" = vertical fin
      // Optional spanwise sections (half-span, root->tip) for a lofted surface that shows
      // taper AND polyhedral (z rise). When present, overrides the simple trapezoid.
      sections?: { y: number; z: number; chord: number }[];
      mass: number;
    };

// The airframe dimensions every renderer reads out of soln (top-view planform, side view,
// 3D mesh, and the dimension annotations all destructure from here). One place for the soln
// key names + fallbacks so the four views can't drift apart.
export type GeomParams = {
  b: number; cRoot: number; cTip: number;
  boomLen: number; boomY: number; fuseLen: number;
  hstabSpan: number; hstabChord: number;
  vstabSpan: number; vstabRootChord: number;
  dihedralRad: number;
};

export function parseGeometry(soln: Record<string, any>): GeomParams {
  const mw = soln["Main Wing"] || {};
  const g = soln.Geometry || {};
  const h = soln.HStab || {};
  const v = soln["V Stab"] || {};
  const cRoot = mw.chordlen || 0;
  return {
    b: mw.wingspan || mw.b_w || 0,
    cRoot,
    cTip: mw.chord_tip || cRoot, // rev7 taper; falls back to rectangular
    boomLen: g.boom_length || 0,
    boomY: g.boom_y || 0,
    fuseLen: g.fuselage_length || 0.5, // rev7: read from soln, was hardcoded 0.5
    hstabSpan: h.hstab_span || 0,
    hstabChord: h.hstab_chordlen || 0,
    vstabSpan: v.vstab_span || 0,
    vstabRootChord: v.vstab_root_chord || 0,
    dihedralRad: ((mw.dihedral_angle || 0) * Math.PI) / 180,
  };
}

export function reconstructGeometry(
  soln: Record<string, any>,
  mp: MassProperties,
  fuselageRadius: number,
): Record<string, Geom> {
  const power = soln.Power || {};
  const propulsion = soln.Propulsion || {};
  const airfoils = soln.airfoils || {};

  const {
    b: b_w, cRoot: c_root, cTip: c_tip, boomLen: boom_len, boomY: boom_y,
    vstabSpan: vstab_span, vstabRootChord: vstab_root_chord,
    hstabChord: hstab_chord, hstabSpan: hstab_span, fuseLen: fuselage_length,
    dihedralRad: dih_wing,
  } = parseGeometry(soln);

  const t_wing = (airfoils.wing_t_over_c || 0.12) * c_root;
  const t_tail = (airfoils.tail_t_over_c || 0.1) * Math.max(hstab_chord, 1e-6);
  const boom_radius = (soln.Geometry || {}).boom_radius || 0.01;
  const solar_side = power.solar_panel_side_length || 0.125;
  const batt_box = [0.2, 0.1, 0.03];
  const solar_area = (power.solar_panel_n || 0) * solar_side * solar_side;
  const solar_mean_chord = solar_area / Math.max(b_w, 1e-6);
  const prop_radius = 0.5 * (propulsion.propeller_diameter || 0.4);

  const components = mp.components || {};
  const masses = Object.values(components).map((c) => c.mass || 0);
  const maxMass = Math.max(1e-9, ...(masses.length ? masses : [1e-9]));

  const out: Record<string, Geom> = {};
  for (const [name, comp] of Object.entries(components)) {
    const center: [number, number, number] = [comp.x_cg || 0, comp.y_cg || 0, comp.z_cg || 0];
    const mass = comp.mass || 0;

    if (name === "Main wing") {
      // Half-span sections: root -> inboard break (flat, full chord) -> tip (tapered, raised
      // by polyhedral). dihedral is applied only outboard of the break (boom_y).
      const s = b_w / 2;
      const zTip = Math.sin(dih_wing) * Math.max(s - boom_y, 0);
      const sections = [
        { y: 0, z: 0, chord: c_root },
        { y: boom_y, z: 0, chord: c_root },
        { y: s, z: zTip, chord: c_tip },
      ];
      out[name] = {
        type: "wing", center, axis: "y",
        rootChord: Math.max(c_root, 1e-6), tipChord: Math.max(c_tip, 1e-6),
        span: b_w, thickness: Math.max(t_wing, 0.005), sections, mass,
      };
    } else if (name === "Horizontal stabilizer") {
      // Render as a "wing" (rectangular) so it uses the SAME -0.4*chord LE-offset convention as
      // the vertical fins -- otherwise a box centered on the CG sits ~0.06 m off the fin tops.
      out[name] = {
        type: "wing", center, axis: "y",
        rootChord: Math.max(hstab_chord, 1e-6), tipChord: Math.max(hstab_chord, 1e-6),
        span: hstab_span, thickness: Math.max(t_tail, 0.003), mass,
      };
    } else if (name.includes("Vertical stabilizer")) {
      // Tapered fin: root chord at the boom, tip chord = hstab chord at the top.
      out[name] = {
        type: "wing", center, axis: "z",
        rootChord: Math.max(vstab_root_chord, 1e-6), tipChord: Math.max(hstab_chord, 1e-6),
        span: Math.max(vstab_span, 1e-6), thickness: Math.max(t_tail, 0.003), mass,
      };
    } else if (name === "Boom L" || name === "Boom R") {
      out[name] = { type: "cylinder", center, radius: boom_radius, length: boom_len, axis: "x", mass };
    } else if (name === "Fuselage L" || name === "Fuselage R") {
      out[name] = { type: "cylinder", center, radius: fuselageRadius, length: fuselage_length, axis: "x", mass };
    } else if (name === "Solar cells") {
      // Thin tapered layer over the wing planform (so it doesn't read as a rectangular wing).
      const cov = c_root > 0 ? Math.min(1, solar_mean_chord / c_root) : 0.8;
      out[name] = {
        type: "wing", center, axis: "y",
        rootChord: Math.max(c_root * cov, 0.03), tipChord: Math.max(c_tip * cov, 0.03),
        span: b_w, thickness: 0.002, mass,
      };
    } else if (name === "Batteries") {
      out[name] = { type: "box", center, Lx: batt_box[0], Ly: 2 * boom_y, Lz: batt_box[2], mass };
    } else if (name === "Prop L" || name === "Prop R") {
      out[name] = { type: "cylinder", center, radius: Math.max(prop_radius, 0.01), length: 0.001, axis: "x", mass };
    } else {
      out[name] = { type: "point", center, radius: 0.05 * (mass / maxMass) + 0.02, mass };
    }
  }
  return out;
}

// Solar-cell layout, ported from solar-uav-design's _cell_quads(): pack 125 mm cells in
// columns from the centerline out to the tip; each column carries floor(local_chord / pitch)
// chordwise rows, so taper sheds rows as the chord shrinks. Aileron TE is kept clear outboard.
// The wing fills first; any cells left in the soln's budget spill onto the hstab (which the
// optimizer also treats as panel-able, hstab_panel_frac = 0.6). Cell centers are returned in
// each surface's LOCAL frame (wing: LE at x = -0.4*rootChord; hstab: box-centered), so the
// caller anchors them at that component's center. Also returns the cell side length.
export function solarCells(soln: Record<string, any>): {
  wing: [number, number, number][];
  hstab: [number, number, number][];
  side: number;
} {
  const power = soln.Power || {};
  const hstab = soln.HStab || {};
  const airfoils = soln.airfoils || {};
  const { b, cRoot, cTip, boomY: yBreak, dihedralRad: dih } = parseGeometry(soln);
  const s = b / 2;
  const pitch = 0.13;
  const side = 0.125;
  const x0 = -0.4 * cRoot; // match the lofted wing mesh origin (LE offset)
  // Use the optimizer's EXACT wing/hstab split (Power.panels_on_wing / panels_on_hstab) as the
  // per-surface budgets, so the viewer matches the report. Fall back to all-on-wing for old runs.
  const nWing = Math.round(power.panels_on_wing ?? power.solar_panel_n ?? 0);
  const nHstab = Math.round(power.panels_on_hstab ?? 0);
  let budget = nWing; // wing budget first

  const chordAt = (y: number) =>
    y <= yBreak ? cRoot : cRoot + (cTip - cRoot) * Math.min(Math.max((y - yBreak) / Math.max(s - yBreak, 1e-6), 0), 1);
  const zAt = (y: number) => (y <= yBreak ? 0 : Math.sin(dih) * (y - yBreak)) + 0.012;

  const wing: [number, number, number][] = [];
  for (let col = 0; (col + 0.5) * pitch <= s && budget > 0; col++) {
    const yc = (col + 0.5) * pitch;
    const c = chordAt(yc);
    const rows = Math.floor(c / pitch);
    const z = zAt(yc);
    for (let r = 0; r < rows && budget > 0; r++) {
      // aileron keep-out: skip the aft rows in the outboard aileron span (55%-95% semi-span)
      if (yc > 0.55 * s && yc < 0.95 * s && (r + 1) * pitch > 0.75 * c) continue;
      const xc = x0 + (r + 0.5) * pitch;
      wing.push([xc, yc, z]);   // right wing
      budget--;
      if (budget <= 0) break;
      wing.push([xc, -yc, z]);  // left wing (mirror)
      budget--;
    }
  }

  // Spill the remaining budget onto the hstab. Cells sit on its forward panel-able band: a
  // ~10% LE margin then panels back to ~70% chord (the aft ~30% is the elevator keep-out),
  // matching hstab_panel_frac = 0.6. Wing-local frame: LE at x = -0.4*chord (same convention
  // as the hstab mesh), top surface at +t/2.
  const hChord = hstab.hstab_chordlen || 0;
  const hSpan = hstab.hstab_span || 0;
  const tTail = (airfoils.tail_t_over_c || 0.1) * Math.max(hChord, 1e-6);
  const hZ = tTail / 2 + 0.002;
  const hstabCells: [number, number, number][] = [];
  budget = nHstab; // hstab gets exactly the optimizer's overflow count
  if (budget > 0 && hChord > 0 && hSpan > 0) {
    const xLE = -0.4 * hChord + 0.1 * hChord;      // LE at -0.4c (mesh convention) + 10% margin
    const nRows = Math.floor((0.6 * hChord) / pitch); // panel-able band is 60% of chord
    const nCols = Math.floor(hSpan / pitch);
    const y0 = -((nCols - 1) * pitch) / 2;          // centered spanwise
    for (let col = 0; col < nCols && budget > 0; col++) {
      const yc = y0 + col * pitch;
      for (let r = 0; r < nRows && budget > 0; r++) {
        hstabCells.push([xLE + (r + 0.5) * pitch, yc, hZ]);
        budget--;
      }
    }
  }

  return { wing, hstab: hstabCells, side };
}

const CATEGORY_COLORS: Record<string, number> = {
  structure: 0x4a90e2,
  power: 0xe74c3c,
  avionics: 0x27ae60,
  propulsion: 0xe67e22,
};

export function colorFor(name: string, category?: string): number {
  if (category && CATEGORY_COLORS[category]) return CATEGORY_COLORS[category];
  if (name.includes("Boom") || name.includes("Fuselage")) return 0x8b7355;
  if (name.includes("Solar")) return 0xf1c40f;
  return 0x9b59b6;
}
